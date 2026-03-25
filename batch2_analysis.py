#!/usr/bin/env python3
"""
Batch 2 funnel analysis + Batch 1 vs Batch 2 comparison.
"""
import csv, re, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# ── Load data ─────────────────────────────────────────────────────────────────

# Batch 2
with open("rfd_batch2_results.csv") as f:
    b2_rows = list(csv.DictReader(f))

b2_iptm  = [float(r["avg_iptm"])  for r in b2_rows if r["avg_iptm"]]
b2_ptm   = [float(r["mono_ptm"])  for r in b2_rows if r["mono_ptm"]]
b2_pass  = [r for r in b2_rows
            if r["avg_iptm"] and r["mono_ptm"]
            and float(r["avg_iptm"]) >= 0.70
            and float(r["mono_ptm"]) >= 0.70]

# Terminal poly-Ala detection
terminal_ala = set()
for yaml in glob.glob("boltz_rfdiffusion/monomer/*.yaml"):
    with open(yaml) as f:
        for line in f:
            if "sequence:" in line:
                seq  = line.split("sequence:")[-1].strip()
                name = yaml.split("/")[-1].replace("_monomer.yaml", "")
                if re.match("^A{6,}", seq) or re.search("A{6,}$", seq):
                    terminal_ala.add(name)
                break

b2_clean = [r for r in b2_pass if r["seq_id"] not in terminal_ala]

# Batch 1
with open("master_sequences.csv") as f:
    b1_rows = list(csv.DictReader(f))

glmn_rows = [r for r in b1_rows if r["scaffold"] == "GLMN"]
whb_rows  = [r for r in b1_rows if r["scaffold"] == "CUL1_WHB"]
b1_submitted = [r for r in b1_rows
                if r["boltz_complex_iptm"]
                and r["boltz_monomer_ptm"]
                and float(r["boltz_complex_iptm"]) >= 0.70
                and float(r["boltz_monomer_ptm"])  >= 0.70]

glmn_iptm = [float(r["boltz_complex_iptm"]) for r in glmn_rows if r["boltz_complex_iptm"]]
whb_iptm  = [float(r["boltz_complex_iptm"]) for r in whb_rows  if r["boltz_complex_iptm"]]
glmn_ptm  = [float(r["boltz_monomer_ptm"])  for r in glmn_rows if r["boltz_monomer_ptm"]]
whb_ptm   = [float(r["boltz_monomer_ptm"])  for r in whb_rows  if r["boltz_monomer_ptm"]]

# ── Funnel numbers ────────────────────────────────────────────────────────────
funnel = [
    ("RFdiffusion\nbackbones", 200),
    ("Length filter\n(65-95 AA)", 151),
    ("ProteinMPNN\nbest-per-backbone", 151),
    ("Boltz-2 monomer\ncompleted", 118),
    ("Pass ipTM ≥ 0.70\n& pTM ≥ 0.70", 32),
    ("Clean\n(no poly-Ala)", 31),
]

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 14))
fig.patch.set_facecolor("#0f0f0f")
gs = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

DARK  = "#0f0f0f"
PANEL = "#1a1a1a"
GRID  = "#2a2a2a"
B1C   = "#4fc3f7"   # Batch 1 blue
B2C   = "#81c784"   # Batch 2 green
WARN  = "#ef9a9a"   # warning red
GOLD  = "#ffd54f"

def style_ax(ax, title=""):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors="#aaaaaa", labelsize=9)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID)
    ax.xaxis.label.set_color("#cccccc")
    ax.yaxis.label.set_color("#cccccc")
    if title:
        ax.set_title(title, color="#eeeeee", fontsize=10, fontweight="bold", pad=8)
    ax.grid(axis="y", color=GRID, linewidth=0.5, alpha=0.7)

# ── Panel 1: Funnel (spans top row, col 0-1) ─────────────────────────────────
ax_funnel = fig.add_subplot(gs[0, :2])
style_ax(ax_funnel, "Batch 2 Selection Funnel: 200 → 31")

labels = [f[0] for f in funnel]
values = [f[1] for f in funnel]
colors_f = [B2C if v > 50 else GOLD if v > 30 else WARN for v in values]
colors_f[-1] = B2C

bars = ax_funnel.bar(range(len(funnel)), values, color=colors_f, edgecolor=GRID, linewidth=0.8, width=0.6)
ax_funnel.set_xticks(range(len(funnel)))
ax_funnel.set_xticklabels(labels, fontsize=8.5, color="#cccccc")
ax_funnel.set_ylabel("Sequences", color="#cccccc")
ax_funnel.set_ylim(0, 230)

# Annotate bars
for i, (bar, val) in enumerate(zip(bars, values)):
    ax_funnel.text(bar.get_x() + bar.get_width()/2, val + 4,
                   str(val), ha="center", va="bottom", color="#eeeeee", fontsize=10, fontweight="bold")
    if i > 0:
        drop = values[i-1] - val
        if drop > 0:
            ax_funnel.text(bar.get_x() + bar.get_width()/2, val/2,
                           f"−{drop}", ha="center", va="center", color="#ff7043", fontsize=8, fontweight="bold")

# Annotate key reasons
reasons = [
    (1, "Binder length\nout of range"),
    (3, "Monomer run\nstopped early"),
    (4, "Failed dual\nconfidence filter"),
    (5, "Poly-Ala\ntermini"),
]
for idx, reason in reasons:
    ax_funnel.annotate(reason,
        xy=(idx, values[idx] + 15), xytext=(idx, 190),
        arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8),
        ha="center", fontsize=7.5, color="#aaaaaa",
        bbox=dict(boxstyle="round,pad=0.3", fc="#1a1a1a", ec=GRID, lw=0.8))

# ── Panel 2: Why 32 pass — scatter ipTM vs pTM ───────────────────────────────
ax_sc = fig.add_subplot(gs[0, 2])
style_ax(ax_sc, "Batch 2: ipTM vs mono pTM (n=118)")

pass_mask = np.array([float(r["avg_iptm"]) >= 0.70 and float(r["mono_ptm"]) >= 0.70
                      for r in b2_rows if r["avg_iptm"] and r["mono_ptm"]])
x_all = np.array(b2_iptm)
y_all = np.array(b2_ptm)

ax_sc.scatter(x_all[~pass_mask], y_all[~pass_mask], c=WARN, s=22, alpha=0.7, label="Fail (n=86)", zorder=3)
ax_sc.scatter(x_all[pass_mask],  y_all[pass_mask],  c=B2C,  s=30, alpha=0.9, label="Pass (n=32)", zorder=4)

ax_sc.axvline(0.70, color=GOLD, lw=1.2, ls="--", alpha=0.8)
ax_sc.axhline(0.70, color=GOLD, lw=1.2, ls="--", alpha=0.8)
ax_sc.set_xlabel("avg ipTM (complex)", fontsize=9)
ax_sc.set_ylabel("mono pTM", fontsize=9)
ax_sc.set_xlim(0.3, 0.95)
ax_sc.set_ylim(0.3, 1.0)
ax_sc.text(0.72, 0.33, "ipTM cutoff →", fontsize=7, color=GOLD)
ax_sc.text(0.32, 0.71, "pTM cutoff", fontsize=7, color=GOLD)
ax_sc.legend(fontsize=8, loc="upper left",
             facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc")
# Shade passing quadrant
ax_sc.fill_between([0.70, 0.95], 0.70, 1.0, alpha=0.07, color=B2C)
ax_sc.grid(axis="both", color=GRID, linewidth=0.5, alpha=0.7)
ax_sc.xaxis.grid(True)

# ── Panel 3: ipTM distribution comparison ────────────────────────────────────
ax_iptm = fig.add_subplot(gs[1, :2])
style_ax(ax_iptm, "ipTM Distribution: Batch 1 vs Batch 2")

bins = np.linspace(0.3, 1.0, 22)
ax_iptm.hist(glmn_iptm, bins=bins, alpha=0.75, color=B1C,  label="Batch 1 — GLMN (n=48)", edgecolor=DARK)
ax_iptm.hist(whb_iptm,  bins=bins, alpha=0.75, color="#ce93d8", label="Batch 1 — CUL1_WHB (n=48)", edgecolor=DARK)
ax_iptm.hist(b2_iptm,   bins=bins, alpha=0.75, color=B2C,  label="Batch 2 — RFdiffusion (n=118)", edgecolor=DARK)
ax_iptm.axvline(0.70, color=GOLD, lw=1.5, ls="--", alpha=0.9, label="Cutoff 0.70")

ax_iptm.set_xlabel("avg ipTM", fontsize=9)
ax_iptm.set_ylabel("Count", fontsize=9)
ax_iptm.legend(fontsize=8, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc")

# Means
for vals, col, lbl in [(glmn_iptm, B1C, f"GLMN μ={np.mean(glmn_iptm):.3f}"),
                        (b2_iptm,   B2C, f"RFD μ={np.mean(b2_iptm):.3f}")]:
    ax_iptm.axvline(np.mean(vals), color=col, lw=1.2, ls=":", alpha=0.9)

# ── Panel 4: mono pTM distribution comparison ────────────────────────────────
ax_ptm = fig.add_subplot(gs[1, 2])
style_ax(ax_ptm, "Monomer pTM: Batch 1 vs Batch 2")

ax_ptm.hist(glmn_ptm, bins=bins, alpha=0.75, color=B1C, label="GLMN (n=48)", edgecolor=DARK)
ax_ptm.hist(b2_ptm,   bins=bins, alpha=0.75, color=B2C, label="RFD (n=118)", edgecolor=DARK)
ax_ptm.axvline(0.70, color=GOLD, lw=1.5, ls="--", alpha=0.9, label="Cutoff 0.70")
ax_ptm.set_xlabel("mono pTM", fontsize=9)
ax_ptm.set_ylabel("Count", fontsize=9)
ax_ptm.legend(fontsize=8, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc")

# ── Panel 5: Pass rate comparison ────────────────────────────────────────────
ax_rate = fig.add_subplot(gs[2, 0])
style_ax(ax_rate, "Pass Rate by Stage")

categories = ["GLMN\n(B1)", "CUL1_WHB\n(B1)", "RFdiffusion\n(B2)"]
iptm_rates  = [
    sum(1 for v in glmn_iptm if v >= 0.70) / len(glmn_iptm) * 100,
    sum(1 for v in whb_iptm  if v >= 0.70) / len(whb_iptm)  * 100,
    sum(1 for v in b2_iptm   if v >= 0.70) / len(b2_iptm)   * 100,
]
ptm_rates = [
    sum(1 for v in glmn_ptm if v >= 0.70) / len(glmn_ptm) * 100,
    sum(1 for v in whb_ptm  if v >= 0.70) / len(whb_ptm)  * 100,
    sum(1 for v in b2_ptm   if v >= 0.70) / len(b2_ptm)   * 100,
]

x = np.arange(len(categories))
w = 0.35
bars1 = ax_rate.bar(x - w/2, iptm_rates, w, label="ipTM ≥ 0.70", color=B1C, edgecolor=DARK)
bars2 = ax_rate.bar(x + w/2, ptm_rates,  w, label="pTM ≥ 0.70",  color=B2C, edgecolor=DARK)
ax_rate.set_xticks(x)
ax_rate.set_xticklabels(categories, fontsize=8.5, color="#cccccc")
ax_rate.set_ylabel("Pass rate (%)", color="#cccccc")
ax_rate.set_ylim(0, 115)
ax_rate.legend(fontsize=8, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc")
for bar in list(bars1) + list(bars2):
    ax_rate.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                 f"{bar.get_height():.0f}%", ha="center", va="bottom", color="#eeeeee", fontsize=8)

# ── Panel 6: Mean scores comparison table ────────────────────────────────────
ax_tbl = fig.add_subplot(gs[2, 1])
style_ax(ax_tbl, "Mean Scores Comparison")
ax_tbl.axis("off")

b2_pass_iptm = [float(r["avg_iptm"]) for r in b2_pass]
b2_pass_ptm  = [float(r["mono_ptm"]) for r in b2_pass]

table_data = [
    ["Metric",          "GLMN (B1)",  "CUL1_WHB (B1)", "RFD B2\n(all)", "RFD B2\n(pass)"],
    ["n submitted",     "48",         "7",              "—",             "31"],
    ["mean ipTM",       f"{np.mean(glmn_iptm):.3f}", f"{np.mean(whb_iptm):.3f}",
                        f"{np.mean(b2_iptm):.3f}",  f"{np.mean(b2_pass_iptm):.3f}"],
    ["mean mono pTM",   f"{np.mean(glmn_ptm):.3f}", f"{np.mean(whb_ptm):.3f}",
                        f"{np.mean(b2_ptm):.3f}",  f"{np.mean(b2_pass_ptm):.3f}"],
    ["pass rate",       "100%",        "15%",           f"{len(b2_pass)/len(b2_iptm)*100:.0f}%", "100%"],
    ["length (AA)",     "247",         "72",             "65-95",         "65-95"],
    ["scaffold type",   "GLMN PDB",   "CUL1 PDB",       "De novo",       "De novo"],
]

tbl = ax_tbl.table(cellText=table_data[1:], colLabels=table_data[0],
                   cellLoc="center", loc="center",
                   bbox=[0, 0, 1, 1])
tbl.auto_set_font_size(False)
tbl.set_fontsize(8)
for (row, col), cell in tbl.get_celld().items():
    cell.set_facecolor(PANEL if row % 2 == 0 else "#222222")
    cell.set_edgecolor(GRID)
    cell.set_text_props(color="#eeeeee" if row > 0 else "#ffd54f")
    if row == 0:
        cell.set_facecolor("#1e3a4a")

# ── Panel 7: Top candidates bar ───────────────────────────────────────────────
ax_top = fig.add_subplot(gs[2, 2])
style_ax(ax_top, "Top 15 Batch 2 by ipTM")

b2_pass_sorted = sorted(b2_pass, key=lambda r: -float(r["avg_iptm"]))[:15]
labels_top = [r["seq_id"].replace("_best", "") for r in b2_pass_sorted]
vals_top   = [float(r["avg_iptm"]) for r in b2_pass_sorted]
colors_top = [WARN if r["seq_id"] in terminal_ala else B2C for r in b2_pass_sorted]

ax_top.barh(range(len(labels_top)), vals_top, color=colors_top, edgecolor=GRID, linewidth=0.5)
ax_top.set_yticks(range(len(labels_top)))
ax_top.set_yticklabels(labels_top, fontsize=7.5, color="#cccccc")
ax_top.invert_yaxis()
ax_top.axvline(0.70, color=GOLD, lw=1, ls="--", alpha=0.8)
ax_top.set_xlabel("avg ipTM", fontsize=9)
ax_top.set_xlim(0.60, 0.95)
ax_top.grid(axis="x", color=GRID, linewidth=0.5, alpha=0.7)
ax_top.yaxis.grid(False)

# GLMN mean reference
ax_top.axvline(np.mean(glmn_iptm), color=B1C, lw=1.2, ls=":", alpha=0.8,
               label=f"GLMN mean {np.mean(glmn_iptm):.3f}")
ax_top.legend(fontsize=7, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc", loc="lower right")

# ── Title ─────────────────────────────────────────────────────────────────────
fig.suptitle("RBX1 Binder Design — Batch 2 Funnel Analysis & Batch 1 vs Batch 2 Comparison",
             color="#eeeeee", fontsize=13, fontweight="bold", y=0.98)

plt.savefig("batch2_analysis.png", dpi=150, bbox_inches="tight",
            facecolor=fig.get_facecolor())
print("Saved batch2_analysis.png")
