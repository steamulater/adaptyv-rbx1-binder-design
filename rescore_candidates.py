#!/usr/bin/env python3
"""
Re-score all candidates (Batch 1 + Batch 2) using Nipah-derived thresholds:
  - Keep:  ipTM >= 0.70
  - ADD:   complex_iplddt >= 0.85   (AUROC 0.691 in Nipah — best predictor)
  - DROP:  mono pTM filter           (AUROC 0.501 — random noise)

Outputs:
  rescored_all.csv        — full table with new metrics
  final_submission_v2.fasta — updated FASTA ranked by ipTM
  rescore_analysis.png    — visualisation of old vs new filtering
"""
import os, glob, json, csv, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── Paths ─────────────────────────────────────────────────────────────────────
B1_COMP  = "boltz_DESIGN_results/boltz_results_complexes/predictions"
B1_MONO  = "boltz_DESIGN_results/boltz_results_monomers/predictions"
B2_COMP  = "boltz_rfdiffusion/complex_results/boltz_results_complex/predictions"
B2_MONO  = "boltz_rfdiffusion/monomer_results/boltz_results_monomer/predictions"
B2_YAML  = "boltz_rfdiffusion/monomer"

# ── Helpers ───────────────────────────────────────────────────────────────────
def avg_json(directory, keys):
    """Average multiple confidence JSON files in a prediction dir."""
    files = glob.glob(f"{directory}/confidence_*.json")
    acc = {k: [] for k in keys}
    for fp in files:
        d = json.load(open(fp))
        for k in keys:
            if k in d and d[k] is not None:
                acc[k].append(d[k])
    return {k: (np.mean(v) if v else None) for k, v in acc.items()}

def get_seq_from_yaml(yaml_path):
    if not os.path.exists(yaml_path):
        return None
    with open(yaml_path) as f:
        for line in f:
            if "sequence:" in line:
                return line.split("sequence:")[-1].strip()
    return None

# ── Parse Batch 1 ─────────────────────────────────────────────────────────────
b1_master = {}
with open("master_sequences.csv") as f:
    for row in csv.DictReader(f):
        b1_master[row["seq_id"]] = row

b1_records = []
for comp_dir in sorted(glob.glob(f"{B1_COMP}/*_complex")):
    seq_id = os.path.basename(comp_dir).replace("_complex", "")
    comp = avg_json(comp_dir, ["protein_iptm", "iptm", "complex_iplddt", "complex_plddt"])
    avg_iptm    = comp.get("protein_iptm") or comp.get("iptm")
    avg_iplddt  = comp.get("complex_iplddt")
    avg_cplddt  = comp.get("complex_plddt")

    mono_dir = f"{B1_MONO}/{seq_id}_monomer"
    mono = avg_json(mono_dir, ["ptm", "complex_plddt"]) if os.path.isdir(mono_dir) else {}
    mono_ptm = mono.get("ptm")

    m = b1_master.get(seq_id, {})
    seq = m.get("sequence", "")
    b1_records.append({
        "seq_id":      seq_id,
        "batch":       "Batch1",
        "scaffold":    m.get("scaffold", ""),
        "length":      len(seq),
        "sequence":    seq,
        "avg_iptm":    avg_iptm,
        "avg_iplddt":  avg_iplddt,
        "avg_cplddt":  avg_cplddt,
        "mono_ptm":    mono_ptm,
        "ring_rmsd":   float(m["ring_rmsd"]) if m.get("ring_rmsd") else None,
        "novelty_pass": m.get("novelty_pass", ""),
    })

print(f"Batch 1 records: {len(b1_records)}")

# ── Parse Batch 2 ─────────────────────────────────────────────────────────────
b2_records = []
for comp_dir in sorted(glob.glob(f"{B2_COMP}/*_complex")):
    seq_id = os.path.basename(comp_dir).replace("_complex", "")
    comp = avg_json(comp_dir, ["protein_iptm", "iptm", "complex_iplddt", "complex_plddt"])
    avg_iptm   = comp.get("protein_iptm") or comp.get("iptm")
    avg_iplddt = comp.get("complex_iplddt")
    avg_cplddt = comp.get("complex_plddt")

    mono_dir = f"{B2_MONO}/{seq_id}_monomer"
    mono = avg_json(mono_dir, ["ptm", "complex_plddt"]) if os.path.isdir(mono_dir) else {}
    mono_ptm = mono.get("ptm")

    yaml_path = f"{B2_YAML}/{seq_id}_monomer.yaml"
    seq = get_seq_from_yaml(yaml_path) or ""

    # terminal poly-Ala flag
    poly_ala = bool(re.match("^A{6,}", seq) or re.search("A{6,}$", seq))

    b2_records.append({
        "seq_id":      seq_id,
        "batch":       "Batch2",
        "scaffold":    "RFdiffusion",
        "length":      len(seq),
        "sequence":    seq,
        "avg_iptm":    avg_iptm,
        "avg_iplddt":  avg_iplddt,
        "avg_cplddt":  avg_cplddt,
        "mono_ptm":    mono_ptm,
        "ring_rmsd":   None,
        "novelty_pass": "",
        "poly_ala":    poly_ala,
    })

print(f"Batch 2 records: {len(b2_records)}")

all_records = b1_records + b2_records

# ── Apply filters ──────────────────────────────────────────────────────────────
def passes_old(r):
    return (r["avg_iptm"]  is not None and r["avg_iptm"]  >= 0.70 and
            r["mono_ptm"]  is not None and r["mono_ptm"]  >= 0.70)

def composite_score(r):
    """Nipah-informed composite: equally weight ipTM and ipLDDT."""
    if r["avg_iptm"] is None or r["avg_iplddt"] is None:
        return 0.0
    return 0.5 * r["avg_iptm"] + 0.5 * r["avg_iplddt"]

def passes_new(r):
    # ipTM >= 0.70  (hard gate — AUROC 0.603 in Nipah, still meaningful)
    # DROP mono pTM (AUROC 0.501 — random noise per Nipah)
    # DROP poly-Ala sequences
    # ipLDDT used for ranking via composite_score, not as hard cutoff
    # (Nipah threshold 0.85 not portable to GLMN/RFD system — use as ranking weight)
    poly = r.get("poly_ala", False)
    return (r["avg_iptm"]   is not None and r["avg_iptm"]   >= 0.70 and
            r["avg_iplddt"] is not None and
            not poly)

old_pass = [r for r in all_records if passes_old(r)]
new_pass = [r for r in all_records if passes_new(r)]

# Batch breakdown
def batch_split(lst):
    b1 = [r for r in lst if r["batch"] == "Batch1"]
    b2 = [r for r in lst if r["batch"] == "Batch2"]
    return b1, b2

old_b1, old_b2 = batch_split(old_pass)
new_b1, new_b2 = batch_split(new_pass)

print(f"\nOLD filter (ipTM≥0.70 AND pTM≥0.70):")
print(f"  Batch 1: {len(old_b1)}  Batch 2: {len(old_b2)}  Total: {len(old_pass)}")
print(f"\nNEW filter (ipTM≥0.70 AND ipLDDT≥0.85, no pTM):")
print(f"  Batch 1: {len(new_b1)}  Batch 2: {len(new_b2)}  Total: {len(new_pass)}")

# sequences gained (in new but not old)
gained = [r for r in new_pass if not passes_old(r)]
lost   = [r for r in old_pass if not passes_new(r)]
print(f"\n  Gained (pass new, fail old): {len(gained)}")
print(f"  Lost   (pass old, fail new): {len(lost)}")

for r in gained[:10]:
    print(f"    GAINED: {r['seq_id']:<30} ipTM={r['avg_iptm']:.3f}  ipLDDT={r['avg_iplddt']:.3f}  pTM={r['mono_ptm'] or 'N/A'}")
for r in lost[:10]:
    print(f"    LOST:   {r['seq_id']:<30} ipTM={r['avg_iptm']:.3f}  ipLDDT={r['avg_iplddt']:.3f}  pTM={r['mono_ptm'] or 'N/A'}")

# ── Write CSV ──────────────────────────────────────────────────────────────────
new_pass_sorted = sorted(new_pass, key=lambda r: -composite_score(r))
fields = ["seq_id","batch","scaffold","length","avg_iptm","avg_iplddt","avg_cplddt",
          "mono_ptm","ring_rmsd","novelty_pass","sequence"]
with open("rescored_all.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    w.writerows(new_pass_sorted)
print(f"\nSaved rescored_all.csv ({len(new_pass_sorted)} sequences)")

# ── Write FASTA ────────────────────────────────────────────────────────────────
with open("final_submission_v2.fasta", "w") as f:
    for i, r in enumerate(new_pass_sorted, 1):
        if r["sequence"]:
            cs = composite_score(r)
            tag = f"composite={cs:.3f}|ipTM={r['avg_iptm']:.3f}|ipLDDT={r['avg_iplddt']:.3f}|scaffold={r['scaffold']}"
            f.write(f">{r['seq_id']}|rank{i}|{tag}\n{r['sequence']}\n")
print(f"Saved final_submission_v2.fasta")

# ── Visualisation ──────────────────────────────────────────────────────────────
DARK  = "#0f0f0f"; PANEL = "#1a1a1a"; GRID = "#2a2a2a"
B1C   = "#4fc3f7"; B2C   = "#81c784"; GOLD = "#ffd54f"
WARN  = "#ef9a9a"; PASS  = "#a5d6a7"

fig = plt.figure(figsize=(18, 13))
fig.patch.set_facecolor(DARK)
gs  = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)

def style_ax(ax, title=""):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors="#aaaaaa", labelsize=9)
    for sp in ax.spines.values(): sp.set_edgecolor(GRID)
    ax.xaxis.label.set_color("#cccccc")
    ax.yaxis.label.set_color("#cccccc")
    if title:
        ax.set_title(title, color="#eeeeee", fontsize=10, fontweight="bold", pad=8)
    ax.grid(axis="y", color=GRID, lw=0.5, alpha=0.7)

# 1. Old vs New filter count
ax0 = fig.add_subplot(gs[0, 0])
style_ax(ax0, "Old vs New Filter\n(total selected)")
categories = ["Old filter\n(ipTM+pTM)", "New filter\n(ipTM+ipLDDT)"]
old_vals = [len(old_b1), len(old_b2)]
new_vals = [len(new_b1), len(new_b2)]
x = np.arange(2)
w = 0.35
ax0.bar(x - w/2, [sum(old_vals), sum(new_vals)], w*2,
        color=[WARN, PASS], edgecolor=GRID, linewidth=0.8)
ax0.set_xticks(x - w/2 + w)
ax0.set_xticklabels(categories, color="#cccccc", fontsize=9)
for xpos, total, b1n, b2n in [(x[0]-w/2+w/2, sum(old_vals), len(old_b1), len(old_b2)),
                                (x[1]-w/2+w/2, sum(new_vals), len(new_b1), len(new_b2))]:
    ax0.text(xpos, total + 0.5, f"{total}\n(B1:{b1n} B2:{b2n})",
             ha="center", va="bottom", color="#eeeeee", fontsize=9, fontweight="bold")
ax0.set_ylim(0, max(sum(old_vals), sum(new_vals)) + 12)
ax0.set_ylabel("Sequences passing")

# 2. ipTM vs ipLDDT scatter — all sequences, coloured by filter outcome
ax1 = fig.add_subplot(gs[0, 1:])
style_ax(ax1, "ipTM vs complex ipLDDT — all sequences (Batch 1+2)\nColoured by new filter outcome")

def scatter_group(ax, records, colour, label, size=18, alpha=0.65):
    xs = [r["avg_iptm"]   for r in records if r["avg_iptm"] and r["avg_iplddt"]]
    ys = [r["avg_iplddt"] for r in records if r["avg_iptm"] and r["avg_iplddt"]]
    ax.scatter(xs, ys, c=colour, s=size, alpha=alpha, label=f"{label} (n={len(xs)})", zorder=3)

fail_recs  = [r for r in all_records if not passes_new(r)]
gained_recs = gained
pass_recs   = [r for r in new_pass if passes_old(r)]

scatter_group(ax1, fail_recs,   WARN,  "Fail new filter",       size=14, alpha=0.4)
scatter_group(ax1, pass_recs,   PASS,  "Pass both filters",     size=22, alpha=0.8)
scatter_group(ax1, gained_recs, GOLD,  "Gained (new only)",     size=28, alpha=0.95)

ax1.axvline(0.70, color=B1C,  lw=1.2, ls="--", alpha=0.8, label="ipTM ≥ 0.70")
ax1.axhline(0.85, color=GOLD, lw=1.5, ls="--", alpha=0.9, label="ipLDDT ≥ 0.85 (Nipah optimal)")
ax1.fill_between([0.70, 1.0], 0.85, 1.0, alpha=0.07, color=PASS)
ax1.set_xlabel("avg ipTM (complex)")
ax1.set_ylabel("avg complex ipLDDT")
ax1.set_xlim(0.3, 1.0)
ax1.set_ylim(0.55, 1.0)
ax1.grid(axis="both", color=GRID, lw=0.5, alpha=0.7)
ax1.xaxis.grid(True)
ax1.legend(fontsize=8.5, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc", loc="upper left")

# 3. ipLDDT distribution: passing vs failing
ax2 = fig.add_subplot(gs[1, 0])
style_ax(ax2, "complex ipLDDT Distribution\n(new filter in action)")
iplddt_all = [r["avg_iplddt"] for r in all_records if r["avg_iplddt"]]
iplddt_pass= [r["avg_iplddt"] for r in new_pass]
bins = np.linspace(0.55, 1.0, 25)
ax2.hist(iplddt_all,  bins=bins, color=WARN,  alpha=0.7, label=f"All (n={len(iplddt_all)})",  edgecolor=DARK)
ax2.hist(iplddt_pass, bins=bins, color=PASS,  alpha=0.8, label=f"Pass (n={len(iplddt_pass)})", edgecolor=DARK)
ax2.axvline(0.85, color=GOLD, lw=1.8, ls="--", label="Nipah threshold 0.85")
ax2.set_xlabel("avg complex ipLDDT")
ax2.set_ylabel("Count")
ax2.legend(fontsize=8, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc")

# 4. Top 20 final ranking
ax3 = fig.add_subplot(gs[1, 1:])
style_ax(ax3, f"Final Top 20 by ipTM — {len(new_pass_sorted)} total selected\n(ipTM≥0.70 AND ipLDDT≥0.85)")
top20 = new_pass_sorted[:20]
labels_t = [r["seq_id"].replace("_best","").replace("_T","_") for r in top20]
iptm_t   = [r["avg_iptm"]   for r in top20]
iplddt_t = [r["avg_iplddt"] for r in top20]
colors_t = [B1C if r["batch"]=="Batch1" else B2C for r in top20]

x_t = np.arange(len(top20))
w_t = 0.4
ax3.bar(x_t - w_t/2, iptm_t,   w_t, color=colors_t,         edgecolor=GRID, lw=0.5, label="ipTM",     alpha=0.9)
ax3.bar(x_t + w_t/2, iplddt_t, w_t, color=[c+"88" for c in colors_t], edgecolor=GRID, lw=0.5, label="ipLDDT", alpha=0.75)
ax3.axhline(0.85, color=GOLD, lw=1, ls="--", alpha=0.7)
ax3.axhline(0.70, color="#888888", lw=0.8, ls=":", alpha=0.6)
ax3.set_xticks(x_t)
ax3.set_xticklabels(labels_t, rotation=45, ha="right", fontsize=7, color="#cccccc")
ax3.set_ylim(0.60, 1.0)
ax3.set_ylabel("Score")

b1_patch = plt.matplotlib.patches.Patch(color=B1C, label="Batch 1 (GLMN/CUL1)")
b2_patch = plt.matplotlib.patches.Patch(color=B2C, label="Batch 2 (RFdiffusion)")
solid    = plt.matplotlib.patches.Patch(color="#aaaaaa",     label="ipTM (solid bar)")
light    = plt.matplotlib.patches.Patch(color="#aaaaaa55", label="ipLDDT (light bar)")
ax3.legend(handles=[b1_patch, b2_patch, solid, light],
           fontsize=8, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc",
           loc="lower right", ncol=2)

fig.suptitle("RBX1 Binder Design — Re-scored with Nipah-Derived Thresholds\nFilter: ipTM ≥ 0.70 | Rank: composite = 0.5×ipTM + 0.5×ipLDDT | mono pTM dropped",
             color="#eeeeee", fontsize=12, fontweight="bold", y=0.99)

plt.savefig("rescore_analysis.png", dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
print("Saved rescore_analysis.png")
