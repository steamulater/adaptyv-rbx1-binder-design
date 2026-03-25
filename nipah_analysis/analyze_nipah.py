#!/usr/bin/env python3
"""
Nipah binder competition retrospective analysis.
Parses all computational + experimental metrics, correlates predictors with
experimental binding outcome, generates summary + visualisation.
"""
import csv, json, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import defaultdict
from scipy import stats

# ── Parse ─────────────────────────────────────────────────────────────────────
DATA = os.path.join(os.path.dirname(__file__),
                    "../old_proteinbase_collection_nipah-binder-competition-results.csv")

def parse_row(row):
    evals = json.loads(row["evaluations"])
    r = {
        "id":     row["id"],
        "name":   row["name"],
        "author": row["author"],
        "method": row["designMethod"],
        "seq":    row["sequence"],
        "length": len(row["sequence"]),
    }
    for e in evals:
        m = e["metric"]
        t = e.get("target", "")
        v = e.get("value")
        vt = e.get("valueType", "")
        key = f"{m}:{t}" if t else m

        # skip JSON blob values we can't use as numerics/labels
        if vt == "json":
            continue

        # last-write wins (handles repeated measurements — we'll aggregate separately)
        if vt == "numeric" and isinstance(v, (int, float)):
            r[key] = v
        elif vt in ("boolean", "label") and not isinstance(v, dict):
            r[key] = v

    return r

with open(DATA, newline="", encoding="utf-8-sig") as f:
    rows = [parse_row(row) for row in csv.DictReader(f)]

# ── Experimental outcome ───────────────────────────────────────────────────────
# binding:nipah-glycoprotein-g can appear multiple times — use majority vote
def get_binding_vote(row_raw):
    evals = json.loads(row_raw["evaluations"])
    votes = [e["value"] for e in evals
             if e["metric"] == "binding" and e.get("target") == "nipah-glycoprotein-g"
             and isinstance(e["value"], bool)]
    if not votes:
        return None
    return sum(votes) > len(votes) / 2   # majority True = binder

def get_min_kd(row_raw):
    evals = json.loads(row_raw["evaluations"])
    kds = [e["value"] for e in evals
           if e["metric"] == "kd" and e.get("target") == "nipah-glycoprotein-g"
           and isinstance(e["value"], (int, float))]
    return min(kds) if kds else None

def get_expressed(row_raw):
    evals = json.loads(row_raw["evaluations"])
    votes = [e["value"] for e in evals
             if e["metric"] == "expressed" and isinstance(e["value"], bool)]
    if not votes:
        return None
    return sum(votes) > len(votes) / 2

with open(DATA, newline="", encoding="utf-8-sig") as f:
    raw_rows = list(csv.DictReader(f))

for i, rr in enumerate(raw_rows):
    rows[i]["binding"]   = get_binding_vote(rr)
    rows[i]["min_kd"]    = get_min_kd(rr)
    rows[i]["expressed"] = get_expressed(rr)

# ── Summary numbers ────────────────────────────────────────────────────────────
total       = len(rows)
has_exp     = [r for r in rows if r["binding"] is not None]
binders     = [r for r in has_exp if r["binding"]]
non_binders = [r for r in has_exp if not r["binding"]]
expressed   = [r for r in rows if r.get("expressed") is True]
no_express  = [r for r in rows if r.get("expressed") is False]
kd_binders  = [r["min_kd"] for r in binders if r["min_kd"]]

print(f"Total sequences:        {total}")
print(f"  Experimentally tested:{len(has_exp)}")
print(f"  Binders:              {len(binders)} ({100*len(binders)/len(has_exp):.1f}%)")
print(f"  Non-binders:          {len(non_binders)}")
print(f"  Expressed:            {len(expressed)}")
print(f"  Not expressed:        {len(no_express)}")
print(f"  KD data:              {len(kd_binders)}")
print(f"  Median KD (M):        {np.median(kd_binders):.2e}" if kd_binders else "  No KD data")
print()

# design class breakdown
dc = defaultdict(lambda: {"total":0,"bind":0})
for r in has_exp:
    cls = r.get("design_class","unknown")
    dc[cls]["total"] += 1
    if r["binding"]: dc[cls]["bind"] += 1
print("Binding by design class:")
for cls, d in sorted(dc.items(), key=lambda x: -x[1]["total"]):
    rate = d["bind"]/d["total"]*100 if d["total"] else 0
    print(f"  {cls:<20} {d['bind']:>3}/{d['total']:>3} ({rate:.0f}%)")

# ── Computational predictors vs outcome ───────────────────────────────────────
COMP_METRICS = [
    ("boltz2_iptm:nipah-glycoprotein-g",           "Boltz2 ipTM"),
    ("boltz2_ptm:nipah-glycoprotein-g",             "Boltz2 pTM"),
    ("boltz2_plddt:nipah-glycoprotein-g",           "Boltz2 pLDDT"),
    ("boltz2_complex_plddt:nipah-glycoprotein-g",   "Boltz2 complex pLDDT"),
    ("boltz2_complex_iplddt:nipah-glycoprotein-g",  "Boltz2 complex ipLDDT"),
    ("boltz2_pdockq:nipah-glycoprotein-g",          "Boltz2 pDockQ"),
    ("boltz2_pdockq2:nipah-glycoprotein-g",         "Boltz2 pDockQ2"),
    ("boltz2_lis:nipah-glycoprotein-g",             "Boltz2 LIS"),
    ("boltz2_ipsae:nipah-glycoprotein-g",           "Boltz2 ipSAE"),
    ("boltz2_min_ipsae:nipah-glycoprotein-g",       "Boltz2 min ipSAE"),
    ("shape_complimentarity_boltz2_binder_ss:nipah-glycoprotein-g", "Shape complementarity"),
    ("esmfold_plddt",                               "ESMFold pLDDT"),
    ("proteinmpnn_score",                           "ProteinMPNN score"),
    ("redesigned_proteinmpnn_score",                "Redesigned MPNN score"),
    ("isoelectric_point",                           "Isoelectric point"),
    ("molecular_weight",                            "Molecular weight"),
]

print("\nPredictor analysis (binders vs non-binders):")
print(f"{'Metric':<40} {'Binders':>10} {'Non-bind':>10} {'p-val':>10} {'AUROC':>8}")
print("-"*82)

from sklearn.metrics import roc_auc_score

predictor_stats = []
for key, label in COMP_METRICS:
    b_vals  = [r[key] for r in binders     if key in r and r[key] is not None]
    nb_vals = [r[key] for r in non_binders if key in r and r[key] is not None]
    if len(b_vals) < 5 or len(nb_vals) < 5:
        continue

    _, pval = stats.mannwhitneyu(b_vals, nb_vals, alternative="two-sided")

    # AUROC: higher value = binder
    all_vals = b_vals + nb_vals
    labels   = [1]*len(b_vals) + [0]*len(nb_vals)
    try:
        auroc = roc_auc_score(labels, all_vals)
        if auroc < 0.5:
            auroc = 1 - auroc   # flip direction
    except:
        auroc = float("nan")

    predictor_stats.append((label, key, np.mean(b_vals), np.mean(nb_vals), pval, auroc,
                            b_vals, nb_vals))
    star = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
    print(f"{label:<40} {np.mean(b_vals):>10.3f} {np.mean(nb_vals):>10.3f} "
          f"{pval:>10.2e} {auroc:>8.3f} {star}")

predictor_stats.sort(key=lambda x: -x[5])   # sort by AUROC

# ── KD analysis ───────────────────────────────────────────────────────────────
print(f"\nKD distribution (binders with KD data, n={len(kd_binders)}):")
for thresh, label in [(1e-9,"<1 nM"),(5e-9,"1-5 nM"),(10e-9,"5-10 nM"),(50e-9,"10-50 nM"),(1,"≥50 nM")]:
    n = sum(1 for v in kd_binders if v < thresh)
    print(f"  {label}: {n}")

# ── Optimal thresholds ─────────────────────────────────────────────────────────
print("\nOptimal thresholds (Youden index) for top predictors:")
for label, key, bm, nbm, pval, auroc, b_vals, nb_vals in predictor_stats[:6]:
    all_v  = sorted(set(b_vals + nb_vals))
    best_j, best_thresh = 0, None
    for t in all_v:
        tp = sum(1 for v in b_vals  if v >= t) / len(b_vals)   # sensitivity
        tn = sum(1 for v in nb_vals if v <  t) / len(nb_vals)  # specificity
        j  = tp + tn - 1
        if j > best_j:
            best_j, best_thresh = j, t
    if best_thresh is not None:
        print(f"  {label:<40} threshold={best_thresh:.3f}  J={best_j:.3f}")

# ── VISUALISATION ─────────────────────────────────────────────────────────────
DARK  = "#0f0f0f"
PANEL = "#1a1a1a"
GRID  = "#2a2a2a"
BIND  = "#81c784"
NOBIND= "#ef9a9a"
GOLD  = "#ffd54f"
BLUE  = "#4fc3f7"

fig = plt.figure(figsize=(20, 16))
fig.patch.set_facecolor(DARK)
gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=0.50, wspace=0.38)

def style_ax(ax, title=""):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors="#aaaaaa", labelsize=8.5)
    for sp in ax.spines.values(): sp.set_edgecolor(GRID)
    ax.xaxis.label.set_color("#cccccc")
    ax.yaxis.label.set_color("#cccccc")
    if title:
        ax.set_title(title, color="#eeeeee", fontsize=9.5, fontweight="bold", pad=7)
    ax.grid(axis="y", color=GRID, lw=0.5, alpha=0.7)

# 1. Outcome pie
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_facecolor(PANEL)
ax0.set_title("Experimental Outcomes\n(n=1030)", color="#eeeeee", fontsize=9.5, fontweight="bold", pad=7)
pie_vals = [len(binders), len(non_binders), total - len(has_exp)]
pie_lbls = [f"Binders\n{len(binders)}", f"Non-binders\n{len(non_binders)}",
            f"Not tested\n{total-len(has_exp)}"]
colors_p = [BIND, NOBIND, "#555555"]
wedges, texts, autotexts = ax0.pie(pie_vals, labels=pie_lbls, colors=colors_p,
                                    autopct="%1.0f%%", startangle=90,
                                    textprops={"color":"#cccccc","fontsize":8},
                                    pctdistance=0.7)
for at in autotexts: at.set_color("#111111"); at.set_fontsize(8)

# 2. Binding by design class
ax1 = fig.add_subplot(gs[0, 1])
style_ax(ax1, "Binding Rate by Design Class")
classes = sorted(dc.keys(), key=lambda x: -dc[x]["bind"]/max(dc[x]["total"],1))[:7]
rates  = [dc[c]["bind"]/dc[c]["total"]*100 for c in classes]
totals = [dc[c]["total"] for c in classes]
bars = ax1.bar(range(len(classes)), rates,
               color=[BIND if r >= 50 else NOBIND for r in rates],
               edgecolor=GRID, linewidth=0.7)
ax1.set_xticks(range(len(classes)))
ax1.set_xticklabels([c.replace("_","\n") for c in classes], fontsize=7, color="#cccccc")
ax1.set_ylabel("Binding rate (%)")
ax1.set_ylim(0, 110)
for i, (bar, r, n) in enumerate(zip(bars, rates, totals)):
    ax1.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1.5,
             f"{r:.0f}%\n(n={n})", ha="center", va="bottom", color="#eeeeee", fontsize=7)

# 3. KD distribution
ax2 = fig.add_subplot(gs[0, 2])
style_ax(ax2, f"KD Distribution (binders, n={len(kd_binders)})")
kd_nm = [v*1e9 for v in kd_binders]
ax2.hist(kd_nm, bins=30, color=BIND, edgecolor=DARK, alpha=0.85)
ax2.axvline(np.median(kd_nm), color=GOLD, lw=1.5, ls="--",
            label=f"Median {np.median(kd_nm):.1f} nM")
ax2.set_xlabel("KD (nM)")
ax2.set_ylabel("Count")
ax2.set_xlim(0, 100)
ax2.legend(fontsize=8, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc")

# 4. AUROC bar — top predictors
ax3 = fig.add_subplot(gs[0, 3])
style_ax(ax3, "Predictor AUROC\n(binders vs non-binders)")
top = predictor_stats[:8]
labels_a = [x[0].replace("Boltz2 ","B2 ").replace("ProteinMPNN","MPNN") for x in top]
aurocs    = [x[5] for x in top]
colors_a  = [BIND if a >= 0.70 else GOLD if a >= 0.60 else NOBIND for a in aurocs]
ax3.barh(range(len(top)), aurocs, color=colors_a, edgecolor=GRID, linewidth=0.6)
ax3.set_yticks(range(len(top)))
ax3.set_yticklabels(labels_a, fontsize=8, color="#cccccc")
ax3.invert_yaxis()
ax3.axvline(0.70, color=GOLD, lw=1.2, ls="--", alpha=0.8, label="AUROC 0.70")
ax3.axvline(0.50, color="#888888", lw=0.8, ls=":", alpha=0.6)
ax3.set_xlabel("AUROC")
ax3.set_xlim(0.4, 1.0)
ax3.yaxis.grid(False)
ax3.legend(fontsize=7.5, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc")
for i, (bar, a) in enumerate(zip(ax3.patches, aurocs)):
    ax3.text(a + 0.005, i, f"{a:.3f}", va="center", color="#eeeeee", fontsize=7.5)

# 5–8. Top 4 predictors — violin/box binders vs non
top4 = predictor_stats[:4]
axes_v = [fig.add_subplot(gs[1, i]) for i in range(4)]
for ax, (label, key, bm, nbm, pval, auroc, b_vals, nb_vals) in zip(axes_v, top4):
    style_ax(ax, label.replace("Boltz2 ","Boltz2\n"))
    vp = ax.violinplot([nb_vals, b_vals], positions=[0,1],
                       showmedians=True, showextrema=False)
    vp["cmedians"].set_color(GOLD)
    vp["cmedians"].set_linewidth(2)
    for i, body in enumerate(vp["bodies"]):
        body.set_facecolor(NOBIND if i==0 else BIND)
        body.set_alpha(0.75)
    ax.set_xticks([0,1])
    ax.set_xticklabels(["Non-binder","Binder"], fontsize=8, color="#cccccc")
    ax.set_ylabel("Score")
    star = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
    ax.text(0.5, 0.97, f"p={pval:.1e} {star}  AUROC={auroc:.3f}",
            ha="center", va="top", transform=ax.transAxes,
            color=GOLD, fontsize=7.5)

# 9. ipTM vs LIS coloured by binding
ax_sc1 = fig.add_subplot(gs[2, :2])
style_ax(ax_sc1, "Boltz2 ipTM vs LIS (coloured by experimental binding)")
iptm_b  = [r.get("boltz2_iptm:nipah-glycoprotein-g") for r in binders  if "boltz2_iptm:nipah-glycoprotein-g" in r]
lis_b   = [r.get("boltz2_lis:nipah-glycoprotein-g")  for r in binders  if "boltz2_lis:nipah-glycoprotein-g"  in r and "boltz2_iptm:nipah-glycoprotein-g" in r]
iptm_nb = [r.get("boltz2_iptm:nipah-glycoprotein-g") for r in non_binders if "boltz2_iptm:nipah-glycoprotein-g" in r]
lis_nb  = [r.get("boltz2_lis:nipah-glycoprotein-g")  for r in non_binders if "boltz2_lis:nipah-glycoprotein-g"  in r and "boltz2_iptm:nipah-glycoprotein-g" in r]

# make sure same length
min_b  = min(len(iptm_b), len(lis_b))
min_nb = min(len(iptm_nb), len(lis_nb))
ax_sc1.scatter(iptm_nb[:min_nb], lis_nb[:min_nb], c=NOBIND, s=12, alpha=0.5, label=f"Non-binder (n={min_nb})", zorder=3)
ax_sc1.scatter(iptm_b[:min_b],   lis_b[:min_b],   c=BIND,   s=15, alpha=0.7, label=f"Binder (n={min_b})",    zorder=4)
ax_sc1.set_xlabel("Boltz2 ipTM")
ax_sc1.set_ylabel("Boltz2 LIS")
ax_sc1.legend(fontsize=8.5, facecolor=PANEL, edgecolor=GRID, labelcolor="#cccccc")
ax_sc1.grid(axis="both", color=GRID, lw=0.5, alpha=0.7)
ax_sc1.xaxis.grid(True)

# 10. Recommended thresholds table
ax_tbl = fig.add_subplot(gs[2, 2:])
style_ax(ax_tbl, "Recommended Thresholds for RBX1 Submission\n(derived from Nipah results)")
ax_tbl.axis("off")

# compute Youden thresholds for top metrics
thresholds = {}
for label, key, bm, nbm, pval, auroc, b_vals, nb_vals in predictor_stats:
    all_v = sorted(set(b_vals + nb_vals))
    best_j, best_t = 0, None
    for t in all_v:
        sens = sum(1 for v in b_vals  if v >= t) / len(b_vals)
        spec = sum(1 for v in nb_vals if v <  t) / len(nb_vals)
        j    = sens + spec - 1
        if j > best_j:
            best_j, best_t = j, t
    thresholds[label] = (best_t, best_j, auroc)

def fmt_t(label, direction, decimals, thresholds):
    t, _, auroc = thresholds.get(label, (None, 0, 0))
    if t is None:
        return "n/a", "n/a"
    fmt = f"{direction} {t:.{decimals}f}"
    return fmt, f"{auroc:.3f}"

iptm_t,    iptm_a    = fmt_t("Boltz2 ipTM",           "≥", 2, thresholds)
lis_t,     lis_a     = fmt_t("Boltz2 LIS",            "≥", 2, thresholds)
dq2_t,     dq2_a     = fmt_t("Boltz2 pDockQ2",        "≥", 3, thresholds)
iplddt_t,  iplddt_a  = fmt_t("Boltz2 complex ipLDDT", "≥", 2, thresholds)
ipsae_t,   ipsae_a   = fmt_t("Boltz2 ipSAE",          "≤", 2, thresholds)
sc_t,      sc_a      = fmt_t("Shape complementarity",  "≥", 1, thresholds)

table_data = [
    ["Metric",                "Current cutoff", "Nipah optimal", "AUROC",  "Priority"],
    ["Boltz2 ipTM",           "≥ 0.70",         iptm_t,          iptm_a,   "HIGH"],
    ["Boltz2 complex ipLDDT", "not used",        iplddt_t,        iplddt_a, "HIGH"],
    ["Shape complementarity", "not used",        sc_t,            sc_a,     "HIGH"],
    ["Boltz2 LIS",            "not used",        lis_t,           lis_a,    "MEDIUM"],
    ["Boltz2 ipSAE",          "not used",        ipsae_t,         ipsae_a,  "MEDIUM"],
    ["Boltz2 pDockQ2",        "not used",        dq2_t,           dq2_a,    "LOW"],
    ["Monomer pTM",           "≥ 0.70",          "≥ 0.70",        "—",      "HIGH"],
]

tbl = ax_tbl.table(cellText=table_data[1:], colLabels=table_data[0],
                   cellLoc="center", loc="center", bbox=[0, 0, 1, 1])
tbl.auto_set_font_size(False)
tbl.set_fontsize(8)
for (row, col), cell in tbl.get_celld().items():
    cell.set_edgecolor(GRID)
    if row == 0:
        cell.set_facecolor("#1e3a4a")
        cell.set_text_props(color=GOLD)
    else:
        cell.set_facecolor(PANEL if row % 2 == 0 else "#222222")
        cell.set_text_props(color="#eeeeee")
        # highlight priority
        if col == 4:
            txt = table_data[row][col]
            if txt == "HIGH":   cell.set_text_props(color="#81c784")
            elif txt == "MEDIUM": cell.set_text_props(color=GOLD)

fig.suptitle("Nipah Binder Competition — Retrospective Analysis\nWhat predicted experimental binding?",
             color="#eeeeee", fontsize=13, fontweight="bold", y=0.99)

out = os.path.join(os.path.dirname(__file__), "nipah_analysis.png")
plt.savefig(out, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
print(f"\nSaved {out}")
