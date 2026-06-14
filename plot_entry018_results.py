"""
Entry 018 — RBX1 Binder Panel: Binding, Expression & Combined plots
Data from Adaptyv Foundry API retrieved 2026-06-14
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ── Data ──────────────────────────────────────────────────────────────────────

designs = [
    "PDFA_Cterm\nmpnn10",
    "PDFA_Nterm\nmpnn3",
    "GLMN\nT0.1_s11",
    "GLMN\nT0.3_s8",
    "RFD_167\nbest",
    "CUL1_WHB\nT0.2_s16",
    "PDFA_Nterm\nmpnn6",
]

short = [
    "PDFA_Cterm_mpnn10",
    "PDFA_Nterm_mpnn3",
    "GLMN_T0.1_s11",
    "GLMN_T0.3_s8",
    "RFD_167_best",
    "CUL1_WHB_T0.2_s16",
    "PDFA_Nterm_mpnn6",
]

boltz_iptm = [0.918, 0.910, 0.887, 0.878, 0.848, 0.761, 0.927]
source      = ["PDFA", "PDFA", "Steamulater", "Steamulater", "Steamulater", "Steamulater", "PDFA"]
scaffold    = ["BindCraft", "BindCraft", "GLMN", "GLMN", "RFdiffusion", "CUL1_WHB", "BindCraft"]

# BLI kinetics — only mpnn10 has values; others NaN
kd_mean  = [185e-9,  np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]  # M
kon_mean = [54591,   np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
koff_mean= [9.03e-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

# Replicate KDs for hit only (nM)
kd_reps  = [144, 126, 286]

# Expression classification (Adaptyv)
expr_label = ["high", "high", "high", "high", "medium", "high", "high"]
expr_score = {"high": 1.0, "medium": 0.6, "low": 0.2}
expr_val   = [expr_score[e] for e in expr_label]

# Binding
binding    = [True, False, False, False, False, False, False]

# Loading concentration (nM) — use 84.4/2 for "<84.4"
loading_nm = [582.1, 711.8, 442.8, 129.7, 330.0, 131.8, 42.2]

# Colours
SOURCE_COLORS = {"PDFA": "#4C72B0", "Steamulater": "#DD8452"}
EXPR_COLORS   = {"high": "#2ecc71", "medium": "#f39c12", "low": "#e74c3c"}
BIND_COLORS   = {True: "#e74c3c", False: "#bdc3c7"}

c_source  = [SOURCE_COLORS[s] for s in source]
c_expr    = [EXPR_COLORS[e]   for e in expr_label]
c_binding = [BIND_COLORS[b]   for b in binding]

x = np.arange(len(designs))
w = 0.6


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Binding data
# ══════════════════════════════════════════════════════════════════════════════
fig1, axes = plt.subplots(1, 2, figsize=(13, 5))
fig1.suptitle("RBX1 Binder Panel — Binding Data (SUL-001-002)", fontweight="bold", fontsize=13)

# Panel A: Boltz-2 iptm coloured by binding outcome
ax = axes[0]
bars = ax.bar(x, boltz_iptm, color=c_binding, edgecolor="white", linewidth=0.8, width=w)
ax.set_xticks(x); ax.set_xticklabels(designs, fontsize=8.5)
ax.set_ylabel("Boltz-2 ipTM"); ax.set_ylim(0.6, 1.0)
ax.set_title("A  Boltz-2 ipTM by Binding Outcome")
ax.axhline(0.85, ls="--", lw=1, color="#7f8c8d", label="ipTM = 0.85 threshold")
for i, (b, v) in enumerate(zip(boltz_iptm, boltz_iptm)):
    ax.text(i, v + 0.005, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
bind_patch  = mpatches.Patch(color=BIND_COLORS[True],  label="Binds (K_D 185 nM)")
nobind_patch= mpatches.Patch(color=BIND_COLORS[False], label="No binding")
ax.legend(handles=[bind_patch, nobind_patch, ax.get_lines()[0]], fontsize=8.5)

# Panel B: replicate KDs for hit
ax = axes[1]
reps = [1, 2, 3]
ax.bar(reps, kd_reps, color="#4C72B0", edgecolor="white", linewidth=0.8, width=0.5)
ax.axhline(185, ls="--", lw=1.5, color="#e74c3c", label=f"Mean K_D = 185 nM")
for i, v in enumerate(kd_reps):
    ax.text(i+1, v+4, f"{v} nM", ha="center", va="bottom", fontsize=9)
ax.set_xticks(reps); ax.set_xticklabels(["Rep 1", "Rep 2", "Rep 3"])
ax.set_ylabel("K_D (nM)"); ax.set_ylim(0, 350)
ax.set_title("B  PDFA_Cterm_mpnn10 — Replicate K_D (BLI)")
ax.legend(fontsize=9)

plt.tight_layout()
fig1.savefig("rbx1_binding_data_entry018.png", dpi=180, bbox_inches="tight")
plt.close(fig1)
print("Saved: rbx1_binding_data_entry018.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Expression data
# ══════════════════════════════════════════════════════════════════════════════
fig2, axes = plt.subplots(1, 2, figsize=(13, 5))
fig2.suptitle("RBX1 Binder Panel — Expression Data (SUL-001-002)", fontweight="bold", fontsize=13)

# Panel A: expression classification bar
ax = axes[0]
ax.bar(x, expr_val, color=c_expr, edgecolor="white", linewidth=0.8, width=w)
ax.set_xticks(x); ax.set_xticklabels(designs, fontsize=8.5)
ax.set_ylabel("Expression level"); ax.set_ylim(0, 1.25)
ax.set_yticks([0.2, 0.6, 1.0]); ax.set_yticklabels(["Low", "Medium", "High"])
ax.set_title("A  Adaptyv Expression Classification (all replicates)")
for i, (e, v) in enumerate(zip(expr_label, expr_val)):
    ax.text(i, v + 0.03, e.capitalize(), ha="center", va="bottom", fontsize=8.5, fontweight="bold")
high_p  = mpatches.Patch(color=EXPR_COLORS["high"],   label="High")
med_p   = mpatches.Patch(color=EXPR_COLORS["medium"], label="Medium")
low_p   = mpatches.Patch(color=EXPR_COLORS["low"],    label="Low")
ax.legend(handles=[high_p, med_p, low_p], fontsize=8.5, title="Expression")

# Panel B: BLI tip loading concentration
ax = axes[1]
bars = ax.bar(x, loading_nm, color=c_source, edgecolor="white", linewidth=0.8, width=w)
ax.set_xticks(x); ax.set_xticklabels(designs, fontsize=8.5)
ax.set_ylabel("Tip loading concentration (nM)")
ax.set_title("B  BLI Tip Loading Concentration")
for i, (v, label) in enumerate(zip(loading_nm, short)):
    txt = f"<84" if label == "PDFA_Nterm_mpnn6" else f"{v:.0f}"
    ax.text(i, v + 10, txt, ha="center", va="bottom", fontsize=8.5)
pdfa_p  = mpatches.Patch(color=SOURCE_COLORS["PDFA"],        label="PDFA")
steam_p = mpatches.Patch(color=SOURCE_COLORS["Steamulater"], label="Steamulater")
ax.legend(handles=[pdfa_p, steam_p], fontsize=8.5, title="Source")

plt.tight_layout()
fig2.savefig("rbx1_expression_data_entry018.png", dpi=180, bbox_inches="tight")
plt.close(fig2)
print("Saved: rbx1_expression_data_entry018.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Combined overview
# ══════════════════════════════════════════════════════════════════════════════
fig3, axes = plt.subplots(2, 2, figsize=(14, 10))
fig3.suptitle("RBX1 Binder Panel — Full Experimental Summary (SUL-001-002, Entry 018)",
              fontweight="bold", fontsize=13)

# ── Top-left: ipTM coloured by source
ax = axes[0, 0]
ax.bar(x, boltz_iptm, color=c_source, edgecolor="white", linewidth=0.8, width=w)
ax.set_xticks(x); ax.set_xticklabels(designs, fontsize=8)
ax.set_ylabel("Boltz-2 ipTM"); ax.set_ylim(0.6, 1.0)
ax.set_title("Boltz-2 ipTM by Source")
ax.axhline(0.85, ls="--", lw=1, color="#7f8c8d")
for i, v in enumerate(boltz_iptm):
    ax.text(i, v + 0.005, f"{v:.3f}", ha="center", fontsize=7.5)
ax.legend(handles=[pdfa_p, steam_p], fontsize=8)

# ── Top-right: expression classification
ax = axes[0, 1]
ax.bar(x, expr_val, color=c_expr, edgecolor="white", linewidth=0.8, width=w)
ax.set_xticks(x); ax.set_xticklabels(designs, fontsize=8)
ax.set_ylabel("Expression"); ax.set_ylim(0, 1.3)
ax.set_yticks([0.2, 0.6, 1.0]); ax.set_yticklabels(["Low", "Medium", "High"])
ax.set_title("Adaptyv Expression Classification")
for i, e in enumerate(expr_label):
    ax.text(i, expr_score[e] + 0.04, e.capitalize(), ha="center", fontsize=7.5, fontweight="bold")
ax.legend(handles=[high_p, med_p, low_p], fontsize=8)

# ── Bottom-left: loading concentration vs iptm scatter
ax = axes[1, 0]
for i, (iptm, load, bind, src) in enumerate(zip(boltz_iptm, loading_nm, binding, source)):
    marker = "★" if bind else "o"
    ax.scatter(iptm, load, s=140 if bind else 80,
               color=SOURCE_COLORS[src],
               edgecolor="#2c3e50" if bind else "none",
               linewidth=1.5, zorder=3)
    ax.annotate(designs[i].replace("\n", " "), (iptm, load),
                textcoords="offset points", xytext=(5, 3), fontsize=7)
ax.set_xlabel("Boltz-2 ipTM"); ax.set_ylabel("Tip loading concentration (nM)")
ax.set_title("ipTM vs. Tip Loading  (★ = binder)")
ax.legend(handles=[pdfa_p, steam_p], fontsize=8)

# ── Bottom-right: summary heatmap (ipTM / expression / binding / loading)
ax = axes[1, 1]
metrics = np.array([
    boltz_iptm,
    expr_val,
    [1.0 if b else 0.0 for b in binding],
    [min(v / 800, 1.0) for v in loading_nm],
])
metric_labels = ["Boltz-2 ipTM\n(scaled 0.6–1)", "Expression\n(H/M/L)", "Binding\n(1=yes)", "Loading conc\n(norm. /800nM)"]
# Normalise ipTM to 0-1 for display
metrics[0] = (np.array(boltz_iptm) - 0.6) / 0.4

im = ax.imshow(metrics, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
ax.set_xticks(x); ax.set_xticklabels([d.replace("\n", " ") for d in designs], fontsize=7.5, rotation=30, ha="right")
ax.set_yticks(range(4)); ax.set_yticklabels(metric_labels, fontsize=8.5)
ax.set_title("Multi-metric Heatmap")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Score (0–1)")

# Annotate cells
cell_vals_raw = [boltz_iptm, expr_val, [int(b) for b in binding], [min(v/800,1) for v in loading_nm]]
cell_labels   = [
    [f"{v:.3f}" for v in boltz_iptm],
    [e.capitalize()[0] for e in expr_label],
    ["✓" if b else "✗" for b in binding],
    [f"{v:.0f}" if v > 42 else "<84" for v in loading_nm],
]
for r in range(4):
    for c in range(len(designs)):
        ax.text(c, r, cell_labels[r][c], ha="center", va="center",
                fontsize=8, color="black", fontweight="bold")

plt.tight_layout()
fig3.savefig("rbx1_combined_summary_entry018.png", dpi=180, bbox_inches="tight")
plt.close(fig3)
print("Saved: rbx1_combined_summary_entry018.png")
print("\nAll plots complete.")
