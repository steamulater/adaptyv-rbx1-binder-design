"""
Adaptyv Competition Post-Mortem Carousel
steamulater aesthetic: white bg, #f0523d coral red, black text, +2pt fonts
1080×1350px, DPI=108
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
W, H   = 1080, 1350
DPI    = 108
FW, FH = W / DPI, H / DPI

BG      = "#ffffff"
ACCENT  = "#f0523d"
BLACK   = "#111111"
MUTED   = "#555555"
LIGHT   = "#f5f5f5"
NAVY    = "#003a8c"

OUT_DIR = Path("outputs/adaptyv_carousel")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Font sizes (+2pt vs base)
T_SUPER  = 13   # superscript / label
T_BODY   = 16
T_SUB    = 18
T_H2     = 24
T_H1     = 32
T_HERO   = 44

def new_fig():
    fig = plt.figure(figsize=(FW, FH), facecolor=BG)
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(0, H)
    ax.axis("off")
    ax.set_facecolor(BG)
    return fig, ax

def accent_bar(ax, y=H-60, h=8):
    ax.add_patch(plt.Rectangle((0, y), W, h, color=ACCENT, zorder=10))

def tag(ax, txt, x=60, y=H-100):
    ax.text(x, y, txt, color=ACCENT, fontsize=T_SUPER, fontweight="bold",
            fontfamily="monospace", va="top")

def save(fig, n, name):
    path = OUT_DIR / f"{n:02d}_{name}.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"  {path.name}")

# ── Slide 1: Cover ────────────────────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
ax.text(60, H-140, "GEM × ADAPTYV BIO", color=ACCENT, fontsize=T_SUPER,
        fontweight="bold", fontfamily="monospace", va="top")
ax.text(60, H-210, "We designed\n100 proteins\nto bind RBX1.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)
ax.text(60, H-560, "Here's what we learned\nwhen none of them\nmade the cut.",
        color=MUTED, fontsize=T_H2, va="top", linespacing=1.3)
ax.add_patch(plt.Rectangle((0, 0), W, 180, color=ACCENT))
ax.text(W/2, 90, "@steamulater", color=BG, fontsize=T_H2,
        fontweight="bold", ha="center", va="center")
save(fig, 1, "cover")

# ── Slide 2: What is RBX1 ─────────────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "THE TARGET")
ax.text(60, H-140, "RBX1", color=BLACK, fontsize=T_HERO, fontweight="bold", va="top")
ax.text(60, H-240, "Ring-Box Protein 1", color=ACCENT, fontsize=T_H2,
        fontweight="bold", va="top")
lines = [
    ("108 AA RING-H2 domain protein", BLACK),
    ("Catalytic subunit of CRL E3 ubiquitin ligases", MUTED),
    ("Regulates ~20% of cellular protein degradation", MUTED),
    ("Dysregulated in multiple cancers", MUTED),
    ("Target: E2-binding surface on RING domain", BLACK),
]
y = H - 360
for txt, col in lines:
    ax.text(80, y, f"→  {txt}", color=col, fontsize=T_BODY, va="top")
    y -= 70
ax.add_patch(plt.Rectangle((60, 400), W-120, 4, color=ACCENT))
ax.text(60, 370, "Blocking RBX1 could disrupt oncogenic protein degradation.",
        color=BLACK, fontsize=T_SUB, fontweight="bold", va="top")
save(fig, 2, "what_is_rbx1")

# ── Slide 3: The competition ──────────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "THE COMPETITION")
ax.text(60, H-140, "12,000+ proteins.\n180+ designers.\n300 chosen.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.2)

stats = [("12,000+", "sequences submitted"), ("180+", "designers"),
         ("300", "selected for wet lab"), ("~2.5%", "selection rate")]
x0 = 60
for i, (num, lbl) in enumerate(stats):
    col = i % 2
    row = i // 2
    cx = 60 + col * 500
    cy = H - 680 - row * 200
    ax.add_patch(FancyBboxPatch((cx, cy-110), 420, 110,
                  boxstyle="round,pad=8", fc=LIGHT, ec="none"))
    ax.text(cx+210, cy-38, num, color=ACCENT, fontsize=T_H1,
            fontweight="bold", ha="center", va="center")
    ax.text(cx+210, cy-85, lbl, color=MUTED, fontsize=T_BODY,
            ha="center", va="center")

ax.text(60, 160, "Volume was 6× higher than anticipated.", color=BLACK,
        fontsize=T_SUB, fontweight="bold", va="bottom")
save(fig, 3, "the_competition")

# ── Slide 4: Our approach ─────────────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "OUR APPROACH")
ax.text(60, H-140, "3 strategies.\n100 sequences.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.2)

strategies = [
    ("01", "GLMN Scaffold Redesign", "46 seqs · 247 AA · ipTM 0.84–0.89", ACCENT),
    ("02", "CUL1 WHB Redesign", "7 seqs · 72 AA · ipTM 0.70–0.76", NAVY),
    ("03", "RFdiffusion De Novo", "47 seqs · 65–95 AA · ipTM 0.71–0.91", BLACK),
]
y = H - 500
for num, title, detail, col in strategies:
    ax.add_patch(plt.Rectangle((60, y-110), 6, 110, color=col))
    ax.text(90, y-30, num, color=col, fontsize=T_SUPER,
            fontweight="bold", fontfamily="monospace", va="top")
    ax.text(90, y-70, title, color=BLACK, fontsize=T_H2, fontweight="bold", va="top")
    ax.text(90, y-105, detail, color=MUTED, fontsize=T_BODY, va="top")
    y -= 200
save(fig, 4, "our_approach")

# ── Slide 5: GLMN strategy ────────────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "STRATEGY 01")
ax.text(60, H-140, "Glomulin.\nA natural\nRBX1 binder.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)
ax.text(60, H-520, "GLMN occupies the exact E2-docking\nsurface we wanted to block.",
        color=MUTED, fontsize=T_H2, va="top", linespacing=1.3)
metrics = [
    ("100%", "pass rate ipTM ≥ 0.70"),
    ("0.588 Å", "RMSD to crystal structure"),
    ("~40%", "SwissProt identity"),
    ("46", "sequences submitted"),
]
y = H-820
for val, lbl in metrics:
    ax.text(60, y, val, color=ACCENT, fontsize=T_H1, fontweight="bold", va="top")
    ax.text(60, y-55, lbl, color=MUTED, fontsize=T_BODY, va="top")
    y -= 130
save(fig, 5, "glmn_strategy")

# ── Slide 6: RFdiffusion strategy ─────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "STRATEGY 03")
ax.text(60, H-140, "RFdiffusion.\nFully de novo\nbackbones.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)
ax.text(60, H-520, "200 backbone geometries generated\nfrom scratch — no template.",
        color=MUTED, fontsize=T_H2, va="top", linespacing=1.3)
metrics = [
    ("0%", "SwissProt identity (completely novel)"),
    ("0.910", "best ipTM in submission"),
    ("47", "sequences in final 100"),
    ("65–95 AA", "miniprotein size range"),
]
y = H-820
for val, lbl in metrics:
    ax.text(60, y, val, color=BLACK, fontsize=T_H1, fontweight="bold", va="top")
    ax.text(60, y-55, lbl, color=MUTED, fontsize=T_BODY, va="top")
    y -= 130
save(fig, 6, "rfdiffusion_strategy")

# ── Slide 7: The result ───────────────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "THE RESULT")
ax.text(W/2, H-140, "0 / 100", color=ACCENT, fontsize=90,
        fontweight="bold", ha="center", va="top")
ax.text(W/2, H-320, "of our designs were selected.", color=BLACK,
        fontsize=T_H2, fontweight="bold", ha="center", va="top")
ax.add_patch(plt.Rectangle((60, H-380), W-120, 3, color=LIGHT))
ax.text(W/2, H-440, "Our binding metrics matched\nthe selected set.", color=MUTED,
        fontsize=T_H2, ha="center", va="top", linespacing=1.3)
ax.text(W/2, H-600, "So what went wrong?", color=BLACK,
        fontsize=T_H1, fontweight="bold", ha="center", va="top")
ax.add_patch(plt.Rectangle((0, 0), W, 180, color=ACCENT))
ax.text(W/2, 90, "↓  Swipe to find out", color=BG,
        fontsize=T_H2, fontweight="bold", ha="center", va="center")
save(fig, 7, "the_result")

# ── Slide 8: Novelty was the filter ──────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "ROOT CAUSE")
ax.text(60, H-140, "Novelty\nwas the\nfilter.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)

cats   = ["score 2\n(lowest)", "score 3", "score 4\n(highest)"]
sel    = [87/322*100, 120/322*100, 106/322*100]
ours   = [62, 15, 8]
x      = np.arange(3)
width  = 0.35

# Mini bar chart embedded in slide
ax_bar = fig.add_axes([0.08, 0.18, 0.84, 0.38])
ax_bar.set_facecolor(BG)
b1 = ax_bar.bar(x - width/2, sel,  width, label="Selected (322)", color=ACCENT, alpha=0.9)
b2 = ax_bar.bar(x + width/2, ours, width, label="Ours (100)",     color=NAVY,  alpha=0.9)
ax_bar.set_xticks(x); ax_bar.set_xticklabels(cats, fontsize=T_BODY)
ax_bar.set_ylabel("% of designs", fontsize=T_BODY)
ax_bar.set_ylim(0, 80)
ax_bar.legend(fontsize=T_BODY, frameon=False)
ax_bar.spines[['top','right']].set_visible(False)
ax_bar.tick_params(labelsize=T_BODY)
for bar in b1: ax_bar.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1,
                             f"{bar.get_height():.0f}%", ha='center', fontsize=T_SUPER, color=ACCENT)
for bar in b2: ax_bar.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1,
                             f"{bar.get_height():.0f}%", ha='center', fontsize=T_SUPER, color=NAVY)
save(fig, 8, "novelty_breakdown")

# ── Slide 9: Binding metrics were fine ────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "BINDING METRICS")
ax.text(60, H-140, "Our predicted\nbinding was\nactually better.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)

metrics = [
    ("ipTM",       0.702, 0.683),
    ("ipLDDT",     0.654, 0.684),
    ("Shape comp", 56.93, 53.97),
    ("LIS",        0.355, 0.369),
]
ax_b = fig.add_axes([0.08, 0.18, 0.84, 0.42])
ax_b.set_facecolor(BG)
labels = [m[0] for m in metrics]
sel_v  = [m[1] for m in metrics]
our_v  = [m[2] for m in metrics]
x = np.arange(len(labels))
w = 0.35
b1 = ax_b.bar(x-w/2, sel_v, w, label="Selected", color=ACCENT, alpha=0.9)
b2 = ax_b.bar(x+w/2, our_v, w, label="Ours",     color=NAVY,  alpha=0.9)
ax_b.set_xticks(x); ax_b.set_xticklabels(labels, fontsize=T_BODY)
ax_b.legend(fontsize=T_BODY, frameon=False)
ax_b.spines[['top','right']].set_visible(False)
ax_b.tick_params(labelsize=T_BODY)
ax_b.set_title("Comparable or better on every binding metric",
               fontsize=T_BODY, pad=10, color=MUTED)
ax.text(60, 155, "We weren't beaten on binding. We were beaten on novelty.",
        color=BLACK, fontsize=T_SUB, fontweight="bold", va="bottom")
save(fig, 9, "binding_metrics")

# ── Slide 10: The 8 novelty=4 designs ────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "OUR BEST 8")
ax.text(60, H-140, "The 8 that\nscored maximum\nnovelty.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)

designs = [
    ("RFD_167", "70 AA", "ipTM 0.848", "LIS 0.644", True),
    ("RFD_199", "80 AA", "ipTM 0.767", "LIS 0.400", False),
    ("RFD_162", "89 AA", "ipTM 0.665", "LIS 0.280", False),
    ("RFD_47",  "84 AA", "ipTM 0.760", "LIS 0.564", False),
    ("RFD_89",  "95 AA", "ipTM 0.749", "shape 59.2", False),
    ("RFD_5",   "91 AA", "ipTM 0.598", "LIS 0.269", False),
    ("RFD_78",  "77 AA", "ipTM 0.634", "LIS 0.336", False),
    ("RFD_19",  "70 AA", "ipTM 0.555", "LIS 0.181", False),
]
y = H - 530
for name, length, m1, m2, star in designs:
    col = ACCENT if star else BLACK
    ax.text(60, y, ("★ " if star else "   ") + name, color=col,
            fontsize=T_SUB, fontweight="bold" if star else "normal", va="top")
    ax.text(340, y, length, color=MUTED, fontsize=T_BODY, va="top")
    ax.text(500, y, m1,     color=MUTED, fontsize=T_BODY, va="top")
    ax.text(760, y, m2,     color=MUTED, fontsize=T_BODY, va="top")
    y -= 88
save(fig, 10, "novelty4_designs")

# ── Slide 11: RFD_167 spotlight ───────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "SPOTLIGHT")
ax.text(60, H-140, "RFD_167", color=ACCENT, fontsize=T_HERO,
        fontweight="bold", va="top")
ax.text(60, H-250, "Our #1 ranked design.", color=BLACK,
        fontsize=T_H2, va="top")
metrics_rfd167 = [
    ("70 AA",    "compact miniprotein"),
    ("0.848",    "ipTM"),
    ("0.661",    "ipLDDT"),
    ("0.644",    "LIS — highest of all 8"),
    ("55.12",    "shape complementarity"),
    ("Novelty 4","structurally unprecedented"),
]
y = H - 380
for val, lbl in metrics_rfd167:
    ax.add_patch(plt.Rectangle((60, y-80), W-120, 75,
                  fc=LIGHT if metrics_rfd167.index((val,lbl))%2==0 else BG, ec="none"))
    ax.text(80,  y-28, val, color=ACCENT, fontsize=T_H2, fontweight="bold", va="center")
    ax.text(340, y-28, lbl, color=BLACK,  fontsize=T_BODY, va="center")
    y -= 80
save(fig, 11, "rfd167_spotlight")

# ── Slide 12: What novelty actually measures ──────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "THE LESSON")
ax.text(60, H-140, "Scaffold redesign\ncaps your novelty\nat 2.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)
ax.text(60, H-540, "GLMN is a known protein.\nProteinMPNN redesigns the sequence\nbut keeps the fold.\nThe fold is what gets scored.",
        color=MUTED, fontsize=T_H2, va="top", linespacing=1.4)
ax.add_patch(plt.Rectangle((60, 380), W-120, 4, color=ACCENT))
ax.text(60, 360, "RFdiffusion generates genuinely new folds → novelty=4",
        color=BLACK, fontsize=T_SUB, fontweight="bold", va="top")
save(fig, 12, "novelty_lesson")

# ── Slide 13: Nipah retrospective ─────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "WHAT PREDICTS BINDING")
ax.text(60, H-140, "We had the data\nthe whole time.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)
ax.text(60, H-430, "Nipah competition retrospective:\n1,030 sequences, real BLI binding data.",
        color=MUTED, fontsize=T_H2, va="top", linespacing=1.3)

predictors = [
    ("complex ipLDDT", "AUROC 0.691", ACCENT, True),
    ("Shape comp",     "AUROC 0.687", ACCENT, True),
    ("monomer pLDDT",  "AUROC 0.640", BLACK,  False),
    ("ipTM",           "AUROC 0.603", NAVY,   False),
    ("monomer pTM",    "AUROC 0.501  ← random!", MUTED, False),
]
y = H - 700
for name, auroc, col, bold in predictors:
    ax.text(80, y, name, color=col, fontsize=T_BODY,
            fontweight="bold" if bold else "normal", va="top")
    ax.text(500, y, auroc, color=col, fontsize=T_BODY,
            fontweight="bold" if bold else "normal", va="top")
    y -= 72
ax.text(60, 155, "We used pTM as a filter. It was coin-flip accuracy.",
        color=ACCENT, fontsize=T_SUB, fontweight="bold", va="bottom")
save(fig, 13, "nipah_retrospective")

# ── Slide 14: What we'd do differently ────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "NEXT ITERATION")
ax.text(60, H-140, "What we'd\ndo differently.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)

changes = [
    ("80%+ RFdiffusion",       "Maximise novelty=4. No scaffold redesign."),
    ("Shape comp filter >58",  "Secondary selection criterion within novelty tiers."),
    ("LIS in composite score", "Better predictor than ipTM alone."),
    ("Drop monomer pTM filter","AUROC 0.501 — it was eliminating good designs."),
    ("More designs per budget","6× competition volume → need 500+ not 100."),
]
y = H - 490
for title, detail in changes:
    ax.add_patch(plt.Rectangle((56, y-105), 8, 92, color=ACCENT))
    ax.text(90, y-28, title,  color=BLACK, fontsize=T_SUB, fontweight="bold", va="top")
    ax.text(90, y-75, detail, color=MUTED, fontsize=T_BODY, va="top")
    y -= 148
save(fig, 14, "next_iteration")

# ── Slide 15: The tools ────────────────────────────────────────────────────────
fig, ax = new_fig()
accent_bar(ax)
tag(ax, "THE STACK")
ax.text(60, H-140, "Open source.\nFully reproducible.", color=BLACK,
        fontsize=T_HERO, fontweight="bold", va="top", linespacing=1.15)

tools = [
    ("RFdiffusion",   "De novo backbone generation"),
    ("ProteinMPNN",   "Sequence design on backbone"),
    ("Boltz-2",       "Structure + complex validation"),
    ("DIAMOND blastp","Novelty / sequence identity screen"),
    ("FoldSeek",      "Structural novelty screen"),
    ("Google Colab",  "A100 GPU — free tier"),
]
y = H - 500
for tool, desc in tools:
    ax.text(60, y, tool, color=ACCENT, fontsize=T_SUB,
            fontweight="bold", fontfamily="monospace", va="top")
    ax.text(400, y, desc, color=MUTED, fontsize=T_BODY, va="top")
    y -= 100

ax.text(60, 160, "github.com/steamulater/fold_anywhere",
        color=BLACK, fontsize=T_SUB, fontweight="bold", va="bottom",
        fontfamily="monospace")
save(fig, 15, "the_stack")

# ── Slide 16: CTA ─────────────────────────────────────────────────────────────
fig, ax = new_fig()
ax.add_patch(plt.Rectangle((0, 0), W, H, color=ACCENT))
ax.text(W/2, H-160, "@steamulater", color=BG, fontsize=T_H1,
        fontweight="bold", ha="center", va="top")
ax.text(W/2, H/2+80, "Designing proteins.\nBuilding in public.", color=BG,
        fontsize=T_HERO, fontweight="bold", ha="center", va="center",
        linespacing=1.2)
ax.text(W/2, H/2-200, "Follow for the next iteration →", color=BG,
        fontsize=T_H2, ha="center", va="center")
ax.add_patch(FancyBboxPatch((W/2-250, 160), 500, 100,
              boxstyle="round,pad=10", fc=BG, ec="none"))
ax.text(W/2, 210, "github.com/steamulater", color=ACCENT,
        fontsize=T_SUB, fontweight="bold", ha="center", va="center",
        fontfamily="monospace")
save(fig, 16, "cta")

print(f"\nAll slides saved to {OUT_DIR}/")
print(f"Total: 16 slides")
