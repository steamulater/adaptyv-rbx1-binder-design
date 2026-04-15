# Social — RBX1 Binder Design Competition

Post-competition content for Instagram / TikTok.
All assets generated from analysis in `lab_notebook.md` (Entry 014).

---

## Reels

| File | Title | Duration | Hook |
|------|-------|----------|------|
| `transcripts/reel_01_postmortem.md` | "100 Proteins, 0 Selected" | ~90s | We had better binding metrics than the winners |
| `transcripts/reel_02_five_tools.md` | "5 Tools, 1 Target, 0 Chill" | ~90s | Pros/cons of every tool in the pipeline |

Each transcript includes: full script, visual cut guide, key facts to verify.

---

## Carousel (16 slides · 1080×1350px)

Located in `carousel/`. Steamulater aesthetic: white bg, #f0523d coral red.

| Slide | Title |
|-------|-------|
| 01 | Cover — "100 proteins, none made the cut" |
| 02 | What is RBX1 |
| 03 | The competition stats |
| 04 | Our 3 strategies |
| 05 | GLMN scaffold metrics |
| 06 | RFdiffusion strategy |
| 07 | The result — 0/100 |
| 08 | Novelty breakdown bar chart |
| 09 | Binding metrics comparison |
| 10 | The 8 novelty=4 designs |
| 11 | RFD_167 spotlight |
| 12 | The lesson — scaffold redesign caps novelty |
| 13 | Nipah retrospective (what actually predicts binding) |
| 14 | Next iteration — 5 changes |
| 15 | The stack + fold_anywhere repo |
| 16 | CTA — coral red full bleed |

Regenerate with: `python3 make_carousel_adaptyv.py` (root of this repo)

---

## Visuals

Raw assets in `visuals/` — used as B-roll or on-screen graphics in reels.

| File | What it shows |
|------|--------------|
| `Boltz_GLMN_4F52_gaps_filled_alignment.png` | Boltz-2 GLMN prediction vs crystal (RMSD 0.588Å) |
| `hotspots_for_rfdiffusion.png` | RBX1 RING domain with 14 hotspot residues (red sticks) |
| `latentx_results.png` | LatentX run summary |
| `batch2_analysis.png` | RFdiffusion batch funnel + Batch 1 vs 2 comparison |
| `rescore_analysis.png` | Old vs new filter (Nipah-derived composite) |
| `ipsae_analysis.png` | ipSAE by scaffold — why RFD_167 leads |
| `novelty_screen_results.png` | DIAMOND SwissProt identity distribution |
| `proteworks_screenshot.png` | Proteworks UI |
| `novelty4_metrics.md` | Per-design metrics for the 8 novelty=4 designs |
| `00_RBX1_target_RING_domain.pdb` | RBX1 RING domain PDB for Nanome/PyMOL |

---

## Posting notes

- Tag: `#ProteinDesign #ComputationalBiology #RFdiffusion #OpenScience #BuildingInPublic`
- Link in bio → `github.com/steamulater/fold_anywhere`
- Reel 01 is the hero post. Reel 02 is the follow-up (tools deep dive).
- Carousel can be posted independently or linked from either reel.
