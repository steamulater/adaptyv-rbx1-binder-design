# Reel 02 — "5 Tools, 1 Target, 0 Chill"
**Format:** Instagram Reel / TikTok
**Duration:** ~90 seconds (~235 words, trim table narration if tight)
**Tone:** Tongue-in-cheek, tool-by-tool breakdown with inline pros/cons

---

## SCRIPT

I didn't just use one protein design tool.
I used five.
Because apparently I have a problem.

**[CUT — grid of 5 tool names]**

---

**Tool one. ProteinMPNN on native scaffolds.**

I took two proteins that naturally bind RBX1 —
Glomulin and a fragment of Cullin-1 —
and told ProteinMPNN to redesign their sequences from scratch.

Pro: 100% pass rate. Every sequence predicted to bind.
Con: you inherit the parent protein's fold.
Structural novelty score? Two out of four.
Ceiling hit before I even started.

**[CUT — `visuals/Boltz_GLMN_4F52_gaps_filled_alignment.png` or carousel slide 05]**

---

**Tool two. RFdiffusion.**

This one generates protein backbones from pure noise.
No template. No scaffold. Just diffusion and vibes.

Pro: genuinely unprecedented topologies. Novelty score four.
Con: x86 CUDA only.
My M-chip laptop couldn't even install it.
Had to borrow a brain from Google. At 2am.

**[CUT — `visuals/hotspots_for_rfdiffusion.png` or carousel slide 06]**

---

**Tool three. LatentX.**

Out-of-the-box generative design.
No pipeline. No YAML files. Just submit and see what comes back.

Pro: zero setup time. Genuinely plug-and-play.
Con: best ipTM we got was 0.51.
Our competition threshold was 0.70.
Beautiful idea. Wrong target, maybe.

**[CUT — `visuals/latentx_results.png`]**

---

**Tool four. Nanome.**

This is the VR one.
I put on a headset and walked *inside* the protein complex.
Stood at the RBX1 binding interface.
Counted the hotspot residues by hand.

Pro: spatial intuition that no flat screen gives you.
You feel the shape of the problem.
Con: you cannot fix your designs in virtual reality.
I tried.

**[CUT — Nanome screen record or screenshot]**

---

**Tool five. Proteworks.**

Beautiful UI. Genuinely fun to use.
Think Figma — but for proteins.

Pro: incredible for exploring sequence space and visualising design intent.
Con: we had a deadline, a Colab notebook already running,
and 100 sequences to validate.
So Proteworks became the mood board, not the engine.
No shame. It slaps.

**[CUT — `visuals/proteworks_screenshot.png`]**

---

**[CUT — pros/cons summary card or carousel slide 14]**

Five tools.
One target.
One hundred sequences.
Zero selected.

But I now know which tool does what,
which one to double down on,
and which one I'll revisit
when I have more compute and less deadline.

We're going again.

**[CUT — carousel slide 16 / CTA]**

---

## PROS/CONS SUMMARY CARD
*(use as on-screen graphic or final cut)*

| Tool | Pro | Con |
|------|-----|-----|
| ProteinMPNN | Reliable · 100% pass rate | Novelty ceiling = 2 |
| RFdiffusion | Novelty=4 · de novo topology | x86+CUDA only · no M-chip |
| LatentX | Zero setup · plug-and-play | ipTM max 0.51 (threshold 0.70) |
| Nanome VR | Spatial interface intuition | Can't fix bad science in VR |
| Proteworks | Beautiful · fun to explore | Became a vibe board, not an engine |

---

## VISUAL CUT GUIDE

| Cue | Suggested visual |
|-----|-----------------|
| "ProteinMPNN" | `social/visuals/Boltz_GLMN_4F52_gaps_filled_alignment.png` |
| "RFdiffusion" | `social/visuals/hotspots_for_rfdiffusion.png` |
| "LatentX" | `social/visuals/latentx_results.png` |
| "Nanome" | Nanome screen record from `nanome_vr/` session |
| "Proteworks" | `social/visuals/proteworks_screenshot.png` |
| Pros/cons card | Generate from table above or `social/carousel/14_next_iteration.png` |
| CTA | `social/carousel/16_cta.png` |

---

## KEY FACTS (verify before posting)

- ProteinMPNN: 100% ipTM ≥ 0.70 pass rate on GLMN scaffold; novelty=2
- RFdiffusion: 47/100 submitted; 8 scored novelty=4; top ipTM 0.910
- LatentX: 2 runs ("supreme_korsakov" 150AA, "spellbound_carter" 100AA); best ipTM 0.508
- Nanome: VR session with top 9 Boltz-2 complex structures loaded (see `nanome_vr/`)
- Proteworks: design exploration only; no sequences from Proteworks in final submission
- RFdiffusion requires x86 + CUDA — does not run on Apple Silicon (M1/M2/M3/M4)
