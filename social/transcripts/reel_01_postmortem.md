# Reel 01 — "100 Proteins, 0 Selected"
**Format:** Instagram Reel / TikTok
**Duration:** ~90 seconds (~215 words at comfortable pace)
**Tone:** Tongue-in-cheek, scientifically accurate, self-aware

---

## SCRIPT

So I entered a global protein design competition.

The task: design a binder for RBX1 —
a protein that quietly controls the degradation
of *twenty percent* of everything inside your cells.
Cancer loves this thing.

**[CUT — RBX1 structure / carousel slide 02]**

My pipeline?
RFdiffusion to hallucinate brand new protein backbones from scratch.
ProteinMPNN to design sequences that fold into them.
Boltz-2 to validate the predicted binding —
five independent simulations per design.
All of it running on a free Google Colab A100.
At 2am.

**[CUT — Colab notebook or terminal / slide 06]**

I submitted 100 sequences.
My top design scored in the 91st percentile
for predicted binding confidence.
My interface metrics were *literally higher*
than the average of the 300 they actually selected.

And they selected zero of mine.

**[CUT — "0 / 100" slide 07]**

Here's the thing —
they weren't filtering on binding.
They were filtering on *structural novelty*.

Sixty-two of my designs were ProteinMPNN redesigns
of a natural binder.
Great predicted binding. Already-seen fold.
The algorithm clocked it immediately.

**[CUT — novelty breakdown bar chart / slide 08]**

Lesson:
if you redesign a known protein scaffold,
you inherit its novelty ceiling.
Next time — ninety percent RFdiffusion,
de novo, no templates.

I didn't get a protein in the wet lab.
But I did get a paper-worthy post-mortem.

**[CUT — CTA slide 16]**

Full breakdown in the carousel.
Follow for the next iteration — we're going again.

---

## VISUAL CUT GUIDE

| Cue | Suggested visual |
|-----|-----------------|
| "RBX1" | `social/carousel/02_what_is_rbx1.png` |
| "RFdiffusion… ProteinMPNN… Boltz-2" | `social/carousel/06_rfdiffusion_strategy.png` |
| "Colab A100 at 2am" | Screen record of Colab notebook running |
| "zero of mine" | `social/carousel/07_the_result.png` |
| "structural novelty" | `social/carousel/08_novelty_breakdown.png` |
| "novelty ceiling" | `social/carousel/12_novelty_lesson.png` |
| CTA | `social/carousel/16_cta.png` |

---

## KEY FACTS (verify before posting)

- RBX1 regulates ~20% of cellular protein degradation via CRL E3 ligase complexes
- 12,000+ sequences submitted; 322 selected (~2.5% rate)
- Our top design (RFD_167): ipTM 0.848, LIS 0.644, novelty=4
- Selected set mean ipTM: 0.702 vs ours: 0.683
- Selected set mean ipLDDT: 0.654 vs ours: 0.684 (ours higher)
- 62% of our designs scored novelty=2 vs 33% of selected set at novelty=4
