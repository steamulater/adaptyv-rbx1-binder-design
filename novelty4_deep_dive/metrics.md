# Novelty=4 Designs — Deep Dive
**8 designs from steamulater submission that scored maximum structural novelty**

All are RFdiffusion de novo miniproteins. These were our highest-novelty designs
but were not selected for wet lab (300/12,000+ chosen overall).

---

| Design | Rank | Length | ipTM | ipLDDT | ipSAE | Shape Comp | LIS | MW (Da) |
|--------|------|--------|------|--------|-------|-----------|-----|---------|
| RFD_167 | 1 | 70 AA | 0.848 | 0.661 | 0.619 | 55.12 | 0.644 | 8,504 |
| RFD_199 | 3 | 80 AA | 0.767 | 0.524 | 0.315 | 51.65 | 0.400 | 7,841 |
| RFD_162 | 8 | 89 AA | 0.665 | 0.651 | 0.270 | 50.52 | 0.280 | 10,263 |
| RFD_47  | 52 | 84 AA | 0.760 | 0.791 | 0.631 | 47.83 | 0.564 | 9,388 |
| RFD_89  | 45 | 95 AA | 0.749 | 0.597 | 0.450 | 59.16 | 0.503 | 10,354 |
| RFD_5   | 65 | 91 AA | 0.598 | 0.616 | 0.275 | 52.30 | 0.269 | 8,911 |
| RFD_78  | 82 | 77 AA | 0.634 | 0.625 | 0.292 | 54.15 | 0.336 | 7,148 |
| RFD_19  | 46 | 70 AA | 0.555 | 0.606 | 0.111 | 54.93 | 0.181 | 6,738 |

---

## Metric definitions (Proteinbase / Boltz-2)

- **ipTM**: interface predicted Template Modeling score (0–1). Higher = more confident interface geometry.
- **ipLDDT**: interface predicted Local Distance Difference Test (0–1). Confidence in interface residue positions.
- **ipSAE**: interface Structural Assessment Error from PAE matrix (Ångstroms, normalised 0–1 here). Lower raw = better; Proteinbase reports normalised.
- **Shape Comp**: shape complementarity score at the binder–RBX1 interface (Proteinbase metric, higher = better fit).
- **LIS**: Ligand Interface Score. Composite interface quality metric.
- **MW**: molecular weight in Daltons.

---

## Standouts

**RFD_167** — best overall. 70 AA, ipTM=0.848, highest LIS (0.644). Ranked #1 in our submission.
  → Structure: `RFD_167_complex.cif`

**RFD_47** — best ipLDDT (0.791) and highest LIS among mid-rankers. 84 AA.
  → Structure: `RFD_47_complex.cif`

**RFD_89** — highest shape complementarity (59.16). 95 AA.
  → Structure: `RFD_89_complex.cif`

**RFD_19** — smallest (70 AA, 6.7 kDa). Lowest LIS and ipSAE — weakest of the 8 but most compact.
  → Structure: `RFD_19_complex.cif`

---

## Why they weren't selected despite novelty=4

The selected set (322 designs) scored comparably on binding metrics:
- Selected mean ipTM: 0.702 vs ours: 0.683
- Selected mean ipLDDT: 0.654 vs ours: 0.684 (we were HIGHER)
- Selected mean shape comp: 56.93 vs ours: 53.97

The novelty=4 subset of ours (these 8) likely competed directly with the
~106 novelty=4 designs in the selected set. Selection among equally novel
designs appears to favour **shape complementarity** and **LIS** — metrics
we did not optimise for during design.

Hypothesis for next iteration: after RFdiffusion + MPNN, add a shape
complementarity filter (target >58) before final ranking.
