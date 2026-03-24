# RBX1 Binder Design — Submission Write-up
## GEM x Adaptyv Bio RBX1 Binder Design Challenge
**Submission date:** March 23, 2026 (Batch 1 of 2)
**Sequences submitted:** 55
**Designer:** Tamuka Martin

---

## Overview

We submit 55 de novo designed protein binders targeting the RBX1 RING domain, generated through a structure-guided computational pipeline combining scaffold selection, ProteinMPNN sequence design, and Boltz-2 structure validation. All sequences are novel (<75% identity to any SwissProt entry) and achieve high predicted interface confidence (avg ipTM 0.867 for the primary scaffold class).

---

## Target

**RBX1 (RING-box protein 1, UniProt P62877)**
- Full length: 108 AA
- Structured RING-H2 domain: residues 32–108 (77 AA)
- Residues 1–31: intrinsically disordered, excluded from interface analysis
- Function: E3 ubiquitin ligase scaffold; recruits E2 ubiquitin-conjugating enzymes via the RING domain
- Design objective: bind the RING domain face that is exposed in the free form (same face engaged by natural RING-domain interactors)

---

## Design Strategy

### Scaffold Selection — Structure-Guided Mimic Approach

Rather than generating binders from random noise (pure de novo), we selected two structurally validated scaffolds that have been experimentally observed to engage the RBX1 RING domain:

**Scaffold 1 — GLMN (Glomulin, PDB 4F52)**
- GLMN is a 2386 AA scaffold protein in the CRL2A ubiquitin ligase complex
- Residues 336–582 (247 AA) constitute the minimal RBX1-contacting domain extracted from PDB 4F52
- This fragment packs directly against the RBX1 RING domain in the native complex
- Rationale: use the natural binding geometry as a structural template; ProteinMPNN redesigns the sequence to generate novel binders that retain the native fold and interface geometry

**Scaffold 2 — CUL1 WHB domain (PDB 1LDJ, residues 705–776)**
- The winged-helix B (WHB) domain of Cullin-1 is a 72 AA fragment that cradles RBX1 in the SCF complex
- Rationale: compact alternative scaffold with a different contact geometry to GLMN
- Limitation discovered post-validation: without the full Cullin scaffold, the isolated WHB fragment cannot constrain the RBX1 geometry and requires RBX1 rearrangement to bind (see Validation section)

### Sequence Design — ProteinMPNN

ProteinMPNN (Dauparas et al. 2022) was used to redesign sequences on each scaffold backbone. The interface residues (contacts within 5 Å of RBX1 in the native complex) were fixed as partial constraints; all other positions were redesigned.

**Design parameters:**
- Temperatures: T = 0.1 (conservative), 0.2 (balanced), 0.3 (diverse)
- Samples per temperature: 16
- Total sequences generated: 48 per scaffold × 2 scaffolds = 96 unique sequences

---

## Computational Validation

### Boltz-2 Structure Prediction

All 96 sequences were validated using Boltz-2 (single-sequence mode, no MSA, 5 diffusion samples per prediction):

- **Monomer prediction**: assess whether the designed sequence folds into the target scaffold structure (pTM, pLDDT)
- **Complex prediction with RBX1**: assess binding interface quality (ipTM)

**Boltz-2 results summary:**

| Scaffold | n | avg ipTM | avg pTM | avg pLDDT | Sequences ≥ 0.70 ipTM |
|----------|---|----------|---------|-----------|----------------------|
| GLMN | 48 | 0.867 | 0.889 | 0.823 | 48/48 (100%) |
| CUL1_WHB | 48 | 0.595 | 0.942 | 0.941 | 7/48 (15%) |

**Key finding:** Every GLMN redesign achieved ipTM ≥ 0.844, demonstrating robust predicted binding to RBX1 across all temperatures and samples. CUL1_WHB showed highly variable confidence attributable to induced-fit binding (see below).

### RBX1 RING Domain RMSD Analysis

To distinguish *locked-in* (preformed interface) from *induced-fit* (target rearrangement) binding modes, we computed the RMSD of RBX1 RING domain Cα atoms (residues 32–108) in each complex prediction versus the Boltz-2 native reference structure:

| Scaffold | Mean RING RMSD | Interpretation |
|----------|---------------|----------------|
| GLMN | 1.09 ± 0.28 Å | **Locked-in** — RBX1 maintains native geometry |
| CUL1_WHB | 5.11 ± 1.05 Å | **Induced-fit** — RBX1 deforms to accommodate binder |

**Implication:** GLMN-based binders engage RBX1 without requiring structural rearrangement. This is the mechanistically preferred mode: the binder captures the native RBX1 conformation, lowering the entropic cost of binding and improving experimental success probability.

**Correlation:** ipTM and RING RMSD show r = −0.84 across all 96 sequences, confirming that better predicted binding is directly linked to preservation of RBX1 native geometry.

### Novelty Screen

All sequences were screened using DIAMOND blastp (sensitive mode) against UniProt/SwissProt (574,627 sequences):

| Scaffold | Mean identity to best hit | All pass (<75%)? |
|----------|--------------------------|-----------------|
| GLMN | 40.0 ± 1.8% | YES (48/48) |
| CUL1_WHB | 42.1 ± 3.3% | YES (7/7) |

The ~40% SwissProt identity confirms that ProteinMPNN generates genuinely novel sequences. Despite using GLMN structure as a scaffold, the redesigned sequences are no more similar to natural GLMN than two distantly related proteins — well below the 75% novelty threshold.

---

## Submitted Sequences (55 total)

### Primary class — GLMN scaffold (48 sequences)

All 48 GLMN-based sequences are submitted. Selection criteria:
- ipTM ≥ 0.844 (all pass)
- RING RMSD < 2.0 Å (all pass, mean 1.09 Å)
- SwissProt identity < 75% (all pass, mean 40.0%)
- Length: 247 AA

**Top 10 by ipTM:**

| Rank | Seq ID | ipTM | pTM | RING RMSD |
|------|--------|------|-----|-----------|
| 1 | GLMN_T0.1_s11 | 0.887 | 0.889 | 0.96 Å |
| 2 | GLMN_T0.1_s4 | 0.886 | 0.891 | 1.04 Å |
| 3 | GLMN_T0.3_s12 | 0.882 | 0.879 | 0.91 Å |
| 4 | GLMN_T0.1_s15 | 0.881 | 0.895 | 1.12 Å |
| 5 | GLMN_T0.2_s5 | 0.879 | 0.891 | 0.89 Å |
| 6 | GLMN_T0.2_s14 | 0.879 | 0.891 | 0.92 Å |
| 7 | GLMN_T0.1_s13 | 0.878 | 0.888 | 1.03 Å |
| 8 | GLMN_T0.3_s8 | 0.878 | 0.896 | 0.60 Å |
| 9 | GLMN_T0.3_s4 | 0.878 | 0.888 | 1.18 Å |
| 10 | GLMN_T0.1_s10 | 0.877 | 0.873 | 1.22 Å |

**Temperature diversity:** 16 sequences each at T=0.1, T=0.2, T=0.3, ensuring chemical diversity across the 48 submitted sequences.

### Secondary class — CUL1 WHB scaffold (7 sequences)

The 7 CUL1_WHB sequences with average ipTM ≥ 0.70 are included to provide scaffold diversity and represent a fundamentally different approach to RBX1 binding:

| Seq ID | avg ipTM | Best single-model ipTM | RING RMSD |
|--------|----------|----------------------|-----------|
| CUL1_WHB_T0.2_s16 | 0.761 | ~0.93 | 5.2 Å |
| CUL1_WHB_T0.1_s8 | 0.759 | ~0.92 | 5.3 Å |
| CUL1_WHB_T0.1_s4 | 0.748 | ~0.90 | 5.1 Å |
| CUL1_WHB_T0.3_s11 | 0.727 | ~0.88 | 5.8 Å |
| CUL1_WHB_T0.1_s6 | 0.725 | ~0.87 | 5.4 Å |
| CUL1_WHB_T0.1_s7 | 0.711 | ~0.85 | 5.7 Å |
| CUL1_WHB_T0.2_s3 | 0.701 | ~0.84 | 5.4 Å |

These compact 72 AA binders represent a distinct binding mode and may be valuable as alternative confirmatory binders or as starting points for further optimization.

---

## Design Philosophy

This submission uses **structure-guided scaffold mimicry** — identifying natural proteins that already occupy the target binding site, then using ProteinMPNN to explore the sequence landscape around that structural solution. This approach:

1. Guarantees the starting point is a structurally viable binding geometry (unlike random RFdiffusion, which may generate backbones that are hard to sequence-design)
2. Produces binders with experimentally interpretable interfaces (based on known PDB structures)
3. Generates high diversity (>40 sequences per scaffold class, spanning 3 temperature regimes) while maintaining high confidence scores

**Limitations of current batch:**
- All 48 GLMN binders share the same backbone scaffold — sequence diversity exists but structural diversity is limited to one binding mode
- The CUL1_WHB binders show induced-fit binding that may be harder to validate experimentally

**Batch 2 (coming):** We will supplement this submission with 45 sequences from a RFdiffusion de novo backbone generation strategy, providing true structural diversity that is independent of the GLMN scaffold. This will cover a broader region of design space and reduce the risk of the entire submission failing for a single scaffold-specific reason.

---

## Files

| File | Description |
|------|-------------|
| `final_submission.fasta` | All 55 sequences in FASTA format, ranked by ipTM |
| `master_sequences.csv` | Complete metadata: ipTM, pTM, pLDDT, RING RMSD, novelty |
| `boltz_confidence_results.json` | Raw Boltz-2 confidence scores (all 96 sequences × 5 models) |
| `novelty_swissprot.tsv` | DIAMOND blastp results vs SwissProt |
| `novelty_screen_results.csv` | Parsed novelty screen with pass/fail |
| `boltz_results_overview.png` | Boltz-2 validation summary figures |
| `rbx1_rmsd_analysis.png` | RBX1 RING RMSD analysis figures |
| `novelty_screen_results.png` | Novelty screen identity distribution |

---

## References

1. Dauparas J et al. (2022). Robust deep learning-based protein sequence design using ProteinMPNN. *Science* 378:49–56.
2. Wohlgemuth N et al. (2025). Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction. *bioRxiv*.
3. Duda DM et al. (2012). Structure of a glomulin-RBX1-CUL1 complex: inhibition of a RING E3 ligase through masking of its E2-binding surface. *Mol Cell* 47:371–382. (PDB: 4F52)
4. Zheng N et al. (2002). Structure of the Cul1-Rbx1-Skp1-F boxSkp2 SCF ubiquitin ligase complex. *Nature* 416:703–709. (PDB: 1LDJ)
