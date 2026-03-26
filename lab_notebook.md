# Lab Notebook: RBX1 Binder Design
**Competition:** GEM x Adaptyv Bio — RBX1 Binder Design Challenge
**Submission deadline:** March 26, 2026
**Experimenter:** Tamuka Martin Chidyausiku @steamulater

---

## Project Overview

Design up to 100 *de novo* protein binder sequences targeting RBX1 (Ring-Box Protein 1) for experimental validation via bio-layer interferometry at Adaptyv's Foundry.

**Constraints:**
- Max 250 amino acids per sequence
- Min 25% edit distance from UniRef50 (or SAbDab for sdAbs)
- Up to 100 sequences

---

## Target: RBX1

**UniProt:** P62877 | **PDB:** 2LGV | **Organism:** *Homo sapiens*

**Full sequence (108 AA):**
```
GGGGTNSGAGKKRFEVKKSNASAQSAWDIVVDNCAICRNHIMDLCIECQANQASATSEECTVAWGVCNHAFHFHCISRWLKTRQVCPLDNREWEFQKYGH
```

**Domain architecture:**
| Region | Residues | Description |
|--------|----------|-------------|
| Disordered N-terminus | 1–31 | Intrinsically disordered, poor binder target |
| RING-H2 domain | 32–90 | Coordinates 3 Zn2+ ions; primary target |
| RING catalytic core | 45–90 | E2 ubiquitin-conjugating enzyme binding surface |
| C-terminal tail | 91–108 | — |

**Function:** Catalytic subunit of Cullin-RING E3 ubiquitin ligase (CRL) complexes. Recruits E2 enzymes to drive ubiquitination of ~20% of cellular proteins, regulating cell cycle, DNA repair, and signal transduction.

**Therapeutic rationale:** CRL activity is frequently dysregulated in cancer. Blocking the RBX1 E2-binding surface on the RING-H2 domain disrupts ubiquitin transfer and could inhibit oncogenic CRL-mediated protein degradation.

**Primary epitope target:** The E2-binding surface on the RING-H2 domain (residues ~45–90), particularly the zinc-coordinating residues and the E2-docking loops.

---

## Design Strategies

Two parallel strategies are being pursued to maximise the diversity and quality of the 100 submissions.

---

### Strategy 1: De Novo Backbone Generation (RFdiffusion + ProteinMPNN)

**Status: COMPLETE — 47 sequences submitted**

Generate entirely new binder backbones conditioned on the RBX1 RING-H2 surface, then design sequences with ProteinMPNN.

```
RBX1 RING domain (rbx1_ring_renumbered.pdb, chain A res 1-77)
        |
        v
[1] RFdiffusion  [Colab A100]  ✓ DONE
    - Hotspot residues: A4,A6,A12,A14,A15,A23,A24,A26,A28,A52,A56,A60,A64,A66
    - Contigs: "A1-77/0 60-100"
    - Generated 200 backbones → 151 pass length filter (65-95 AA)
        |
        v
[2] ProteinMPNN  [Colab A100]  ✓ DONE
    - Chain A = binder (designed), Chain B = RBX1 (fixed)
    - Temperatures 0.1, 0.2 | 2 sequences per backbone
    - 604 sequences → best-per-backbone → 151 candidates
        |
        v
[3] Boltz-2 validation  [Colab A100]  ✓ DONE
    - Monomer: pTM, pLDDT (118/151 completed)
    - Complex: ipTM, complex_iplddt, ipSAE (151/151 completed)
        |
        v
[4] Filter & rank  ✓ DONE
    - ipTM ≥ 0.70 → 57/151 pass
    - Poly-Ala terminal filter → 56 clean
    - Rank by composite: 0.4×ipTM + 0.3×ipLDDT + 0.3×norm_ipSAE
        |
        v
[5] Novelty screen  ✓ DONE
    - DIAMOND blastp vs SwissProt: 0/47 hits — completely novel
        |
        v
[6] Top 47 submitted (of 56 passing, capped at 100 total)
```

**Results summary:**
- 200 backbones → 47 submitted sequences
- ipTM range: 0.710–0.910 | composite range: 0.752–0.845
- Mean ipSAE: 8.6 Å (best: RFD_167 at 5.56 Å — more confident interface than any Batch 1 sequence)
- Novelty: 0% SwissProt identity (no hits at e-value 1e-3)

**Rationale:** Maximally novel binders; not constrained to known binding modes. RFdiffusion generates completely different backbone geometries from GLMN, providing structural diversity that de-risks the overall submission.

---

### Strategy 2a: Scaffold-Based Redesign — Glomulin (4F52)

**Status: COMPLETE — 48 sequences generated, 46 submitted (top 46 of 48 by composite score)**

Use the experimentally validated Glomulin backbone — a natural inhibitor of CRL activity that occupies the E2-docking surface of RBX1's RING domain — as the starting scaffold for ProteinMPNN sequence redesign.

**Scaffold:** PDB 4F52, chain F (Glomulin), residues 336–582 | 247 AA | UniProt Q92990

```
PDB 4F52 (Glomulin-RBX1-CUL1 crystal structure)
        |
        v
[1] PyMOL  [local]  ✓ DONE
    - Extract Glomulin chain F, residues 336-582 (247 AA)
    - Alignment RMSD to native: 0.588 Å (near-perfect)
        |
        v
[2] ProteinMPNN  [local, M3 MPS]  ✓ DONE
    - Fixed receptor: RBX1 | Redesign Glomulin sequence
    - Temperatures 0.1, 0.2, 0.3 | 16 sequences per temperature → 48 total
        |
        v
[3] Boltz-2 validation  [Colab A100]  ✓ DONE
    - 5 diffusion samples per sequence (monomer + complex)
    - All 48/48 pass ipTM ≥ 0.70 | mean ipTM = 0.867
    - RBX1 RING RMSD: mean 1.09 Å (locked-in binding, not induced-fit)
        |
        v
[4] Novelty screen  ✓ DONE
    - DIAMOND blastp vs SwissProt: all 48 pass (<75% identity)
    - Mean SwissProt identity: 40.0% (genuinely redesigned)
        |
        v
[5] 46 submitted (2 displaced by higher-scoring RFdiffusion sequences in final ranking)
```

**Results summary:**
- 48 sequences generated → 46 submitted
- ipTM: 0.844–0.887 (100% pass rate — every GLMN redesign binds)
- RBX1 RING RMSD: 0.60–1.22 Å (locked-in geometry, not induced-fit)
- Mean ipSAE: 9.06 Å (7.72–10.79 Å)
- Novelty: ~40% SwissProt identity (all pass <75% threshold)

**Key structural references:**
- 4F52 scaffold: chain F, residues 336–582
- Aligned to 2LGV (RBX1 NMR): RMSD = 7.692 Å on 88 atoms

![4F52 Glomulin scaffold in PyMOL — orange: Glomulin chain F res 336-582; green: CRL complex](Scaffold_4F52_582-336.png)

---

### Strategy 2b: Scaffold-Based Redesign — CUL1 WHB domain (1LDJ)

**Status: COMPLETE — 48 sequences generated, 7 submitted**

Use the C-terminal winged-helix B (WHB) domain of Cullin-1 that directly cradles the RBX1 RING domain as a compact (72 AA) scaffold for ProteinMPNN redesign.

**Scaffold:** PDB 1LDJ, chain A (CUL1), residues 705–776 | 72 AA | UniProt Q13616

```
PDB 1LDJ (CUL1-RBX1-SKP1-SKP2 SCF complex)
        |
        v
[1] PyMOL  [local]  ✓ DONE
    - Extract CUL1 chain A, residues 705-776 (72 AA)
    - PyMOL alignment RMSD to native 1LDJ: 4.844 Å (bimodal — induced fit)
        |
        v
[2] ProteinMPNN  [local, M3 MPS]  ✓ DONE
    - Temperatures 0.1, 0.2, 0.3 | 16 sequences per temperature → 48 total
        |
        v
[3] Boltz-2 validation  [Colab A100]  ✓ DONE
    - 5 diffusion samples per sequence (monomer + complex)
    - Only 7/48 pass ipTM ≥ 0.70 (15% pass rate)
    - RBX1 RING RMSD: mean 5.11 Å (induced-fit — WHB cannot bind without full Cullin scaffold)
        |
        v
[4] Novelty screen  ✓ DONE
    - All 7 pass DIAMOND blastp vs SwissProt (<75% identity)
        |
        v
[5] 7 submitted
```

**Results summary:**
- 48 sequences generated → 7 submitted (15% pass rate vs 100% for GLMN)
- ipTM: 0.701–0.761 (marginal — near threshold)
- RBX1 RING RMSD: 4.8–5.8 Å — induced-fit binding mode, structurally problematic
- Mean ipSAE: 12.47 Å (8.82–17.66 Å — worst of all scaffolds)
- Included for scaffold diversity; lower experimental success probability than GLMN

**Key lesson:** Without the full Cullin-1 scaffold, the isolated WHB fragment cannot constrain the RBX1 geometry. The 7 sequences included represent an alternative binding mode hypothesis, not a high-confidence scaffold class.

**Key structural references:**
- 1LDJ scaffold: chain A, residues 705–776

![1LDJ CUL1 WHB domain scaffold in PyMOL — yellow: CUL1 chain A res 705-776; red: RBX1](Scaffold_1LDJ_705-776.png)

---

## Experiment Log

---

### Entry 001 — 2026-03-22

**Status:** Project initialized

**Work completed:**
- Reviewed RBX1 biology and domain architecture via BioReason analysis
- Identified RING-H2 domain (res 32–90) as primary binding epitope
- Identified E2-binding surface (res 45–90) as key functional hotspot
- Retrieved full-length structure PDB: 2LGV, inspected in PyMOL
- Confirmed disordered N-terminus (res 1–31) is a poor binder target
- Established computational pipeline: RFdiffusion → ProteinMPNN → ColabFold → filter → rank

**Key structural observations from 2LGV:**
- RING-H2 domain is well-structured and amenable to binder design
- Three zinc coordination sites create a rigid scaffold
- E2-docking surface is solvent-exposed and flat — will need binders with concave interface or extended loops
- N-terminal ~31 residues are disordered in solution (NMR structure)

**Next steps:**
- [ ] Prepare RBX1 input structure for RFdiffusion (extract RING-H2 domain, define hotspot residues)
- [ ] Run RFdiffusion backbone generation on Colab
- [ ] Run ProteinMPNN sequence design on generated backbones
- [ ] Validate with ColabFold AF2-multimer
- [ ] Filter for novelty and rank top 100

---

### Entry 002 — 2026-03-22

**Status:** Complete

**Work completed:**
- Identified confirmed PDB structures containing human RBX1 (4F52, 2HYE, 5N4W, 8Q7H, 1LDJ)
- Selected 4F52 (Glomulin-RBX1-CUL1) as scaffold source for Strategy 2a
- Loaded 4F52 in PyMOL, aligned against 2LGV (RMSD = 7.692 Å on 88 atoms)
- Extracted Glomulin chain F residues 336–582 (247 AA) → `Scaffold_4F52_336-582.pdb`
- Found 2 gaps in crystal structure: res 433–438 (6 AA) and 533–549 (17 AA)
- Selected 1LDJ CUL1 WHB domain (chain A, res 705–776, 72 AA) as Strategy 2b scaffold → `Scaffold_1LDJ_705-776.pdb`
- Decided to use Boltz-1 on ColabFold to model complete structures (fills missing loops better than ESMFold)
- Set up local ProteinMPNN environment (conda env: `proteinmpnn`, PyTorch 2.10, MPS enabled)
- Created 4 Boltz-1 YAML input files in `boltz_inputs/`

**Key sequences confirmed:**
- RBX1 (P62877): 108 AA — `MAAAMDVDTPSGTNSGAGKKRFEVKKWNAVALWAW...`
- Glomulin 336–582 (Q92990): 247 AA
- CUL1 WHB 705–776 (Q13616): 72 AA — `EDRKLLIQAAIVRIMKMRKVLKHQQLLGEVLTQ...`

**Next steps:**
- [ ] Run 4 Boltz-1 jobs on ColabFold (5 diffusion samples each)
- [ ] Assess ipTM in complex predictions — confirm both scaffolds form PPI with RBX1
- [ ] Use best Boltz monomer outputs as ProteinMPNN scaffolds
- [ ] Run ProteinMPNN redesign locally (Strategies 2a and 2b)
- [ ] Set up RFdiffusion Colab run (Strategy 1)
- [ ] Validate all designs with ColabFold AF2-multimer
- [ ] Filter for novelty and rank top 100

---

### Entry 003 — 2026-03-22

**Status:** Complete

**Work completed:**
- BLAST verified all three sequences against SwissProt (NCBI blastp)
- All sequences confirmed 100% identity to correct human proteins at correct positions

**BLAST results:**

| Sequence | BLAST Top Hit | Identity | Position | NCBI RID |
|----------|--------------|----------|----------|----------|
| RBX1 | P62877 — E3 ubiquitin-protein ligase RBX1, *H. sapiens* | 108/108 (100%) | 1–108 | W1DWTW13014 |
| Glomulin 336–582 | Q92990 — Glomulin (FAP), *H. sapiens* | 247/247 (100%) | 336–582 | W1DWUCET014 |
| CUL1 705–776 | Q13616 — Cullin-1, *H. sapiens* | 72/72 (100%) | 705–776 | W1DWU7Y7014 |

All sequences verified correct. Boltz-1 YAML inputs confirmed ready to submit.

**Next steps:**
- [x] Submit 4 Boltz-1 jobs on ColabFold (5 diffusion samples each)
- [x] Assess complex ipTM — confirm both scaffolds form PPI with RBX1

---

### Entry 004 — 2026-03-23

**Status:** Complete

**Work completed:**
- Ran all 4 Boltz-2 jobs on ColabFold A100 (5 diffusion samples each, `--no_kernels`)
- Fixed `cuequivariance_ops_torch` dependency error by using `--no_kernels` flag
- Fixed incorrect `--num_diffusion_samples` flag → correct flag is `--diffusion_samples`
- Analysed confidence metrics across all 20 models (4 jobs × 5 samples)

**Boltz-2 results summary:**

| Scaffold | Monomer pTM (mean ± sd) | Monomer pLDDT (mean ± sd) | Complex ipTM (mean ± sd) | Complex ipTM range | Verdict |
|----------|------------------------|--------------------------|--------------------------|-------------------|---------|
| GLMN 336–582 | 0.877 ± 0.002 | 0.885 ± 0.002 | **0.847 ± 0.020** | 0.814 – 0.870 | **Green light** |
| CUL1 WHB 705–776 | 0.925 ± 0.003 | 0.939 ± 0.002 | 0.534 ± 0.222 | 0.286 – 0.845 | Proceed with caution |

**Per-model complex ipTM:**

| Model | GLMN ipTM | CUL1 WHB ipTM |
|-------|-----------|---------------|
| 0 | 0.848 | 0.734 |
| 1 | 0.861 | 0.845 |
| 2 | 0.840 | 0.341 |
| 3 | 0.814 | 0.286 |
| 4 | **0.870** | 0.462 |

**Interpretation:**
- **GLMN fragment:** All 5 models predict a confident, consistent RBX1 interface (ipTM 0.814–0.870). The isolated Glomulin fragment retains RBX1 binding without Cullin context. Strong scaffold — use model 4 (highest ipTM 0.870) for ProteinMPNN.
- **CUL1 WHB fragment:** Monomer folds exceptionally well (pTM 0.925) but complex ipTM is highly variable (0.286–0.845). Only 2/5 models form a confident interface, suggesting the 72 AA fragment needs the full Cullin scaffold for stable RBX1 docking. Use model 1 (ipTM 0.845) as lower-confidence scaffold.

**Best scaffolds selected:**
- Strategy 2a: `glmn_monomer_model_4.pdb` → ProteinMPNN redesign (high confidence)
- Strategy 2b: `cul1_whb_complex_model_1.pdb` chain A → ProteinMPNN redesign (lower confidence)

**Structural validation in PyMOL — Boltz models aligned to native crystal structures:**

*CUL1 WHB — Boltz model vs RBX1 (2LGV) interface view:*
- 1LDJ alignment RMSD = 2.993 Å (72 atoms)
- Boltz model vs 2LGV RMSD = 4.844 Å (72 atoms)
- Boltz model correctly recapitulates the WHB–RBX1 interface geometry

![CUL1 WHB Boltz model (orange/magenta) aligned to 2LGV RBX1 (green) — interface view](Boltz_CUL1WHB_RBX1_interface_alignment.png)

*CUL1 WHB — Boltz model in context of full native 1LDJ structure:*
- Boltz WHB domain vs native 1LDJ WHB region: RMSD = 0.649 Å (72 atoms) — essentially identical
- Compact coloured region (green/orange/teal) sits correctly within full Cullin scaffold (wheat)

![CUL1 WHB Boltz model (coloured) shown in context of full 1LDJ native structure (wheat)](Boltz_CUL1WHB_1LDJ_native_context.png)

*GLMN — Boltz model vs native 4F52 crystal structure, gaps filled:*
- MatchAlign score: 1249.5, 276 atoms aligned
- RMSD = 0.588 Å on 198 atoms — near-perfect match to crystal structure
- Missing loops (res 433–438, 533–549) now fully modelled (visible as extra density in green Boltz model vs salmon crystal)
- Object `aln_4F52_to_glmn_rbx1_complex_model_4` created

![GLMN Boltz model (green/teal) aligned to 4F52 crystal (salmon/magenta) — missing loops now filled](Boltz_GLMN_4F52_gaps_filled_alignment.png)

**Next steps:**
- [x] Download Boltz output PDBs → `boltz_outputs/`
- [ ] Run ProteinMPNN on GLMN model 4 scaffold (Strategy 2a)
- [ ] Run ProteinMPNN on CUL1 WHB model 1 scaffold (Strategy 2b)
- [ ] Set up RFdiffusion Colab run (Strategy 1)
- [ ] Validate all designs with Boltz complex prediction
- [ ] Filter for novelty (≥25% edit distance) and rank top 100

---

### Entry 005 — 2026-03-23

**Status:** Complete

**Work completed:**
- Ran ProteinMPNN on both Boltz-validated scaffolds (chain A designed, chain B RBX1 fixed as context)
- Fixed chain assignment bug (initially designed RBX1 instead of binder — corrected)
- Generated 48 sequences per scaffold × 3 temperatures (T=0.1, 0.2, 0.3) × 16 samples = **96 total binder sequences**
- Computed edit distance to original scaffold sequences as preliminary novelty check
- Created master tracking CSV (`master_sequences.csv`) with all sequence metadata and placeholder columns for future metrics
- Generated 192 Boltz validation YAMLs (96 monomer + 96 complex) in `boltz_validation/yamls/`

**ProteinMPNN settings:**
- Model: vanilla (v_48_020)
- Fixed chain: B (RBX1, 108 AA) — provides interface context
- Designed chain: A (binder scaffold)
- Temperatures: 0.1 (conservative), 0.2 (balanced), 0.3 (diverse)
- 16 sequences per temperature per scaffold
- Seed: 37

**Sequence statistics:**

| Scaffold | N | Length | MPNN score (mean) | MPNN score (range) |
|----------|---|--------|-------------------|--------------------|
| GLMN (Strategy 2a) | 48 | 247 AA | 0.9828 | 0.899–1.088 |
| CUL1_WHB (Strategy 2b) | 48 | 72 AA | 1.0266 | 0.948–1.163 |

**Preliminary novelty check (edit distance to original scaffold):**

| Scaffold | Temp | Edit dist range | Mean | Passes >25% |
|----------|------|-----------------|------|-------------|
| GLMN | 0.1 | 0.575–0.615 | 0.599 | 16/16 |
| GLMN | 0.2 | 0.571–0.660 | 0.600 | 16/16 |
| GLMN | 0.3 | 0.587–0.627 | 0.607 | 16/16 |
| CUL1_WHB | 0.1 | 0.556–0.639 | 0.612 | 16/16 |
| CUL1_WHB | 0.2 | 0.556–0.653 | 0.605 | 16/16 |
| CUL1_WHB | 0.3 | 0.556–0.681 | 0.632 | 16/16 |

All 96 sequences are ≥55% edit distance from the original scaffold (well above the 25% threshold). Full UniRef50 screening via mmseqs2 required before final submission — recorded as `edit_distance_uniprot50 = pending` in master CSV.

**Master CSV columns:** `seq_id`, `scaffold`, `source_pdb`, `strategy`, `temperature`, `sample`, `length`, `mpnn_score`, `edit_distance_to_original_scaffold`, `identity_to_original_scaffold`, `edit_distance_uniprot50`, `boltz_monomer_ptm`, `boltz_monomer_plddt`, `boltz_complex_iptm`, `boltz_complex_ptm`, `boltz_complex_plddt`, `passes_novelty`, `status`, `notes`, `sequence`

**Next steps:**
- [ ] Run 192 Boltz validation jobs on Colab (96 monomer + 96 complex)
- [ ] Update master CSV with Boltz metrics
- [ ] Run full UniRef50 novelty screen (mmseqs2) before submission
- [ ] Select top sequences by ipTM for final submission
- [ ] Set up RFdiffusion Colab run (Strategy 1 — de novo backbones)

---

### Entry 006 — 2026-03-23

**Status:** Complete

**Work completed:**
- Ran full Boltz-2 validation on all 96 ProteinMPNN-designed sequences
  - 96 monomer predictions (binder-only, for pTM/pLDDT)
  - 96 complex predictions (binder + RBX1, for ipTM)
  - All runs used `--use_msa_server --diffusion_samples 5 --no_kernels`
  - MSA fetched from api.colabfold.com (rate-limited, ~7–8 min per sequence)
  - Total MSA fetch time: ~107 min; GPU prediction time: ~25 min (A100)
  - 0 failed examples across all 192 predictions
- Extracted confidence scores from JSON outputs and updated `master_sequences.csv`

**Boltz-2 validation results — all 96 sequences:**

| Metric | GLMN (n=48) | CUL1_WHB (n=48) | Overall (n=96) |
|--------|-------------|-----------------|----------------|
| ipTM mean ± sd | **0.867 ± 0.011** | 0.595 ± 0.096 | 0.731 ± 0.121 |
| ipTM median | **0.865** | 0.606 | 0.759 |
| ipTM range | 0.845 – 0.887 | 0.344 – 0.692 | 0.344 – 0.887 |
| Sequences ≥ ipTM 0.7 | **48/48 (100%)** | 7/48 (15%) | 55/96 (57%) |
| Monomer pTM mean | 0.891 | 0.941 | 0.916 |
| Monomer pLDDT mean | 0.905 | 0.955 | 0.930 |

**Key findings:**
- **GLMN scaffold is a clear winner**: 100% of 48 sequences exceed ipTM 0.7; tight distribution (sd = 0.011) indicates the scaffold geometry is robustly maintained across all ProteinMPNN temperatures and seeds
- **CUL1_WHB has excellent fold quality** (pTM ~0.94, pLDDT ~0.95) but poor complex binding (only 7/48 above 0.7), consistent with Entry 004 findings that the 72 AA fragment struggles to form a stable interface without the full Cullin scaffold
- **Temperature effect (GLMN):** Minimal variance in ipTM across T=0.1/0.2/0.3 — all temperatures yield viable binders
- **Temperature effect (CUL1_WHB):** No consistent improvement at any temperature; variability is structural, not sequence-sampling related

**Top 10 sequences by ipTM:**

| Rank | Sequence ID | Scaffold | ipTM | pTM | pLDDT |
|------|-------------|----------|------|-----|-------|
| 1 | GLMN_T0.1_s11 | GLMN | 0.8869 | 0.8895 | 0.9091 |
| 2 | GLMN_T0.1_s4 | GLMN | 0.8859 | 0.8907 | 0.9054 |
| 3 | GLMN_T0.3_s12 | GLMN | 0.8816 | 0.8795 | 0.9037 |
| 4 | GLMN_T0.1_s15 | GLMN | 0.8807 | 0.8946 | 0.8983 |
| 5 | GLMN_T0.2_s5 | GLMN | 0.8794 | 0.8911 | 0.8892 |
| 6 | GLMN_T0.2_s14 | GLMN | 0.8790 | 0.8906 | 0.9011 |
| 7 | GLMN_T0.1_s13 | GLMN | 0.8784 | 0.8885 | 0.9027 |
| 8 | GLMN_T0.3_s4 | GLMN | 0.8781 | 0.8880 | 0.9024 |
| 9 | GLMN_T0.3_s8 | GLMN | 0.8781 | 0.8962 | 0.9036 |
| 10 | GLMN_T0.1_s10 | GLMN | 0.8771 | 0.8728 | 0.8765 |

**Visualisations:**

*Figure 1 — Overview: ipTM distribution, temperature breakdown, fold quality vs binding*

![Boltz-2 validation overview: violin plots, temperature breakdown, pTM vs ipTM scatter](boltz_results_overview.png)

*Figure 2 — Top 20 sequences ranked by ipTM*

![Top 20 sequences by Boltz-2 complex ipTM, coloured by scaffold](boltz_top20_ranking.png)

*Figure 3 — ipTM heatmap across all temperature × sample combinations*

![ipTM heatmap: rows = temperature (0.1/0.2/0.3), columns = samples s1–s16, per scaffold](boltz_iptm_heatmap.png)

**Submission strategy:**
- Submit all 48 GLMN sequences (ipTM 0.845–0.887, all above threshold)
- Submit the 7 CUL1_WHB sequences with ipTM ≥ 0.70 to maximise scaffold diversity
- Total planned: 55 sequences (within the 100-sequence limit, leaving headroom for Strategy 1 de novo designs if time permits)
- Pending: UniRef50 novelty screen (mmseqs2) to confirm edit_distance_uniprot50 ≥ 0.25

**Next steps:**
- [ ] Run mmseqs2 UniRef50 novelty screen on all 55 candidate sequences
- [ ] Final ranking and selection (top 100 by ipTM, filtered by novelty)
- [ ] Submit sequences by March 26, 2026 deadline

---

### Entry 007 — 2026-03-23

**Status:** Complete

**Work completed:**
- Computed RBX1 RING domain RMSD for all 96 complex predictions vs native scaffold reference
- Reference: chain B (RBX1) from `glmn_rbx1_complex_model_4.pdb` (native Glomulin–RBX1 scaffold prediction)
- Metric: Cα RMSD restricted to **residues 32–108** (structured RING-H2 domain only; residues 1–31 are disordered and inflate RMSD artificially)
- For each sequence, used the **best single-model** complex (highest single-model ipTM from 5 diffusion samples)
- Question: do our binders bind to a pre-formed RBX1 surface (**locked-in**) or force RBX1 to rearrange (**induced fit**)?

**Results — RBX1 RING RMSD (res 32–108, Cα):**

| Scaffold | n | RMSD mean ± sd | RMSD range | Locked <1Å | Moderate 1–2Å | Induced fit >2Å |
|----------|---|----------------|------------|------------|----------------|-----------------|
| **GLMN** | 48 | **1.09 ± 0.27 Å** | 0.60–1.84 Å | 18 (38%) | 30 (62%) | **0 (0%)** |
| CUL1_WHB | 48 | 5.11 ± 1.05 Å | 2.91–7.39 Å | 0 (0%) | 0 (0%) | **48 (100%)** |

**Correlation (avg ipTM vs RING RMSD, all 96 sequences): r = −0.84**

**Interpretation:**

- **GLMN — locked-in binders:** All 48 sequences have RING RMSD < 2 Å; 18 are below 1 Å. RBX1 barely moves. The ProteinMPNN-redesigned sequences dock onto RBX1's pre-formed surface with minimal receptor reorganisation. This is the ideal binding mode — high ipTM without distorting the target. The strong negative correlation (r = −0.84) confirms: **less RBX1 distortion = better predicted binding**.

- **CUL1_WHB — induced fit binders:** All 48 sequences distort the RING domain by >2 Å (mean 5.1 Å), even in the best single-model predictions. The 72 AA fragment is too short to fill the RBX1 interface without pulling loops out of position. This explains both the low average ipTM (mean 0.595) and the extreme variance across diffusion samples — Boltz is sampling different induced-fit poses rather than converging on a single bound state.

- **Biological significance:** The CUL1 WHB domain naturally contacts RBX1 as part of the ~90 kDa Cullin scaffold. The Cullin framework constrains the WHB–RBX1 geometry from both sides. As a 72 AA isolated fragment, it loses that constraint and must deform RBX1 to bind. The GLMN scaffold, by contrast, was evolved as a standalone inhibitor and maintains the correct geometry independently.

**Top 5 GLMN by RMSD (most locked-in interfaces):**

| Rank | Seq ID | avg ipTM | RING RMSD | Classification |
|------|--------|----------|-----------|----------------|
| 1 | GLMN_T0.3_s8 | 0.878 | 0.603 Å | Rigid — locked in |
| 2 | GLMN_T0.3_s12 | 0.882 | 0.910 Å | Minor flex |
| 3 | GLMN_T0.2_s5 | 0.879 | 0.892 Å | Minor flex |
| 4 | GLMN_T0.2_s14 | 0.879 | 0.915 Å | Minor flex |
| 5 | GLMN_T0.1_s11 | 0.887 | 0.960 Å | Minor flex |

**Visualisations:**

*Figure 1 — Full RMSD analysis (4-panel): distribution, ipTM vs RMSD scatter, histogram, temperature breakdown*

![RBX1 RING RMSD analysis: violin, scatter, histogram, and temperature-grouped bar plots for all 96 sequences](rbx1_rmsd_analysis.png)

*Figure 2 — Top 15 sequences: ipTM and RING RMSD side-by-side*

![Top 15 sequences by avg ipTM showing paired ipTM and RING RMSD bars coloured by scaffold](rbx1_rmsd_top15.png)

**Key design implication:**
The GLMN scaffold is not only producing better ipTM scores — it is doing so through a fundamentally superior binding mechanism. For experimental validation, locked-in binders are far more likely to succeed: they don't require the target to pay an entropic/enthalpic cost of rearrangement. For submission, **GLMN sequences with both high ipTM and low RING RMSD are the primary candidates.**

**Next steps:**
- [ ] Run mmseqs2 UniRef50 novelty screen on all 48 GLMN sequences + 7 CUL1_WHB above avg ipTM 0.7
- [ ] Final submission ranking: prioritise by ipTM, then by RING RMSD (lower = better)
- [ ] Submit by March 26, 2026 deadline

---

### Entry 008 — 2026-03-23

**Status:** Complete  
**Goal:** Novelty screen — verify all 55 candidate sequences are genuinely de novo (not natural protein copies)

#### Method
- Tool: **DIAMOND blastp v2.1.24** (sensitive mode), local search
- Database: **UniProt/SwissProt** (574,627 sequences, Mar 2026 release)
- Downloaded: `uniprot_sprot.fasta.gz` via UniProt FTP
- Threshold: sequences with best-hit identity **< 75%** pass as novel
- Note: mmseqs2 UniRef50 database download from previous session was incomplete (index files 0 bytes); DIAMOND used as equivalent local alternative

#### Results
**All 55/55 candidate sequences passed the novelty filter.**

| Scaffold | n seqs | Mean SwissProt identity | Range | Best hit |
|----------|--------|------------------------|-------|----------|
| GLMN-based | 48 | 40.0 ± 1.8% | 34.0–43.3% | sp\|Q92990\|GLMN_HUMAN |
| CUL1_WHB-based | 7 | 42.1 ± 3.3% | 37.7–46.2% | sp\|Q13616\|CUL1_HUMAN |
| **All** | **55** | **40.3 ± 2.2%** | **34.0–46.2%** | — |

No sequence exceeds the 75% identity threshold. This confirms the ProteinMPNN redesign process generates genuinely novel sequences — not mere copies or single-point mutants of the natural scaffold proteins.

**Why the identity is so low (~40%):** Despite using GLMN and CUL1_WHB as structural scaffolds, ProteinMPNN redesigns the amino acid sequence broadly, not just the surface residues. The DIAMOND alignment covers the matched region only; across the full length, the backbone identity is even lower. This is characteristic of ProteinMPNN-generated sequences — they adopt the target fold with largely unrelated sequences compared to the natural protein.

**Per-sequence detail (top 10 by ipTM):**

| Seq ID | Best SwissProt hit | Identity | Pass |
|--------|-------------------|---------|------|
| GLMN_T0.1_s11 | GLMN_HUMAN | 38.0% | YES |
| GLMN_T0.1_s4 | GLMN_HUMAN | 40.0% | YES |
| GLMN_T0.3_s12 | GLMN_HUMAN | 38.2% | YES |
| GLMN_T0.1_s15 | GLMN_HUMAN | 38.8% | YES |
| GLMN_T0.2_s5 | GLMN_HUMAN | 38.5% | YES |
| GLMN_T0.2_s14 | GLMN_HUMAN | 42.9% | YES |
| GLMN_T0.1_s13 | GLMN_HUMAN | 39.6% | YES |
| GLMN_T0.3_s8 | GLMN_HUMAN | 37.6% | YES |
| GLMN_T0.3_s4 | GLMN_HUMAN | 35.9% | YES |
| GLMN_T0.1_s10 | GLMN_HUMAN | 41.6% | YES |

**Visualisation:**

*Figure 1 — Identity distribution by scaffold and per-sequence ranked chart*

![Novelty screen: DIAMOND SwissProt identity distribution for 55 candidate sequences](novelty_screen_results.png)

#### Output files
- `novelty_swissprot.tsv` — raw DIAMOND output (BLAST-6 format)
- `novelty_screen_results.csv` — parsed per-sequence novelty results
- `final_submission.fasta` — all 55 sequences ranked by ipTM, ready for submission

#### Submission shortlist

All 55 candidates are novel and qualified. Final ranking for submission:

**GLMN-based (48 sequences, all recommended):**
- ipTM range: 0.844–0.887
- All sequences: low RING RMSD (mean 1.09 Å, locked-in binding mode)
- All sequences: novel (<75% SwissProt identity)

**CUL1_WHB-based (7 sequences, secondary candidates):**
- ipTM range: 0.701–0.761 (avg ipTM ≥ 0.70 threshold)
- Best single-model ipTM up to 0.93
- Induced-fit binding mode (RING RMSD mean 5.11 Å)
- Novel: all <75% SwissProt identity

**Recommended submission strategy:**
1. Submit all 48 GLMN-based candidates as primary binders (high ipTM + locked-in)
2. Include top 7 CUL1_WHB candidates as scaffold diversity (different binding mode)
3. Total: 55 sequences (well within 100-sequence limit)
4. If slots remain, extend GLMN ipTM threshold to 0.65 to add up to 45 more

**Next step:** Final submission by March 26, 2026 deadline.


---

### Entry 009 — 2026-03-23

**Status:** Complete  
**Goal:** Submit Batch 1 (55 sequences) + plan Strategy 1 RFdiffusion run for Batch 2 (45 sequences)

#### Batch 1 Submission — 55 Sequences

All 55 sequences submitted to GEM x Adaptyv Bio competition portal.

**Submission composition:**
- 48 × GLMN-based binders (Strategy 2a) — primary candidates
- 7 × CUL1 WHB-based binders (Strategy 2b) — scaffold diversity

**Summary metrics:**

| Class | n | ipTM range | mean RING RMSD | SwissProt id | Length |
|-------|---|-----------|----------------|-------------|--------|
| GLMN | 48 | 0.844–0.887 | 1.09 Å | ~40% | 247 AA |
| CUL1_WHB | 7 | 0.701–0.761 | 5.5 Å | ~42% | 72 AA |

**Files:**
- `final_submission.fasta` — all 55 ranked sequences
- `submission_writeup.md` — full write-up with methods, results, rationale

#### Strategy 1 — RFdiffusion De Novo Setup

Prepared all inputs for the RFdiffusion de novo backbone generation run:

**Input structure:** `rbx1_ring_renumbered.pdb`
- Source: chain B from `boltz_NATIVES_outputs/glmn_rbx1_complex_model_4.pdb`
- Chain renamed A, residues renumbered 1–77 (original: RBX1 32–108)
- 77 CA atoms (structured RING-H2 domain only)

**Hotspot residues identified (24 total, 14 core):**
- Full set: 24 RBX1 RING residues within 5 Å of GLMN in native complex
- Core set (for conditioning): `A4,A6,A12,A14,A15,A23,A24,A26,A28,A52,A56,A60,A64,A66`
  (original RBX1 numbering: W35, I37, A43, C45, R46, I54, E55, Q57, N59, C83, W87, R91, P95, D97)
- Key anchors: W35 and W87 (large hydrophobic), R46 and R91 (charged), C45/C83 (Zn-coordinating)

**RFdiffusion parameters:**
```
contigs = "A1-77/0 60-100"        # fix RBX1, generate binder of 60-100 AA
hotspot_res = "A4,A6,A12,A14,A15,A23,A24,A26,A28,A52,A56,A60,A64,A66"
num_designs = 200
```

**Plan document:** `rfdiffusion_setup.md` — full pipeline instructions for Colab run

**Next steps:**
- [ ] Upload `rbx1_ring_renumbered.pdb` to Colab and run RFdiffusion (200 designs)
- [ ] Filter by binder pLDDT > 0.80
- [ ] Run ProteinMPNN (3 temps × ~100 backbones = ~300 sequences)
- [ ] Boltz-2 monomer + complex validation
- [ ] Select top 45 by ipTM → submit as Batch 2


---

### Entry 010 — 2026-03-24

**Status:** Complete  
**Goal:** PyMOL structural validation — document scaffold extraction and Boltz-2 prediction quality against native crystal structures

#### Scaffold Extraction (PyMOL Sessions)

Two scaffolds were extracted from source PDBs using PyMOL prior to ProteinMPNN redesign.

**Scaffold 2a — GLMN from 4F52**

![PyMOL session: PDB 4F52 (green, full GLMN-RBX1-CUL1 complex) with the GLMN contact region (residues 336-582, orange highlight) selected for extraction as Scaffold_4F52_336-582.pdb](Scaffold_4F52_582-336.png)

PDB 4F52 (Duda et al. 2012) loaded into PyMOL. The GLMN chain residues 336–582 (247 AA) were selected by clicking residue LEU 336 and running `set_name sele, Scaffold_4F52_336-582`. The orange highlighted region shows the selected fragment embedded in the full complex. The 4F52 structure aligns to 2LGV (RBX1 free NMR structure) at RMSD = 7.692 Å over 88 atoms — the expected deviation between the complex form and the free RING domain.

**Scaffold 2b — CUL1 WHB domain from 1LDJ**

![PyMOL session: 1LDJ loaded with CUL1 WHB domain (residues 705-776, red/orange) selected and saved as Scaffold_1LDJ_705-776.pdb](Scaffold_1LDJ_705-776.png)

PDB 1LDJ (Zheng et al. 2002) loaded. The CUL1 WHB domain (residues 705–776, 72 AA) was selected (GLU 705 → end) and saved as `Scaffold_1LDJ_705-776.pdb`. The isolated WHB domain (red/orange compact structure) is visibly small relative to the full Cullin scaffold it was extracted from — foreshadowing the induced-fit binding problem observed later.

---

#### Boltz-2 Prediction Quality — GLMN Scaffold

**Global alignment: Boltz-2 GLMN-RBX1 vs native 4F52 crystal structure**

![PyMOL MatchAlign of Boltz-2 GLMN-RBX1 prediction (green/magenta) against native 4F52 crystal structure (cyan/salmon). RMSD = 0.588 Å over 198 Cα atoms after outlier rejection.](Boltz_GLMN_4F52_gaps_filled_alignment.png)

MatchAlign (iterative outlier rejection) of `glmn_rbx1_complex_model_4` against native 4F52. After 5 rejection cycles:

| Cycle | Atoms rejected | RMSD |
|-------|---------------|------|
| 1 | 15 | 26.35 Å |
| 2 | 2 | 3.05 Å |
| 3 | 13 | 1.03 Å |
| 4 | 11 | 0.73 Å |
| 5 | 6 | **0.62 Å** |
| **Final** | — | **0.588 Å (198 atoms)** |

**RMSD = 0.588 Å over 198 Cα atoms** — essentially perfect agreement with the X-ray crystal structure. This is exceptional: Boltz-2, given only the GLMN sequence (no template), predicted a structure that deviates from the published 2.1 Å resolution crystal structure by less than one hydrogen bond length. This validates:
1. The GLMN scaffold folds correctly as predicted
2. The binding geometry to RBX1 is faithfully captured
3. ProteinMPNN sequence designs on this scaffold are grounded in an experimentally validated structural context

---

#### Boltz-2 Prediction Quality — CUL1 WHB Scaffold

**CUL1 WHB complex vs 1LDJ native context**

![PyMOL session: CUL1_WHB_rbx1_complex_model_1 (orange/green) aligned to the full 1LDJ SCF complex (beige ribbon). RMSD to 1LDJ = 4.844 Å (72 atoms).](Boltz_CUL1WHB_1LDJ_native_context.png)

The Boltz-2 CUL1_WHB-RBX1 complex prediction was aligned to the native 1LDJ SCF complex. Result: **RMSD = 4.844 Å** (72 atoms) — a large deviation compared to the 0.588 Å achieved by the GLMN prediction. The WHB domain (orange/green) is visibly displaced relative to where it sits in the intact Cullin-RBX1 geometry (beige ribbon background). This displacement arises because the isolated 72 AA fragment lacks the scaffolding from the rest of Cullin-1 that normally constrains the WHB-RBX1 geometry.

**CUL1 WHB interface view**

![PyMOL interface view: CUL1_WHB binder chain (magenta) and RBX1 (green) from Boltz-2 prediction, overlaid with native 1LDJ WHB (orange). Interface contact region shows the induced-fit displacement.](Boltz_CUL1WHB_RBX1_interface_alignment.png)

Rotated view of the same alignment focused on the RBX1 binding interface. Magenta = designed WHB binder chain, green = RBX1 (from Boltz-2 prediction), orange = native 1LDJ WHB domain at its canonical position in the SCF complex. The offset between magenta and orange makes the induced-fit rearrangement visually clear: the binder contacts RBX1 from a shifted angle, requiring RBX1 to accommodate the new geometry (explaining the high RING RMSD of 4.844 Å seen in the RMSD analysis). The 1LDJ alignment RMSD of 2.993 Å vs the complex model's 4.844 Å shows that even the rigid-body alignment cannot fully reconcile the two poses.

---

#### Best Individual Model Analysis

**Best CUL1_WHB design: CUL1_WHB_T0.1_s14 — model_2**

![PyMOL: CUL1_WHB_T0.1_s14 five diffusion models (0-4) aligned to 1LDJ. Model_2 achieves RMSD = 0.675 Å — perfect native geometry recovery. Other models: 3.968, 3.897, 6.284, 6.344 Å.](best_1dhj\ result.png)

`CUL1_WHB_T0.1_s14` aligned all 5 diffusion samples to the native 1LDJ reference:

| Model | RMSD to 1LDJ |
|-------|-------------|
| model_0 | 6.284 Å |
| model_1 | 3.968 Å |
| **model_2** | **0.675 Å** |
| model_3 | 3.897 Å |
| model_4 | 6.344 Å |

**Model_2 achieves 0.675 Å RMSD** — it has perfectly recovered the native CUL1 WHB-RBX1 binding geometry from first principles, using only the WHB amino acid sequence. The four other diffusion samples land at 3.9–6.3 Å — entirely different poses. This bimodal distribution (one near-native sample, four diverged samples) mechanistically explains:
- Why the average ipTM across models is moderate (0.70–0.76 for passing sequences)
- Why the best single-model ipTM reaches 0.93
- Why RING RMSD variance is so high for CUL1_WHB sequences

The WHB domain *can* adopt the correct binding mode — but without the Cullin scaffold to provide context, Boltz-2 only finds it in ~1 of 5 diffusion samples. This is a genuine induced-fit system.

**Best GLMN design: GLMN_T0.1_s8 — all 5 models**

![PyMOL: GLMN_T0.1_s8 five diffusion models (magenta/green) aligned to 4F52 native structure (grey). All five models cluster at 1.894–2.188 Å RMSD over 224 atoms — locked-in convergence.](best_design_4F52_pymol_session.png)

`GLMN_T0.1_s8` aligned all 5 diffusion samples to 4F52:

| Model | RMSD to 4F52 |
|-------|-------------|
| model_0 | 2.038 Å |
| model_1 | 2.188 Å |
| model_2 | 2.135 Å |
| model_3 | 2.055 Å |
| **model_4** | **1.894 Å** |

All five models fall within a **0.294 Å window** (1.894–2.188 Å). This tight clustering across independent diffusion samples is the hallmark of a locked-in binding mode: Boltz-2 has no uncertainty about the bound pose and converges to the same solution regardless of the random diffusion trajectory. The 2.0 Å offset from the 4F52 crystal structure (vs 0.588 Å for the native GLMN sequence) is expected — the ProteinMPNN-redesigned sequence introduces mutations that shift the side-chain packing slightly while maintaining the overall backbone geometry.

For comparison, the GLMN_T0.1_s8 avg ipTM is 0.851 and RING RMSD is 1.21 Å — solidly in the locked-in regime.

---

#### Summary: What the PyMOL Sessions Confirm

| Finding | Evidence |
|---------|---------|
| GLMN scaffold prediction quality is excellent | 0.588 Å RMSD to 4F52 crystal structure |
| GLMN designs converge consistently | All 5 diffusion samples within 0.3 Å of each other |
| CUL1_WHB is a genuine induced-fit system | 4.844 Å avg RMSD to 1LDJ; bimodal per-model distribution |
| CUL1_WHB CAN recover native geometry | Model_2 of s14 = 0.675 Å to native |
| Design rationale is structurally validated | GLMN predictions match X-ray to sub-Ångstrom accuracy |


---

### Entry 011 — 2026-03-24

**Status:** Complete  
**Goal:** Visualise RFdiffusion hotspot residues and contig in PyMOL — confirm the conditioning surface is coherent before launching the Colab run

#### PyMOL Visualisation of RFdiffusion Inputs

![PyMOL: rbx1_ring_renumbered.pdb (blue surface/cartoon, chain A residues 1-77) with core hotspot residues shown as red sticks. Pink helical structure and green chain in background are native complex objects from prior session.](hotspots_for_rfdiffusion.png)

**Input file:** `rbx1_ring_renumbered.pdb` (chain A, res 1–77; original RBX1 numbering 32–108)

**What is shown:**
- Blue surface + cartoon = RBX1 RING domain (the fixed contig for RFdiffusion: `A1-77`)
- Red sticks = 14 core hotspot residues (`A4,A6,A12,A14,A15,A23,A24,A26,A28,A52,A56,A60,A64,A66`)
- Pink helical + green objects in background = GLMN and RBX1 from the native complex (previous session objects, not part of the RFdiffusion input)

**Key observation — hotspots form a coherent concave patch:**

The 14 core hotspot residues are not scattered across the RBX1 surface — they cluster tightly on a single concave face of the RING domain. This is critical for RFdiffusion: conditioning on a well-defined surface patch lets the diffusion model converge on a coherent binder backbone that approaches from a consistent direction. Dispersed hotspots would produce backbones approaching from multiple angles that cannot all be satisfied simultaneously.

The concave geometry of the hotspot pocket (visible in the surface representation) is also favourable: it provides shape complementarity that a designed helix or loop can dock into, giving enthalpic driving force for binding beyond just hydrophobic burial.

**Validation of hotspot selection:**
The hotspot face in the renumbered structure aligns with the direction from which the GLMN backbone (pink, background) approaches in the native complex — confirming that the conditioning surface correctly identifies the biologically relevant binding face. RFdiffusion will explore new backbone geometries that contact this same face, potentially from different angles or with different secondary structure compositions than GLMN.

**Ready to run RFdiffusion.** Upload `rbx1_ring_renumbered.pdb` to Colab and use:
```
contigs = "A1-77/0 60-100"
hotspot_res = "A4,A6,A12,A14,A15,A23,A24,A26,A28,A52,A56,A60,A64,A66"
num_designs = 200
```


---

### Entry 012 — 2026-03-25

**Status:** Complete
**Goal:** Analyse Batch 2 RFdiffusion results, run Nipah competition retrospective, apply evidence-based metrics to re-score all candidates, finalise top 100 for submission.

---

#### Batch 2 Boltz-2 Results — RFdiffusion Sequences

Boltz-2 complex predictions completed on Colab A100 after ~2.5 hours (including one restart due to missing `cuequivariance-torch` dependency). Results downloaded to `boltz_rfdiffusion/`.

**Funnel: 200 → 31 (original) → 57 (after re-scoring):**

| Stage | n | Loss | Reason |
|---|---|---|---|
| RFdiffusion backbones | 200 | — | — |
| Length filter (65–95 AA) | 151 | −49 | Binder length outside target |
| ProteinMPNN (best-per-backbone) | 151 | 0 | 1 sequence per backbone |
| Boltz-2 monomer completed | 118 | −33 | Monomer run stopped early |
| Pass ipTM ≥ 0.70 (complex) | 57 | −61 | Failed complex confidence filter |
| Clean (no poly-Ala termini) | 56 | −1 | RFD_27 had 11+17 terminal alanines |

The monomer pTM filter (originally applied) eliminated a further 25 sequences from the 57 — but Nipah retrospective showed pTM has AUROC 0.501 (random), so this filter was dropped (see below).

**Score ranges (Batch 2, n=151 with complex data):**
- ipTM: 0.33–0.91 (mean 0.70), 57/151 pass ≥ 0.70
- complex ipLDDT: 0.55–0.85 (mean 0.68)

![Batch 2 funnel analysis and Batch 1 vs Batch 2 score comparison](batch2_analysis.png)

---

#### Nipah Competition Retrospective Analysis

Analysed 1,030 sequences from the Adaptyv Nipah binder competition with full experimental outcomes (BLI/SPR binding data).

**Key findings:**

| Finding | Value |
|---|---|
| Overall binding rate | 10% (103/1030) |
| Expression rate | 86% (884/1030) |
| Best-performing class | scFv (31% binding rate) |
| Median KD of binders | 24.5 nM |
| Sub-nM binders | 8/102 with KD data |

**Predictor AUROC ranking (binders vs non-binders):**

| Metric | AUROC | Significance | Used previously? |
|---|---|---|---|
| Boltz2 complex ipLDDT | **0.691** | p=2.5e-10 *** | No |
| Shape complementarity | **0.687** | p=6.8e-10 *** | No |
| Boltz2 pLDDT (monomer) | 0.640 | p=3.5e-6 *** | No |
| Boltz2 min ipSAE | 0.638 | p=4.8e-6 *** | No |
| Boltz2 ipTM | 0.603 | p=6.9e-4 *** | **Yes (only one)** |
| **Boltz2 pTM (monomer)** | **0.501** | p=0.97 (ns) | **Yes (as filter!)** |

**Critical insight:** The monomer pTM filter we applied (AUROC 0.501) was eliminating sequences at random — it has zero predictive power for experimental binding. `complex_iplddt` is a 15-point better predictor than ipTM. Both are better ranking criteria than pTM.

**Threshold calibration note:** The Nipah-optimal ipLDDT threshold (0.850) is not portable to our GLMN/RFdiffusion sequences — our ipLDDT values range 0.61–0.85 (median 0.71) versus Nipah scFv/nanobodies which scored 0.80–0.95. The Nipah result tells us *ipLDDT ranks binders better*, not that an absolute threshold of 0.85 applies cross-system.

![Nipah retrospective analysis — predictor AUROC, KD distribution, violin plots](nipah_analysis/nipah_analysis.png)

---

#### Re-scoring with Nipah-Derived Approach

Applied evidence-based scoring changes to all 247 candidates (96 Batch 1 + 151 Batch 2):

**New filter and ranking:**
- **Keep:** ipTM ≥ 0.70 (hard gate — AUROC 0.603, still meaningful)
- **Drop:** mono pTM ≥ 0.70 (AUROC 0.501 — random, costs us sequences with no benefit)
- **Drop:** poly-Ala terminal sequences (RFD_27)
- **Rank by composite score:** 0.5 × ipTM + 0.5 × complex_ipLDDT (incorporates Nipah's best predictor)

**Results:**

| | Old filter | New filter |
|---|---|---|
| Batch 1 (GLMN + CUL1_WHB) | 55 | 55 |
| Batch 2 (RFdiffusion) | 32 | 57 |
| **Total** | **87** | **112** |

The 26 sequences gained are Batch 2 designs whose monomers failed or were missing — correctly recovered by dropping the uninformative pTM filter. One sequence lost (RFD_27, poly-Ala termini).

112 > 100 cap → **top 100 by composite score selected for final submission: 53 Batch 1 + 47 Batch 2.**

![Re-scoring with Nipah-derived thresholds](rescore_analysis.png)

**Interim submission file:** `rescored_all.csv` — 112 passing sequences with composite2 score (0.5×ipTM + 0.5×ipLDDT). Final ranking subsequently updated to composite3 after ipSAE computation (see Entry 013).

![Re-scoring with Nipah-derived thresholds — old vs new filter, ipTM vs ipLDDT scatter](rescore_analysis.png)

---

#### Files generated this session

| File | Description |
|---|---|
| `batch2_analysis.py` / `batch2_analysis.png` | Batch 2 funnel + Batch 1 vs Batch 2 comparison |
| `nipah_analysis/analyze_nipah.py` | Nipah retrospective analysis script |
| `nipah_analysis/nipah_analysis.png` | Nipah predictor AUROC, binding rate, KD distribution |
| `rescore_candidates.py` | Re-scoring script with Nipah-derived composite |
| `rescore_analysis.png` | Old vs new filter comparison visualisation |
| `rescored_all.csv` | All 112 passing sequences with scores |

---

---

### Entry 013 — 2026-03-25

**Status:** Complete
**Goal:** Compute ipSAE for all designs, integrate into final composite score, run Batch 2 novelty screen, finalise submission FASTA.

---

#### ipSAE Computation

Computed interface Predicted Structural Assessment Error (ipSAE) from Boltz-2 PAE matrices for all 247 sequences (96 Batch 1 + 151 Batch 2).

**What ipSAE measures:** Mean PAE (Predicted Aligned Error) of the cross-chain residue block in the Boltz-2 PAE matrix. Entry [i,j] = expected positional error of residue j when structure is aligned on residue i. The cross-chain sub-matrix captures how well Boltz-2 knows the relative positions of binder and RBX1 residues at the interface. Lower = more confident = better. Raw units: Angstroms (0–30 Å range).

**Results by scaffold:**

| Scaffold | n | mean ipSAE | range | Interpretation |
|---|---|---|---|---|
| GLMN | 48 | 9.06 Å | 7.72–10.79 | Moderate, consistent interface confidence |
| CUL1_WHB | 48 | 12.47 Å | 8.82–17.66 | Poor — confirms induced-fit geometry instability |
| RFdiffusion | 151 | 11.31 Å | 4.83–21.46 | Wide range — best and worst both here |

**Key finding:** RFD_167 (5.56 Å) and RFD_114 (6.0 Å) have lower ipSAE than any GLMN sequence. The de novo approach, when it works, produces more precisely defined interfaces than the scaffold-based approach.

**Final composite score (3-metric):**
```
composite3 = 0.4 × ipTM + 0.3 × ipLDDT + 0.3 × norm_ipSAE
```
where `norm_ipSAE = 1 - (ipSAE - min) / (max - min)` (inverted so higher = better).

![ipSAE distribution by scaffold, ipSAE vs ipTM scatter, rank shifts, final top 20](ipsae_analysis.png)

---

#### Batch 2 Novelty Screen

DIAMOND blastp vs SwissProt (574,627 sequences) on the 47 Batch 2 sequences in the final top 100:

- **Result: 0/47 hits at e-value 1e-3 — completely novel**
- Mean SwissProt identity: 0% (not detectable)
- All 47 sequences pass the <75% identity novelty threshold by a vast margin
- Stronger novelty than Batch 1 GLMN (40% identity to natural Glomulin)

RFdiffusion generates sequences with no detectable homology to any known protein — genuinely *de novo* in every sense.

---

#### Final Submission

**File:** `final_submission_v2.fasta` — **100 sequences**

| Batch | Scaffold | Sequences | ipTM range | Novelty |
|---|---|---|---|---|
| Batch 1 | GLMN | 46 | 0.844–0.887 | ~40% SwissProt identity (all pass) |
| Batch 1 | CUL1_WHB | 7 | 0.701–0.761 | ~42% SwissProt identity (all pass) |
| Batch 2 | RFdiffusion | 47 | 0.710–0.910 | 0% (no hits) |
| **Total** | | **100** | **0.701–0.910** | **All pass** |

Ranked by composite3 score. Top 5: RFD_167, RFD_114, RFD_199, RFD_97, RFD_1.

---

#### Files generated this session

| File | Description |
|---|---|
| `compute_ipsae.py` | Extracts cross-chain PAE blocks, computes avg_ipsae + min_ipsae |
| `ipsae_results.csv` | ipSAE values for all 247 sequences |
| `ipsae_analysis.png` | ipSAE distribution, scatter, rank shifts, final top 20 |
| `rescored_all.csv` | Updated — all 112 passing sequences with composite3 |
| `batch2_novelty_screen.fasta` | 47 Batch 2 sequences for DIAMOND screen |
| `batch2_novelty_results.csv` | Novelty screen results (all 47: NO_HIT) |
| `final_submission_v2.fasta` | **Final 100 sequences for submission** |

---

## Sequences Submitted

**Total: 100 | File: `final_submission_v2.fasta` | Ranked by composite3 = 0.4×ipTM + 0.3×ipLDDT + 0.3×norm_ipSAE**

| Rank | Seq ID | Batch | Scaffold | Length (AA) | ipTM | ipLDDT | ipSAE (Å) | Composite |
|---|---|---|---|---|---|---|---|---|
| 1 | RFD_167_best | Batch2 | RFdiffusion | 70 | 0.870 | 0.656 | 5.56 | 0.8451 |
| 2 | RFD_114_best | Batch2 | RFdiffusion | 85 | 0.864 | 0.724 | 5.99 | 0.8426 |
| 3 | RFD_199_best | Batch2 | RFdiffusion | 80 | 0.873 | 0.622 | 5.69 | 0.8299 |
| 4 | RFD_97_best | Batch2 | RFdiffusion | 70 | 0.872 | 0.711 | 6.37 | 0.8236 |
| 5 | RFD_1_best | Batch2 | RFdiffusion | 75 | 0.854 | 0.677 | 6.12 | 0.8180 |
| 6 | RFD_52_best | Batch2 | RFdiffusion | 75 | 0.826 | 0.783 | 6.63 | 0.8152 |
| 7 | RFD_38_best | Batch2 | RFdiffusion | 94 | 0.870 | 0.816 | 7.29 | 0.8113 |
| 8 | RFD_162_best | Batch2 | RFdiffusion | 89 | 0.846 | 0.639 | 6.12 | 0.8039 |
| 9 | RFD_34_best | Batch2 | RFdiffusion | 78 | 0.868 | 0.759 | 7.17 | 0.7990 |
| 10 | RFD_106_best | Batch2 | RFdiffusion | 93 | 0.844 | 0.654 | 6.49 | 0.7904 |
| 11 | RFD_118_best | Batch2 | RFdiffusion | 93 | 0.815 | 0.705 | 6.64 | 0.7867 |
| 12 | RFD_152_best | Batch2 | RFdiffusion | 79 | 0.817 | 0.796 | 7.24 | 0.7865 |
| 13 | RFD_122_best | Batch2 | RFdiffusion | 92 | 0.882 | 0.703 | 7.23 | 0.7851 |
| 14 | RFD_173_best | Batch2 | RFdiffusion | 89 | 0.802 | 0.722 | 6.76 | 0.7810 |
| 15 | RFD_160_best | Batch2 | RFdiffusion | 69 | 0.910 | 0.847 | 8.53 | 0.7779 |
| 16 | GLMN_T0.1_s4 | Batch1 | GLMN | 247 | 0.886 | 0.739 | 7.72 | 0.7743 |
| 17 | GLMN_T0.1_s13 | Batch1 | GLMN | 247 | 0.878 | 0.753 | 8.11 | 0.7568 |
| 18 | RFD_37_best | Batch2 | RFdiffusion | 77 | 0.859 | 0.792 | 8.26 | 0.7537 |
| 19 | RFD_10_best | Batch2 | RFdiffusion | 90 | 0.823 | 0.676 | 7.24 | 0.7531 |
| 20 | GLMN_T0.1_s11 | Batch1 | GLMN | 247 | 0.887 | 0.705 | 7.97 | 0.7523 |
| 21 | RFD_195_best | Batch2 | RFdiffusion | 93 | 0.880 | 0.716 | 8.02 | 0.7510 |
| 22 | GLMN_T0.2_s5 | Batch1 | GLMN | 247 | 0.879 | 0.744 | 8.19 | 0.7509 |
| 23 | GLMN_T0.3_s12 | Batch1 | GLMN | 247 | 0.882 | 0.706 | 7.99 | 0.7496 |
| 24 | GLMN_T0.2_s14 | Batch1 | GLMN | 247 | 0.879 | 0.716 | 8.13 | 0.7452 |
| 25 | RFD_16_best | Batch2 | RFdiffusion | 94 | 0.820 | 0.732 | 7.81 | 0.7414 |
| 26 | GLMN_T0.1_s16 | Batch1 | GLMN | 247 | 0.871 | 0.731 | 8.25 | 0.7404 |
| 27 | RFD_91_best | Batch2 | RFdiffusion | 69 | 0.838 | 0.721 | 7.97 | 0.7380 |
| 28 | GLMN_T0.3_s4 | Batch1 | GLMN | 247 | 0.878 | 0.690 | 8.12 | 0.7374 |
| 29 | GLMN_T0.3_s7 | Batch1 | GLMN | 247 | 0.877 | 0.777 | 8.71 | 0.7351 |
| 30 | GLMN_T0.2_s15 | Batch1 | GLMN | 247 | 0.876 | 0.692 | 8.23 | 0.7321 |
| 31 | GLMN_T0.1_s10 | Batch1 | GLMN | 247 | 0.877 | 0.719 | 8.48 | 0.7289 |
| 32 | GLMN_T0.2_s7 | Batch1 | GLMN | 247 | 0.866 | 0.715 | 8.36 | 0.7286 |
| 33 | GLMN_T0.3_s8 | Batch1 | GLMN | 247 | 0.878 | 0.740 | 8.66 | 0.7273 |
| 34 | GLMN_T0.1_s2 | Batch1 | GLMN | 247 | 0.877 | 0.714 | 8.57 | 0.7231 |
| 35 | GLMN_T0.1_s15 | Batch1 | GLMN | 247 | 0.881 | 0.705 | 8.56 | 0.7221 |
| 36 | RFD_56_best | Batch2 | RFdiffusion | 88 | 0.879 | 0.842 | 9.44 | 0.7208 |
| 37 | GLMN_T0.1_s5 | Batch1 | GLMN | 247 | 0.865 | 0.726 | 8.74 | 0.7138 |
| 38 | RFD_110_best | Batch2 | RFdiffusion | 72 | 0.752 | 0.727 | 7.79 | 0.7137 |
| 39 | GLMN_T0.1_s3 | Batch1 | GLMN | 247 | 0.873 | 0.758 | 9.08 | 0.7103 |
| 40 | GLMN_T0.2_s3 | Batch1 | GLMN | 247 | 0.869 | 0.707 | 8.75 | 0.7092 |
| 41 | GLMN_T0.3_s6 | Batch1 | GLMN | 247 | 0.861 | 0.748 | 9.00 | 0.7068 |
| 42 | GLMN_T0.3_s15 | Batch1 | GLMN | 247 | 0.867 | 0.730 | 8.95 | 0.7061 |
| 43 | GLMN_T0.1_s1 | Batch1 | GLMN | 247 | 0.871 | 0.698 | 8.84 | 0.7033 |
| 44 | GLMN_T0.1_s9 | Batch1 | GLMN | 247 | 0.875 | 0.724 | 9.04 | 0.7032 |
| 45 | RFD_89_best | Batch2 | RFdiffusion | 95 | 0.805 | 0.646 | 7.96 | 0.7023 |
| 46 | RFD_19_best | Batch2 | RFdiffusion | 70 | 0.751 | 0.604 | 7.33 | 0.6981 |
| 47 | GLMN_T0.1_s7 | Batch1 | GLMN | 247 | 0.870 | 0.714 | 9.05 | 0.6972 |
| 48 | GLMN_T0.2_s4 | Batch1 | GLMN | 247 | 0.859 | 0.689 | 8.81 | 0.6966 |
| 49 | GLMN_T0.2_s8 | Batch1 | GLMN | 247 | 0.864 | 0.728 | 9.11 | 0.6963 |
| 50 | GLMN_T0.2_s16 | Batch1 | GLMN | 247 | 0.858 | 0.676 | 8.78 | 0.6939 |
| 51 | GLMN_T0.3_s5 | Batch1 | GLMN | 247 | 0.872 | 0.713 | 9.18 | 0.6918 |
| 52 | RFD_47_best | Batch2 | RFdiffusion | 84 | 0.734 | 0.780 | 8.47 | 0.6901 |
| 53 | GLMN_T0.1_s14 | Batch1 | GLMN | 247 | 0.863 | 0.708 | 9.16 | 0.6874 |
| 54 | RFD_33_best | Batch2 | RFdiffusion | 79 | 0.779 | 0.721 | 8.59 | 0.6850 |
| 55 | GLMN_T0.3_s13 | Batch1 | GLMN | 247 | 0.858 | 0.758 | 9.52 | 0.6837 |
| 56 | GLMN_T0.2_s10 | Batch1 | GLMN | 247 | 0.870 | 0.713 | 9.37 | 0.6820 |
| 57 | GLMN_T0.3_s10 | Batch1 | GLMN | 247 | 0.852 | 0.716 | 9.24 | 0.6820 |
| 58 | RFD_7_best | Batch2 | RFdiffusion | 80 | 0.780 | 0.676 | 8.42 | 0.6802 |
| 59 | RFD_26_best | Batch2 | RFdiffusion | 75 | 0.744 | 0.759 | 8.66 | 0.6786 |
| 60 | GLMN_T0.2_s9 | Batch1 | GLMN | 247 | 0.864 | 0.694 | 9.29 | 0.6778 |
| 61 | GLMN_T0.3_s16 | Batch1 | GLMN | 247 | 0.861 | 0.612 | 8.74 | 0.6775 |
| 62 | GLMN_T0.2_s12 | Batch1 | GLMN | 247 | 0.854 | 0.738 | 9.50 | 0.6770 |
| 63 | GLMN_T0.2_s6 | Batch1 | GLMN | 247 | 0.875 | 0.687 | 9.39 | 0.6754 |
| 64 | GLMN_T0.3_s11 | Batch1 | GLMN | 247 | 0.844 | 0.753 | 9.59 | 0.6729 |
| 65 | RFD_5_best | Batch2 | RFdiffusion | 91 | 0.770 | 0.634 | 8.24 | 0.6715 |
| 66 | CUL1_WHB_T0.2_s16 | Batch1 | CUL1_WHB | 72 | 0.761 | 0.731 | 8.82 | 0.6699 |
| 67 | GLMN_T0.3_s14 | Batch1 | GLMN | 247 | 0.854 | 0.684 | 9.36 | 0.6680 |
| 68 | RFD_172_best | Batch2 | RFdiffusion | 95 | 0.782 | 0.831 | 9.78 | 0.6632 |
| 69 | GLMN_T0.2_s13 | Batch1 | GLMN | 247 | 0.860 | 0.695 | 9.57 | 0.6629 |
| 70 | RFD_90_best | Batch2 | RFdiffusion | 85 | 0.795 | 0.642 | 8.72 | 0.6615 |
| 71 | GLMN_T0.1_s8 | Batch1 | GLMN | 247 | 0.851 | 0.675 | 9.50 | 0.6567 |
| 72 | GLMN_T0.1_s6 | Batch1 | GLMN | 247 | 0.864 | 0.725 | 9.95 | 0.6557 |
| 73 | RFD_120_best | Batch2 | RFdiffusion | 75 | 0.774 | 0.691 | 8.99 | 0.6545 |
| 74 | GLMN_T0.2_s1 | Batch1 | GLMN | 247 | 0.859 | 0.724 | 9.94 | 0.6540 |
| 75 | RFD_185_best | Batch2 | RFdiffusion | 71 | 0.744 | 0.801 | 9.47 | 0.6533 |
| 76 | GLMN_T0.2_s2 | Batch1 | GLMN | 247 | 0.853 | 0.665 | 9.57 | 0.6513 |
| 77 | GLMN_T0.3_s2 | Batch1 | GLMN | 247 | 0.856 | 0.709 | 9.92 | 0.6492 |
| 78 | RFD_58_best | Batch2 | RFdiffusion | 84 | 0.758 | 0.754 | 9.38 | 0.6491 |
| 79 | CUL1_WHB_T0.1_s4 | Batch1 | CUL1_WHB | 72 | 0.748 | 0.690 | 8.90 | 0.6488 |
| 80 | GLMN_T0.3_s3 | Batch1 | GLMN | 247 | 0.864 | 0.731 | 10.17 | 0.6469 |
| 81 | CUL1_WHB_T0.3_s11 | Batch1 | CUL1_WHB | 72 | 0.727 | 0.704 | 8.85 | 0.6467 |
| 82 | RFD_78_best | Batch2 | RFdiffusion | 66 | 0.751 | 0.631 | 8.60 | 0.6463 |
| 83 | RFD_20_best | Batch2 | RFdiffusion | 77 | 0.794 | 0.713 | 9.54 | 0.6436 |
| 84 | RFD_2_best | Batch2 | RFdiffusion | 89 | 0.728 | 0.750 | 9.25 | 0.6416 |
| 85 | GLMN_T0.1_s12 | Batch1 | GLMN | 247 | 0.865 | 0.684 | 10.01 | 0.6410 |
| 86 | RFD_133_best | Batch2 | RFdiffusion | 93 | 0.755 | 0.700 | 9.21 | 0.6397 |
| 87 | RFD_21_best | Batch2 | RFdiffusion | 92 | 0.736 | 0.681 | 8.94 | 0.6393 |
| 88 | CUL1_WHB_T0.1_s8 | Batch1 | CUL1_WHB | 72 | 0.759 | 0.736 | 9.56 | 0.6355 |
| 89 | RFD_132_best | Batch2 | RFdiffusion | 92 | 0.791 | 0.744 | 9.89 | 0.6351 |
| 90 | RFD_190_best | Batch2 | RFdiffusion | 76 | 0.731 | 0.634 | 8.72 | 0.6336 |
| 91 | GLMN_T0.3_s9 | Batch1 | GLMN | 247 | 0.862 | 0.719 | 10.49 | 0.6279 |
| 92 | GLMN_T0.3_s1 | Batch1 | GLMN | 247 | 0.851 | 0.725 | 10.48 | 0.6254 |
| 93 | RFD_131_best | Batch2 | RFdiffusion | 81 | 0.787 | 0.830 | 10.82 | 0.6153 |
| 94 | RFD_28_best | Batch2 | RFdiffusion | 82 | 0.723 | 0.733 | 9.73 | 0.6125 |
| 95 | CUL1_WHB_T0.2_s3 | Batch1 | CUL1_WHB | 72 | 0.701 | 0.708 | 9.43 | 0.6099 |
| 96 | RFD_77_best | Batch2 | RFdiffusion | 66 | 0.732 | 0.586 | 9.00 | 0.6064 |
| 97 | GLMN_T0.2_s11 | Batch1 | GLMN | 247 | 0.846 | 0.695 | 10.79 | 0.5999 |
| 98 | RFD_45_best | Batch2 | RFdiffusion | 68 | 0.766 | 0.778 | 10.73 | 0.5952 |
| 99 | RFD_29_best | Batch2 | RFdiffusion | 95 | 0.701 | 0.697 | 9.82 | 0.5882 |
| 100 | RFD_12_best | Batch2 | RFdiffusion | 65 | 0.735 | 0.696 | 10.31 | 0.5783 |

---

## Entry 014 — ipSAE Recalibration & Design Improvement Strategy
**Date:** 2026-03-26

### Background
Proteinbase scores revealed our submitted sequences have ipSAE values of ~0.40–0.73 on their 0–1 scale (higher = better; threshold: >0.6 = likely binder, ~0.8 = strong). Our original ipSAE was raw Å (lower = better) — not directly comparable. Proteinbase re-ran Boltz-2 on the **full-length RBX1 (108 AA)**, not our 77 AA RING domain, so their PAE matrices are fundamentally different. We cannot exactly replicate their scores from our existing data.

### Updated ipSAE metric (normalised)

Computed two new internally-consistent normalised metrics using `recompute_ipsae_norm.py`:

| Metric | Formula | Range | Notes |
|---|---|---|---|
| `ipsae_tm8` | mean(1/(1+(PAE/8)²)) over cross-chain pairs | 0–1, higher=better | TM-score-like, d0=8 Å |
| `ipsae_frac5` | fraction of cross-chain pairs with PAE < 5 Å | 0–1, higher=better | Strict — only high-confidence contacts |

Updated composite: `composite4 = 0.35×ipTM + 0.30×ipLDDT + 0.35×ipsae_tm8` (ipSAE upweighted per Nipah AUROC ranking)

### Results by scaffold

| Scaffold | n | ipSAE_tm8 mean | ipSAE_tm8 range | composite4 mean |
|---|---|---|---|---|
| GLMN | 48 | 0.571 | 0.513–0.626 | 0.717 |
| RFdiffusion | 57 | 0.568 | 0.419–0.726 | 0.691 |
| CUL1_WHB | 7 | 0.491 | 0.430–0.536 | 0.642 |

GLMN has the tightest, most consistent ipSAE — natural scaffold advantage. RFdiffusion has the widest range — top de novo designs (RFD_199: 0.726, RFD_167: 0.715) outperform all GLMN sequences, but the tail is weaker.

### Key rank changes (composite3 → composite4)

| Sequence | c3 rank | c4 rank | Shift | ipSAE_tm8 |
|---|---|---|---|---|
| RFD_160_best | 15 | **1** | +14 | 0.631 |
| RFD_56_best | 36 | **3** | +33 | 0.594 |
| RFD_37_best | 18 | **7** | +11 | 0.623 |
| RFD_167_best | 1 | 10 | −9 | 0.715 |
| RFD_199_best | 3 | 14 | −11 | 0.726 |

### 4 sequences with weak ipSAE in submitted 100

| Sequence | ipSAE_tm8 | ipTM | Flag |
|---|---|---|---|
| RFD_29_best | 0.463 | 0.701 | Low confidence interface |
| CUL1_WHB_T0.2_s3 | 0.490 | 0.701 | Induced-fit risk |
| RFD_45_best | 0.493 | 0.766 | Low interface confidence |
| RFD_28_best | 0.494 | 0.723 | Low interface confidence |

### Why our ipSAE is mild — root causes

1. **Short binders (65–95 AA)** — small contact footprint; fewer cross-chain pairs to average over
2. **De novo backbones** — AI-hallucinated geometries have less precise interface packing than evolved proteins
3. **Single ProteinMPNN sequence per backbone** — only 1–2 sequences explored per backbone
4. **No Rosetta interface minimization** — relaxing the interface would tighten packing and improve ipSAE
5. **No pDockQ filter** — Proteinbase shows pDockQ ~0.01–0.02 for our sequences; pDockQ < 0.10 means interface geometry is uncertain

### Design improvement strategies for next round

**Structural:**
- Longer binders (100–150 AA) — more interface area → more confident cross-chain contacts
- Partial diffusion from GLMN backbone — keep natural binding geometry, generate sequence diversity
- Multi-state design — condition on RBX1 + GLMN to complement the natural binding mode

**Filtering (apply before next submission):**
- pDockQ > 0.10 (hard filter — eliminates wrong interface geometry)
- ipSAE_tm8 > 0.58 (removes bottom quartile)
- min_ipSAE_tm8 > 0.52
- Shape complementarity > 50
- Seed consistency: reject if ipSAE varies >0.15 across 5 seeds

**Sequence design:**
- ProteinMPNN: 5–10 sequences per backbone (currently 1–2)
- ProteinMPNN temperature 0.05–0.10 (lower = more conservative = higher ipSAE)
- Interface-biased ProteinMPNN at hotspot positions

**BoltGen (once hotspot bug fixed):**
- Run 50+ designs conditioned on RING-H2 hotspots
- refolding_rmsd_threshold: 2.0
- Compare BoltGen vs RFdiffusion ipSAE distributions — key benchmark

### Output files
- `recompute_ipsae_norm.py` — recomputation script
- `rescored_v2.csv` — all sequences with ipsae_tm8, ipsae_frac5, composite4

---

## Resources

- Competition page: GEM x Adaptyv RBX1 Binder Design Challenge
- Target PDB: [2LGV](https://www.rcsb.org/structure/2LGV)
- RFdiffusion Colab: [RFdiffusion notebook](https://colab.research.google.com/github/sokrypton/ColabDesign/blob/v1.1.1/rf/examples/diffusion.ipynb)
- ProteinMPNN Colab: [ProteinMPNN notebook](https://colab.research.google.com/github/dauparas/ProteinMPNN/blob/main/colab_notebooks/quickdemo.ipynb)
- ColabFold: [ColabFold notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)
