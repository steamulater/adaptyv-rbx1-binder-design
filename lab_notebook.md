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

Generate entirely new binder backbones conditioned on the RBX1 RING-H2 surface, then design sequences with ProteinMPNN.

```
RBX1 structure (2LGV)
        |
        v
[1] RFdiffusion  [Colab A100]
    - Hotspot residues: RING-H2 surface (res 45-90)
    - Generate ~500 binder backbones
        |
        v
[2] ProteinMPNN  [local, M3 MPS]
    - Design sequences for each backbone
    - 8 sequences per backbone
        |
        v
[3] ColabFold AF2-multimer  [Colab A100]
    - Predict RBX1:binder complex structures
    - Score by ipTM, pLDDT, interface contacts
        |
        v
[4] Novelty filter
    - mmseqs2 vs UniRef50
    - Retain sequences with ≥25% edit distance
        |
        v
[5] Rank & select top N
    - ipTM > 0.7 cutoff
    - Interface score (ΔG, buried SASA)
    - Sequence diversity
```

**Rationale:** Maximally novel binders; not constrained to known binding modes.

---

### Strategy 2a: Scaffold-Based Redesign — Glomulin (4F52)

Use the experimentally validated Glomulin backbone — a natural inhibitor of CRL activity that occupies the E2-docking surface of RBX1's RING domain — as the starting scaffold for ProteinMPNN sequence redesign.

**Scaffold:** PDB 4F52, chain F (Glomulin), residues 336–582 | 247 AA | UniProt Q92990
**Missing residues in crystal structure:** 433–438 (6 AA), 533–549 (17 AA) — to be modelled with Boltz-1

```
PDB 4F52 (Glomulin-RBX1-CUL1 crystal structure)
        |
        v
[1] PyMOL  [local]
    - Align 4F52 against 2LGV (RMSD = 7.692 Å on 88 atoms)
    - Extract Glomulin chain F, residues 336-582 (247 AA)
    - Save as Scaffold_4F52_336-582.pdb
        |
        v
[2] Boltz-1  [ColabFold + tier]
    - Monomer: glmn_monomer.yaml — predict complete structure with missing loops
    - Complex: glmn_rbx1_complex.yaml — validate PPI with RBX1
    - 5 diffusion samples each; assess ipTM consistency across models
        |
        v
[3] ProteinMPNN  [local, M3 MPS]
    - Fixed receptor: RBX1 RING-H2
    - Redesign Glomulin scaffold sequence
    - 8-16 sequences per run, multiple temperature settings
        |
        v
[4] Boltz-1 / ColabFold AF2-multimer  [Colab A100]
    - Validate redesigned sequences in complex with RBX1
    - Score by ipTM, pLDDT
        |
        v
[5] Novelty filter → mmseqs2 vs UniRef50 (≥25% edit distance)
        |
        v
[6] Rank & select top N
```

**Rationale:** Glomulin is a natural CRL inhibitor; backbone geometry is experimentally validated against the exact E2-docking surface. 247 AA — at the competition limit.

**Key structural references:**
- 4F52 scaffold: chain F, residues 336–582
- Aligned to 2LGV (RBX1 NMR): RMSD = 7.692 Å on 88 atoms

![4F52 Glomulin scaffold in PyMOL — orange: Glomulin chain F res 336-582; green: CRL complex](Scaffold_4F52_582-336.png)

---

### Strategy 2b: Scaffold-Based Redesign — CUL1 WHB domain (1LDJ)

Use the C-terminal winged-helix B (WHB) domain of Cullin-1 that directly cradles the RBX1 RING domain as a compact (72 AA) scaffold for ProteinMPNN redesign.

**Scaffold:** PDB 1LDJ, chain A (CUL1), residues 705–776 | 72 AA | UniProt Q13616
**Sequence:** `EDRKLLIQAAIVRIMKMRKVLKHQQLLGEVLTQLSSRFKPRVPVIKKCIDILIEKEYLERVDGEKDTYSYLA`

```
PDB 1LDJ (CUL1-RBX1-SKP1-SKP2 SCF complex)
        |
        v
[1] PyMOL  [local]
    - Extract CUL1 chain A, residues 705-776 (72 AA)
    - Save as Scaffold_1LDJ_705-776.pdb
        |
        v
[2] Boltz-1  [ColabFold + tier]
    - Monomer: cul1_whb_monomer.yaml — predict standalone WHB structure
    - Complex: cul1_whb_rbx1_complex.yaml — validate PPI with RBX1
    - 5 diffusion samples each; assess ipTM consistency across models
        |
        v
[3] ProteinMPNN  [local, M3 MPS]
    - Redesign WHB scaffold sequence
    - 8-16 sequences per run, multiple temperature settings
        |
        v
[4] Boltz-1 / ColabFold AF2-multimer validation → novelty filter → rank
```

**Rationale:** Much more compact than Glomulin (72 vs 247 AA), leaving headroom for linkers or extensions. The WHB domain is the direct structural contact between Cullin and the RBX1 RING domain.

**Key structural references:**
- 1LDJ scaffold: chain A, residues 705–776
- Aligned to 2LGV (RBX1 NMR): RMSD = 7.692 Å on 88 atoms

![1LDJ CUL1 WHB domain scaffold in PyMOL — yellow: CUL1 chain A res 705-776; red: RBX1](Scaffold_1LDJ_705-776.png)

**Compute:** ProteinMPNN runs locally on M3 (MPS-accelerated); Boltz-1 and validation on ColabFold + tier ($50/month)

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

## Sequences Submitted

| # | Sequence ID | Length (AA) | ipTM | pLDDT | Edit dist. | Notes |
|---|-------------|-------------|------|-------|------------|-------|
| — | — | — | — | — | — | — |

---

## Resources

- Competition page: GEM x Adaptyv RBX1 Binder Design Challenge
- Target PDB: [2LGV](https://www.rcsb.org/structure/2LGV)
- RFdiffusion Colab: [RFdiffusion notebook](https://colab.research.google.com/github/sokrypton/ColabDesign/blob/v1.1.1/rf/examples/diffusion.ipynb)
- ProteinMPNN Colab: [ProteinMPNN notebook](https://colab.research.google.com/github/dauparas/ProteinMPNN/blob/main/colab_notebooks/quickdemo.ipynb)
- ColabFold: [ColabFold notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)
