# Lab Notebook: ProteWorks Studio — BoltGen Design Sprint
**Platform:** ProteWorks Studio (studio.proteworks.com) — BoltGen pipeline
**Competition:** GEM x Adaptyv Bio — RBX1 Binder Design Challenge
**Objective:** Generate novel RBX1 binder candidates using BoltGen; concurrently stress-test the ProteWorks platform and log structured feedback for Brandon (founder)
**Experimenter:** Tamuka Martin Chidyausiku @steamulater
**Notebook started:** March 25, 2026

---

## Context

This is a parallel design sprint running alongside (and after) the primary RFdiffusion + ProteinMPNN pipeline documented in `../lab_notebook.md`. That pipeline produced 100 sequences submitted to the competition on March 25, 2026.

This notebook documents a second wave of designs generated entirely within ProteWorks Studio, using their BoltGen workflow. Two goals run in parallel:

1. **Science:** Produce high-quality RBX1 binder candidates for future submission rounds or follow-up experimental work. Apply everything learned from the Nipah retrospective (composite scoring, ipSAE weighting, RING RMSD screening) to evaluate outputs rigorously.

2. **Platform evaluation:** Stress-test ProteWorks Studio end-to-end — input handling, job management, output quality, UI/UX friction — and compile actionable feedback for Brandon. Treat this like a structured beta test with a real design problem.

---

## Target: RBX1 (same as primary campaign)

**UniProt:** P62877 | **PDB:** 2LGV | **Organism:** *Homo sapiens*

**Full sequence (108 AA):**
```
MAAAMDVDTPSGTNSGAGKKRFEVKKWNAVALWAWDIVVDNCAICRNHIMDLCIECQANQASATSEECTVAWGVCNHAFHFHCISRWLKTRQVCPLDNREWEFQKYGH
```

**RING domain (target region, residues 32–90):**
```
NCAICRNHIMDLCIECQANQASATSEECTVAWGVCNHAFHFHCISRWLKTRQVCPLDNREWEFQKYGH
```

**Hotspot residues (renumbered RING domain):** `4, 6, 12, 14, 15, 23, 24, 26, 28, 52, 56, 60, 64, 66`

**Therapeutic rationale:** RBX1 is the catalytic RING subunit of Cullin-RING E3 ligase complexes, directing ubiquitination of ~20% of all cellular proteins. Blocking the E2-binding surface disrupts oncogenic CRL activity. Natural inhibitor: Glomulin (GLMN), which occupies the RING-H2 face with high affinity.

**Key structural files (in `../nanome_vr/`):**
- `00_RBX1_target_RING_domain.pdb` — isolated RING domain
- `01_RFD167_top_denovo_rank1.cif` — best de novo binder (ipTM 0.905)
- `07_GLMN_s4_scaffold_top.pdb` — best GLMN scaffold design (ipTM 0.918)

---

## Scoring Framework (from Nipah retrospective)

Metrics ranked by AUROC on 1,030 experimentally validated binders:

| Metric | AUROC | Use |
|---|---|---|
| complex_ipLDDT | 0.691 | Primary confidence metric |
| min_ipSAE | 0.638 | Interface tightness |
| ipTM | 0.603 | Structural compatibility |
| mono pTM | 0.501 | Not used (random) |

**Composite score:** `composite = 0.4 × ipTM + 0.3 × ipLDDT + 0.3 × norm_ipSAE`

**Hard filters:**
- ipTM >= 0.70
- No poly-Ala or low-complexity sequences
- RING RMSD < 3.0 Å (RBX1 must not deform to accommodate binder)

**Soft thresholds (relative, not absolute):**
- ipLDDT: rank within cohort, do not apply Nipah's 0.85 cutoff
- ipSAE: prefer lower Å; Nipah's 0-1 scale is not directly comparable

---

## Platform: ProteWorks Studio

**URL:** https://studio.proteworks.com/
**Tool:** BoltGen
**Contact:** Brandon (founder)

### What BoltGen is (working understanding)
BoltGen is a web-based protein design tool built on top of Boltz-2. Unlike the Colab-based pipeline used in the primary campaign (RFdiffusion → ProteinMPNN → Boltz-2), BoltGen appears to integrate structure prediction and sequence design into a single interface. To be confirmed through use.

### Evaluation criteria for platform feedback

Track the following for each session:

| Category | Things to note |
|---|---|
| **Input UX** | How is target specified? PDB upload, UniProt ID, sequence? Hotspot specification workflow |
| **Job submission** | Queue times, job naming, batch submission capability |
| **Output quality** | Metrics returned (ipTM? pLDDT? PAE? ipSAE?), structure download format |
| **Output UX** | Result dashboard, sorting/filtering, export options |
| **Reproducibility** | Same inputs → same outputs? Seed control? |
| **Failure modes** | What breaks, what errors look like, recovery path |
| **Speed** | Wall time per job, throughput |
| **Gaps vs Colab pipeline** | What can't be done that we needed before |

---

## Design Plan

### Round 1 — Baseline replication
**Goal:** Reproduce a known-good result from the primary campaign to calibrate BoltGen output quality.
- Input: RBX1 RING domain + hotspot residues (same as primary campaign)
- Design a small set (5–10 sequences) with default BoltGen parameters
- Score with same composite metric
- Compare ipTM/ipLDDT distribution vs primary campaign baseline

### Round 2 — Parameter exploration
**Goal:** Explore BoltGen-specific settings (if any) to understand design space coverage.
- Vary whatever parameters BoltGen exposes (temperature, number of sequences, etc.)
- Document what is and isn't controllable vs Colab pipeline

### Round 3 — Novel design angles
**Goal:** Use BoltGen for design directions not explored in primary campaign.
- Shorter binders (< 60 AA) — if BoltGen supports constrained length
- Beta-sheet dominated scaffolds — primary campaign was helix-heavy
- Mixed alpha/beta if platform generates diverse topologies

---

## Entries

---

### Entry 001 — Platform First Access
**Date:** 2026-03-25
**Status:** IN PROGRESS

#### Session goals
- [ ] Create account / log in
- [ ] Explore dashboard and available tools
- [ ] Understand BoltGen input requirements
- [ ] Submit first test job with RBX1 RING domain

#### Platform notes
*(Fill in as you go)*

**Input interface:**


**Job submission:**


**Parameters available:**


**Queue / runtime:**


#### Feedback for Brandon
*(Log friction points, bugs, missing features, things that worked well)*

| # | Category | Observation | Severity |
|---|---|---|---|
| 1 | | | |

---

## Platform Feedback Log (running)

This section accumulates all feedback across sessions for eventual handoff to Brandon.

| Session | Category | Observation | Severity (Low/Med/High) | Suggested fix |
|---|---|---|---|---|
| | | | | |

---

## Results

*(Populated as designs come in)*

| ID | Length | ipTM | ipLDDT | ipSAE (Å) | composite | RING RMSD | Notes |
|---|---|---|---|---|---|---|---|
| | | | | | | | |

---

## References

- Adaptyv competition: https://design.adaptyv.bio
- Primary campaign notebook: `../lab_notebook.md`
- Nipah retrospective analysis: `../nipah_analysis/`
- Composite scoring rationale: `../submission_writeup.md`
- Boltz-2: Wohlwend et al. 2024, bioRxiv
- RFdiffusion (primary campaign): Watson et al. 2023, Nature
- ProteinMPNN (primary campaign): Dauparas et al. 2022, Science
