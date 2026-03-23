# Lab Notebook: RBX1 Binder Design
**Competition:** GEM x Adaptyv Bio — RBX1 Binder Design Challenge
**Submission deadline:** March 26, 2026
**Experimenter:** Tamuka Martin

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

## Computational Pipeline

```
RBX1 structure (2LGV)
        |
        v
[1] RFdiffusion
    - Hotspot residues: RING-H2 surface (res 45-90)
    - Generate ~500 binder backbones
        |
        v
[2] ProteinMPNN
    - Design sequences for each backbone
    - 8 sequences per backbone
        |
        v
[3] ColabFold (AF2-multimer)
    - Predict RBX1:binder complex structures
    - Score by ipTM, pLDDT, interface contacts
        |
        v
[4] Novelty filter
    - mmseqs2 vs UniRef50
    - Retain sequences with ≥25% edit distance
        |
        v
[5] Rank & select top 100
    - ipTM > 0.7 cutoff
    - Interface score (ΔG, buried SASA)
    - Sequence diversity
```

**Compute:** ColabFold + tier ($50/month) — Google Colab A100 GPUs

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

### Entry 002 — [date]

**Status:** [in progress / complete]

**Work completed:**
-

**Results:**
-

**Next steps:**
-

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
