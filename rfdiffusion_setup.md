# RFdiffusion — De Novo Binder Design for RBX1
## Strategy 1: De Novo Backbone Generation

---

## Overview

Generate 100–200 novel binder backbones conditioned on the RBX1 RING domain interface.
Then sequence-design with ProteinMPNN → validate with Boltz-2 → select top 45.

**Target:** RBX1 RING-H2 domain
**Input PDB:** `rbx1_ring_renumbered.pdb` (chain A, residues 1–77)
**Original numbering:** RBX1 residues 32–108 (structured RING domain; 1–31 are disordered)

---

## Step 1 — RFdiffusion (Colab A100)

### Files needed (upload to Colab)
- `rbx1_ring_renumbered.pdb` — RBX1 RING domain, chain A, residues 1–77

### RFdiffusion parameters

**Notebook:** [RFdiffusion Colab](https://colab.research.google.com/github/sokrypton/ColabDesign/blob/v1.1.1/rf/examples/diffusion.ipynb)

```python
# In RFdiffusion Colab, set these parameters:

name = "rbx1_binder"

# Target chain: chain A (RBX1 RING domain)
# Binder: generate new chain, length 60-100 AA

contigs = "A1-77/0 60-100"
# This keeps the full RBX1 RING domain fixed and generates
# a new chain of 60-100 AA that contacts the hotspot surface

hotspot_res = "A4,A6,A12,A14,A15,A23,A24,A26,A28,A52,A56,A60,A64,A66"
# These are the core RING-H2 surface residues that GLMN contacts
# (original RBX1 numbering: 35,37,43,45,46,54,55,57,59,83,87,91,95,97)

num_designs = 200  # Generate 200 backbones; expect ~50% to pass downstream filters

# Noise schedule (use defaults, or try):
noise_scale_ca = 1.0
noise_scale_frame = 1.0
```

**What RFdiffusion will do:**
- Keep RBX1 RING domain fixed in place
- Diffuse a new protein backbone (60–100 AA) that makes contacts with the hotspot residues
- Generate 200 candidate binder backbones as separate PDB files

### Expected output
200 PDB files named `rbx1_binder_0.pdb`, `rbx1_binder_1.pdb`, etc.
Each contains 2 chains: chain A (RBX1, fixed) + chain B (designed binder backbone).

---

## Step 2 — Filter Backbones by pLDDT

RFdiffusion outputs a per-residue pLDDT for the designed binder chain.
Before running ProteinMPNN, filter to backbones with **mean binder pLDDT > 0.80**
(ensures the backbone is well-defined, not disordered).

```python
# In Colab after RFdiffusion:
# Filter designs by binder chain pLDDT
# Keep designs where binder chain mean pLDDT > 0.80
# Expect ~100-120 designs to pass (50-60% pass rate typical)
```

---

## Step 3 — ProteinMPNN Sequence Design (Local, M3 MPS)

For each passing backbone, run ProteinMPNN to design sequences.

### Settings
```bash
# Run ProteinMPNN on the filtered backbones
# Design only chain B (binder) — keep chain A (RBX1) fixed

python protein_mpnn_run.py \
  --pdb_path rbx1_binder_filtered/ \
  --chain_id_jsonl chain_ids.jsonl \       # {"name": ["B"]} for each design
  --fixed_positions_jsonl fixed_pos.jsonl \ # fix nothing in chain B
  --out_folder mpnn_outputs/ \
  --num_seq_per_target 5 \
  --sampling_temp "0.1 0.2 0.3" \
  --batch_size 1
```

**Target:** ~100 filtered backbones × 3 sequences each = 300 sequences
→ After Boltz-2 filtering (ipTM ≥ 0.70), expect ~45–60 to pass
→ Select top 45 for submission

---

## Step 4 — Boltz-2 Validation (Colab)

Same pipeline as Batch 1:
- Monomer prediction: does the binder fold correctly?
- Complex prediction with RBX1 (full sequence, UniProt P62877): does it bind?
- Filter: avg ipTM ≥ 0.70, monomer pTM ≥ 0.70

---

## Step 5 — Novelty Screen

DIAMOND blastp against SwissProt (same as Batch 1).
Expect RFdiffusion-designed backbones to be highly novel (<30% identity to anything).

---

## Design Rationale

**Why RFdiffusion now (after GLMN)?**

The GLMN scaffold batch (48 sequences) represents a single binding mode — all share the GLMN backbone and contact the same RBX1 surface patch from the same geometry. If this geometry doesn't work experimentally for any scaffold-level reason, the entire batch fails together.

RFdiffusion generates completely new backbones with no structural similarity to GLMN:
- Different secondary structure composition (may not be all-helical like GLMN)
- Different binding angle and contact geometry
- Shorter length (60–100 AA vs 247 AA)
- True de novo sequences: expected <20% identity to any natural protein

These 45 sequences provide independent structural coverage of the RBX1 RING surface,
significantly de-risking the overall submission.

---

## Hotspot Context (for reference)

The 14 core hotspot residues (original RBX1 numbering):

| RBX1 res | Residue | Role |
|----------|---------|------|
| 35 | TRP | Hydrophobic core of RING surface |
| 37 | ILE | Hydrophobic |
| 43 | ALA | Interface packing |
| 45 | CYS | Zn-coordinating (structural) |
| 46 | ARG | Charged contact |
| 54 | ILE | Hydrophobic |
| 55 | GLU | H-bond donor/acceptor |
| 57 | GLN | H-bond |
| 59 | ASN | Polar contact |
| 83 | CYS | Zn-coordinating (structural) |
| 87 | TRP | Large hydrophobic, critical anchor |
| 91 | ARG | Salt bridge candidate |
| 95 | PRO | Structural turn |
| 97 | ASP | Charged contact |

In `rbx1_ring_renumbered.pdb` (chain A, 1–77):
subtract 31 from each residue number.
Hotspot string: `A4,A6,A12,A14,A15,A23,A24,A26,A28,A52,A56,A60,A64,A66`

---

## Timeline

| Step | Tool | Location | Time estimate |
|------|------|----------|--------------|
| RFdiffusion 200 designs | RFdiffusion | Colab A100 | ~30–45 min |
| Filter by pLDDT | Python | Local | <5 min |
| ProteinMPNN | ProteinMPNN | Local (M3 MPS) | ~20 min |
| Boltz-2 monomer (300 seqs) | Boltz-2 | Colab | ~2 hr |
| Boltz-2 complex (300 seqs) | Boltz-2 | Colab | ~2 hr |
| Novelty + selection | DIAMOND | Local | <10 min |
| **Total** | | | **~5–6 hr** |

**Deadline:** March 26, 2026 (3 days from now — feasible)
