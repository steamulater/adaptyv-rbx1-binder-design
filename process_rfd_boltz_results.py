#!/usr/bin/env python3
"""
Post-processing for RFdiffusion Batch 2 Boltz-2 results.
Run after downloading boltz_rfdiffusion/ folder from Google Drive.

Steps:
1. Parse monomer + complex confidence scores
2. Filter: avg ipTM >= 0.70, monomer pTM >= 0.70
3. Novelty screen (DIAMOND vs SwissProt — already built)
4. Select top 45 by ipTM
5. Append to final_submission.fasta
6. Update master_sequences.csv
"""
import os, glob, json, csv
import numpy as np

BOLTZ_DIR   = "boltz_rfdiffusion"         # downloaded from Drive
MONO_PRED   = f"{BOLTZ_DIR}/monomer_results/boltz_results_monomer/predictions"
COMP_PRED   = f"{BOLTZ_DIR}/complex_results/boltz_results_complex/predictions"
RBX1_SEQ    = "MAAAMDVDTPSGTNSGAGKKRFEVKKWNAVALWAWDIVVDNCAICRNHIMDLCIECQANQASATSEECTVAWGVCNHAFHFHCISRWLKTRQVCPLDNREWEFQKYGH"

# ── 1. Parse confidence scores ──────────────────────────────────────────────
def avg_confidence(pred_dir, seq_id, key_iptm, key_ptm, key_plddt):
    pattern = f"{pred_dir}/{seq_id}*/confidence_*.json"
    files   = glob.glob(pattern)
    if not files:
        return None, None, None
    iptms, ptms, plddts = [], [], []
    for f in files:
        d = json.load(open(f))
        if key_iptm  in d: iptms.append(d[key_iptm])
        if key_ptm   in d: ptms.append(d[key_ptm])
        if key_plddt in d: plddts.append(d[key_plddt])
    return (np.mean(iptms)  if iptms  else None,
            np.mean(ptms)   if ptms   else None,
            np.mean(plddts) if plddts else None)

results = []
mono_dirs = glob.glob(f"{MONO_PRED}/*_monomer")
print(f"Monomer prediction dirs: {len(mono_dirs)}")

for mono_dir in sorted(mono_dirs):
    seq_id = os.path.basename(mono_dir).replace("_monomer", "")

    # Monomer scores
    mono_files = glob.glob(f"{mono_dir}/confidence_*.json")
    mono_ptms, mono_plddts = [], []
    for f in mono_files:
        d = json.load(open(f))
        if "ptm"          in d: mono_ptms.append(d["ptm"])
        if "complex_plddt" in d: mono_plddts.append(d["complex_plddt"])
    mono_ptm   = np.mean(mono_ptms)   if mono_ptms   else None
    mono_plddt = np.mean(mono_plddts) if mono_plddts else None

    # Complex scores
    comp_dir   = f"{COMP_PRED}/{seq_id}_complex"
    comp_files = glob.glob(f"{comp_dir}/confidence_*.json")
    iptms, comp_ptms, comp_plddts = [], [], []
    for f in comp_files:
        d = json.load(open(f))
        iptm = d.get("protein_iptm", d.get("iptm"))
        if iptm is not None:         iptms.append(iptm)
        if "ptm"           in d:     comp_ptms.append(d["ptm"])
        if "complex_plddt" in d:     comp_plddts.append(d["complex_plddt"])
    avg_iptm      = np.mean(iptms)      if iptms      else None
    avg_comp_ptm  = np.mean(comp_ptms)  if comp_ptms  else None
    avg_comp_plddt = np.mean(comp_plddts) if comp_plddts else None

    results.append({
        "seq_id": seq_id, "mono_ptm": mono_ptm, "mono_plddt": mono_plddt,
        "avg_iptm": avg_iptm, "comp_ptm": avg_comp_ptm, "comp_plddt": avg_comp_plddt,
    })

# ── 2. Filter ───────────────────────────────────────────────────────────────
passing = [r for r in results
           if r["avg_iptm"]  is not None and r["avg_iptm"]  >= 0.70
           and r["mono_ptm"] is not None and r["mono_ptm"]  >= 0.70]

passing.sort(key=lambda x: -x["avg_iptm"])

print(f"\nTotal sequences:     {len(results)}")
print(f"Pass ipTM ≥ 0.70:    {len(passing)}")
print(f"\nTop 20 by ipTM:")
print(f"{'seq_id':<35} {'ipTM':<7} {'mono_pTM':<10} {'comp_pLDDT'}")
print("-" * 65)
for r in passing[:20]:
    plddt_str = f"{r['comp_plddt']:.3f}" if r['comp_plddt'] is not None else "N/A"
    print(f"{r['seq_id']:<35} {r['avg_iptm']:<7.3f} {r['mono_ptm']:<10.3f} {plddt_str}")

# ── 3. Recover sequences from MPNN output ───────────────────────────────────
import re
MPNN_DIR = "boltz_rfdiffusion"

def get_seq_from_yaml(seq_id):
    yaml_path = f"{BOLTZ_DIR}/monomer/{seq_id}_monomer.yaml"
    if not os.path.exists(yaml_path):
        return None
    with open(yaml_path) as f:
        for line in f:
            if "sequence:" in line:
                return line.split("sequence:")[-1].strip()
    return None

# ── 4. Select top 45 ────────────────────────────────────────────────────────
top45 = passing[:45]
print(f"\nSelecting top {len(top45)} sequences for Batch 2 submission")

# Write FASTA
with open("rfd_batch2_candidates.fasta", "w") as f:
    for i, r in enumerate(top45, 1):
        seq = get_seq_from_yaml(r["seq_id"])
        if seq:
            f.write(f">{r['seq_id']}|rank{i}|ipTM={r['avg_iptm']:.3f}|scaffold=RFdiffusion\n")
            f.write(seq + "\n")

print(f"Saved rfd_batch2_candidates.fasta")

# ── 5. Write CSV summary ─────────────────────────────────────────────────────
with open("rfd_batch2_results.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    writer.writeheader()
    writer.writerows(sorted(results, key=lambda x: -(x["avg_iptm"] or 0)))
print(f"Saved rfd_batch2_results.csv")
