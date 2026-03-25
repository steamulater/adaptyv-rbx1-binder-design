#!/usr/bin/env python3
"""
Compute ipSAE (interface Predicted Structural Assessment Error) for all
Batch 1 and Batch 2 complex predictions.

ipSAE = mean PAE of cross-chain residue pairs, optionally restricted to
        contact-zone residues (within cutoff distance).

We compute two versions (matching Adaptyv's Nipah metrics):
  avg_ipsae     — mean of the full cross-chain PAE block
  min_ipsae     — 5th-percentile of cross-chain PAE (Adaptyv calls this "min")

Both are in Angstroms; lower = better (Boltz-2 is more certain about the interface).
"""
import os, glob, csv
import numpy as np

B1_COMP = "boltz_DESIGN_results/boltz_results_complexes/predictions"
B2_COMP = "boltz_rfdiffusion/complex_results/boltz_results_complex/predictions"

# ── YAML parser to get binder sequence length (= binder chain residue count) ──
def binder_len_from_yaml(seq_id, yaml_dir):
    yaml_path = f"{yaml_dir}/{seq_id}_monomer.yaml"
    if not os.path.exists(yaml_path):
        return None
    with open(yaml_path) as f:
        for line in f:
            if "sequence:" in line:
                return len(line.split("sequence:")[-1].strip())
    return None

# For Batch 1 use master_sequences.csv for binder length
b1_lengths = {}
with open("master_sequences.csv") as f:
    for row in csv.DictReader(f):
        b1_lengths[row["seq_id"]] = int(row["length"])

RBX1_LEN = 108   # full RBX1 sequence used in complex YAML

def compute_ipsae(pred_dir, binder_len):
    """
    Load all PAE matrices in pred_dir, extract cross-chain block,
    return (avg_ipsae, min_ipsae) averaged across models.

    Convention: chain A = binder (first binder_len residues)
                chain B = RBX1  (remaining residues)
    """
    pae_files = sorted(glob.glob(f"{pred_dir}/pae_*.npz"))
    if not pae_files:
        return None, None

    avg_vals, min_vals = [], []
    for fp in pae_files:
        pae = np.load(fp)["pae"]          # shape (N, N), Angstroms
        n   = pae.shape[0]
        bl  = binder_len
        rb  = n - bl                       # RBX1 residue count in this prediction

        if bl <= 0 or rb <= 0 or bl >= n:
            continue

        # Cross-chain blocks:
        # binder→RBX1:  pae[0:bl,  bl:n]
        # RBX1→binder:  pae[bl:n,  0:bl]
        cross = np.concatenate([
            pae[:bl, bl:].flatten(),
            pae[bl:, :bl].flatten()
        ])

        avg_vals.append(cross.mean())
        # "min_ipsae" in Adaptyv = 5th percentile (captures the tightest interface contacts)
        min_vals.append(np.percentile(cross, 5))

    if not avg_vals:
        return None, None
    return float(np.mean(avg_vals)), float(np.mean(min_vals))


# ── Process Batch 1 ───────────────────────────────────────────────────────────
print("Processing Batch 1...")
b1_results = {}
for pred_dir in sorted(glob.glob(f"{B1_COMP}/*_complex")):
    seq_id = os.path.basename(pred_dir).replace("_complex", "")
    blen   = b1_lengths.get(seq_id)
    if blen is None:
        continue
    avg_ip, min_ip = compute_ipsae(pred_dir, blen)
    b1_results[seq_id] = {"avg_ipsae": avg_ip, "min_ipsae": min_ip}

print(f"  Done: {len(b1_results)} sequences")

# ── Process Batch 2 ───────────────────────────────────────────────────────────
print("Processing Batch 2...")
B2_YAML = "boltz_rfdiffusion/monomer"
b2_results = {}
for pred_dir in sorted(glob.glob(f"{B2_COMP}/*_complex")):
    seq_id = os.path.basename(pred_dir).replace("_complex", "")
    blen   = binder_len_from_yaml(seq_id, B2_YAML)
    if blen is None:
        continue
    avg_ip, min_ip = compute_ipsae(pred_dir, blen)
    b2_results[seq_id] = {"avg_ipsae": avg_ip, "min_ipsae": min_ip}

print(f"  Done: {len(b2_results)} sequences")

# ── Summary stats ─────────────────────────────────────────────────────────────
def stats(vals, label):
    vals = [v for v in vals if v is not None]
    if not vals: return
    print(f"  {label}: mean={np.mean(vals):.2f}  p25={np.percentile(vals,25):.2f}  "
          f"median={np.median(vals):.2f}  p75={np.percentile(vals,75):.2f}  min={min(vals):.2f}  max={max(vals):.2f}")

print("\nBatch 1 avg_ipSAE (Å):")
stats([v["avg_ipsae"] for v in b1_results.values()], "all")

glmn_ids = [k for k in b1_results if k.startswith("GLMN")]
whb_ids  = [k for k in b1_results if k.startswith("CUL1")]
stats([b1_results[k]["avg_ipsae"] for k in glmn_ids], "GLMN")
stats([b1_results[k]["avg_ipsae"] for k in whb_ids],  "CUL1_WHB")

print("\nBatch 2 avg_ipSAE (Å):")
stats([v["avg_ipsae"] for v in b2_results.values()], "all")

# ── Nipah reference thresholds (from retrospective) ───────────────────────────
# Nipah: binders avg_ipsae ~0.73 (0-1 normalised) → but Nipah stored 0-1 normalised
# Our PAE is in raw Angstroms (0–30 Å range)
# We need to compare within our own distribution; lower = better
print("\n(ipSAE in raw Å — lower is better. Nipah stored normalised 0-1; not directly comparable.)")
print("Within our data, use as ranking: prefer lower avg_ipsae among passing sequences.")

# ── Write results ─────────────────────────────────────────────────────────────
out_rows = []
for seq_id, vals in {**b1_results, **b2_results}.items():
    out_rows.append({"seq_id": seq_id, **vals})

out_rows.sort(key=lambda r: r["avg_ipsae"] or 999)

with open("ipsae_results.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["seq_id", "avg_ipsae", "min_ipsae"])
    w.writeheader()
    w.writerows(out_rows)

print(f"\nSaved ipsae_results.csv ({len(out_rows)} sequences)")
print("\nTop 15 by avg_ipSAE (lowest = most confident interface):")
print(f"{'seq_id':<35} {'avg_ipSAE':>10} {'min_ipSAE':>10}")
print("-" * 58)
for r in out_rows[:15]:
    print(f"{r['seq_id']:<35} {r['avg_ipsae']:>10.2f} {r['min_ipsae']:>10.2f}")
