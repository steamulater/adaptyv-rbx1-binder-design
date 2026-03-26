#!/usr/bin/env python3
"""
Recompute normalised ipSAE for all 100 submitted sequences.

Proteinbase ipSAE is 0-1, higher = better. Their pipeline uses full-length
RBX1 (108 AA) so we cannot exactly replicate it from our RING-domain PAE
matrices. Instead we compute two internally-consistent normalised metrics:

  ipsae_tm8  : mean of TM-score-like function  1/(1 + (pae_ij/8)^2)
               over all cross-chain residue pairs. Range 0-1, higher = better.
               d0=8 Å is the canonical value used in pTM scoring (AF2/Boltz).

  ipsae_frac5 : fraction of cross-chain pairs with PAE < 5 Å.
                Stricter than the 10 Å cutoff — captures only high-confidence
                contacts. Range 0-1, higher = better.

Both replace the raw-Å metric we computed earlier. Use ipsae_tm8 as the
primary ranking metric (closest in spirit to Proteinbase's approach).

Updated composite score:
  composite4 = 0.35 * ipTM  +  0.30 * ipLDDT  +  0.35 * ipsae_tm8
(upweighted ipSAE per Nipah retrospective: AUROC rank ipSAE > ipTM)
"""

import os, glob, csv, json
import numpy as np

B1_COMP = "boltz_DESIGN_results/boltz_results_complexes/predictions"
B2_COMP = "boltz_rfdiffusion/complex_results/boltz_results_complex/predictions"

# ── Binder lengths from master_sequences.csv ──────────────────────────────────
b1_lengths = {}
with open("master_sequences.csv") as f:
    for row in csv.DictReader(f):
        b1_lengths[row["seq_id"]] = int(row["length"])

def binder_len_from_yaml(seq_id, yaml_dir):
    path = f"{yaml_dir}/{seq_id}_monomer.yaml"
    if not os.path.exists(path):
        return None
    with open(path) as f:
        for line in f:
            if "sequence:" in line:
                return len(line.split("sequence:")[-1].strip())
    return None

# ── Core computation ───────────────────────────────────────────────────────────
def compute_ipsae_norm(pred_dir, binder_len):
    pae_files = sorted(glob.glob(f"{pred_dir}/pae_*.npz"))
    if not pae_files:
        return None, None

    tm8_vals, frac5_vals = [], []
    min_tm8_vals, min_frac5_vals = [], []

    for fp in pae_files:
        pae  = np.load(fp)["pae"]          # shape (N, N), Angstroms
        n, bl = pae.shape[0], binder_len

        if bl <= 0 or n - bl <= 0 or bl >= n:
            continue

        b2r = pae[:bl, bl:]   # binder→RBX1
        r2b = pae[bl:, :bl]   # RBX1→binder

        # TM-like score (d0 = 8 Å) — each pair contributes 0-1
        tm8_b2r  = np.mean(1 / (1 + (b2r / 8) ** 2))
        tm8_r2b  = np.mean(1 / (1 + (r2b / 8) ** 2))
        tm8_vals.append((tm8_b2r + tm8_r2b) / 2)
        min_tm8_vals.append(min(tm8_b2r, tm8_r2b))

        # Strict fraction: PAE < 5 Å
        frac5_b2r = np.mean(b2r < 5.0)
        frac5_r2b = np.mean(r2b < 5.0)
        frac5_vals.append((frac5_b2r + frac5_r2b) / 2)
        min_frac5_vals.append(min(frac5_b2r, frac5_r2b))

    if not tm8_vals:
        return None, None

    return {
        "ipsae_tm8":      float(np.mean(tm8_vals)),
        "min_ipsae_tm8":  float(np.mean(min_tm8_vals)),
        "ipsae_frac5":    float(np.mean(frac5_vals)),
        "min_ipsae_frac5":float(np.mean(min_frac5_vals)),
    }

# ── Process Batch 1 ───────────────────────────────────────────────────────────
print("Processing Batch 1...")
b1_new = {}
for pred_dir in sorted(glob.glob(f"{B1_COMP}/*_complex")):
    seq_id = os.path.basename(pred_dir).replace("_complex", "")
    blen   = b1_lengths.get(seq_id)
    if blen is None:
        continue
    result = compute_ipsae_norm(pred_dir, blen)
    if result:
        b1_new[seq_id] = result
print(f"  Done: {len(b1_new)} sequences")

# ── Process Batch 2 ───────────────────────────────────────────────────────────
print("Processing Batch 2...")
B2_YAML = "boltz_rfdiffusion/monomer"
b2_new = {}
for pred_dir in sorted(glob.glob(f"{B2_COMP}/*_complex")):
    seq_id = os.path.basename(pred_dir).replace("_complex", "")
    blen   = binder_len_from_yaml(seq_id, B2_YAML)
    if blen is None:
        continue
    result = compute_ipsae_norm(pred_dir, blen)
    if result:
        b2_new[seq_id] = result
print(f"  Done: {len(b2_new)} sequences")

all_new = {**b1_new, **b2_new}

# ── Load existing rescored_all.csv and merge ──────────────────────────────────
print("\nMerging with rescored_all.csv...")
rows = []
with open("rescored_all.csv") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames + ["ipsae_tm8", "min_ipsae_tm8",
                                       "ipsae_frac5", "min_ipsae_frac5",
                                       "composite4"]
    for row in reader:
        sid = row["seq_id"]
        new = all_new.get(sid, {})
        row["ipsae_tm8"]       = f"{new.get('ipsae_tm8', ''):.4f}"      if new.get("ipsae_tm8")      else ""
        row["min_ipsae_tm8"]   = f"{new.get('min_ipsae_tm8', ''):.4f}"  if new.get("min_ipsae_tm8")  else ""
        row["ipsae_frac5"]     = f"{new.get('ipsae_frac5', ''):.4f}"    if new.get("ipsae_frac5")    else ""
        row["min_ipsae_frac5"] = f"{new.get('min_ipsae_frac5',''):.4f}" if new.get("min_ipsae_frac5")else ""

        # composite4: upweight ipSAE
        try:
            iptm   = float(row["avg_iptm"])
            iplddt = float(row["avg_iplddt"])
            isae   = float(row["ipsae_tm8"])
            row["composite4"] = f"{0.35*iptm + 0.30*iplddt + 0.35*isae:.4f}"
        except (ValueError, KeyError):
            row["composite4"] = ""

        rows.append(row)

# Sort by composite4 descending
rows.sort(key=lambda r: float(r["composite4"]) if r["composite4"] else 0, reverse=True)

with open("rescored_v2.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

print("Saved rescored_v2.csv")

# ── Summary stats ─────────────────────────────────────────────────────────────
def stats(vals, label):
    vals = [v for v in vals if v]
    if not vals: return
    vals = [float(v) for v in vals]
    print(f"  {label:20s}: mean={np.mean(vals):.3f}  "
          f"p25={np.percentile(vals,25):.3f}  "
          f"median={np.median(vals):.3f}  "
          f"p75={np.percentile(vals,75):.3f}  "
          f"min={min(vals):.3f}  max={max(vals):.3f}")

passing = [r for r in rows if r.get("avg_iptm") and float(r["avg_iptm"]) >= 0.70]
print(f"\nAll sequences with ipTM ≥ 0.70 (n={len(passing)}):")
stats([r["ipsae_tm8"]    for r in passing], "ipSAE_tm8")
stats([r["min_ipsae_tm8"]for r in passing], "min_ipSAE_tm8")
stats([r["ipsae_frac5"]  for r in passing], "ipSAE_frac5")
stats([r["composite4"]   for r in passing], "composite4")

glmn = [r for r in passing if r["seq_id"].startswith("GLMN")]
rfd  = [r for r in passing if r["seq_id"].startswith("RFD")]
whb  = [r for r in passing if r["seq_id"].startswith("CUL")]
for grp, name in [(glmn,"GLMN"),(rfd,"RFdiffusion"),(whb,"CUL1_WHB")]:
    if grp:
        print(f"\n{name} (n={len(grp)}):")
        stats([r["ipsae_tm8"] for r in grp], "ipSAE_tm8")
        stats([r["composite4"]for r in grp], "composite4")

print("\nTop 20 by composite4 (ipTM ≥ 0.70):")
print(f"{'seq_id':<35} {'ipTM':>6} {'ipLDDT':>7} {'ipSAE_tm8':>10} {'composite4':>11}")
print("-" * 72)
shown = 0
for r in rows:
    if r.get("avg_iptm") and float(r["avg_iptm"]) >= 0.70 and shown < 20:
        print(f"{r['seq_id']:<35} {float(r['avg_iptm']):>6.3f} "
              f"{float(r['avg_iplddt']):>7.3f} "
              f"{float(r['ipsae_tm8']) if r['ipsae_tm8'] else 0:>10.3f} "
              f"{float(r['composite4']) if r['composite4'] else 0:>11.4f}")
        shown += 1
