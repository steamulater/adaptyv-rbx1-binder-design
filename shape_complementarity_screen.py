#!/usr/bin/env python3
"""
Screen RFdiffusion backbones for interface geometry quality.

Uses Boltz-2 predicted complex structures (CIF, chain A = binder, chain B = RBX1).
Computes three geometric metrics from Cα/Cβ coordinates — no Rosetta needed:

  n_contacts      : number of inter-chain residue pairs with Cβ distance < 8 Å
                    Proxy for interface size. < 8 = likely broken/no interface.

  contact_density : n_contacts / sqrt(len_A * len_B)
                    Normalised for binder length. Good binders: > 0.15.

  patch_spread    : std of hotspot-side contact positions along 3 PCA axes.
                    Low std = all contacts concentrated at one point (bad).
                    High std = contacts spread across the binding surface (good).

  sc_proxy        : combined score = contact_density * patch_spread_norm
                    Range 0–1, higher = better shape complementarity.

Filters:
  PASS   : n_contacts >= 10  AND  contact_density >= 0.12  AND  sc_proxy >= 0.10
  BORDERLINE : n_contacts >= 6  AND  sc_proxy >= 0.06
  FAIL   : anything else (broken geometry, no real interface)

Output:
  shape_screen_results.csv  — all 151 sequences with metrics + pass/fail
  rerun_mpnn_candidates.txt — seq_ids of PASS sequences (prioritise for more ProteinMPNN)
  rerun_mpnn_backbones.txt  — mapping seq_id → backbone CIF path
"""

import os, glob, csv
import numpy as np

COMP_DIR = "boltz_rfdiffusion/complex_results/boltz_results_complex/predictions"
CB_CUTOFF = 8.0   # Å — residues within this are "in contact"
MIN_CONTACTS = 8  # hard minimum — fewer = broken geometry


# ── CIF parser ────────────────────────────────────────────────────────────────
def parse_cif_coords(cif_path):
    """Return dict: chain -> list of (resnum, x, y, z) for CA atoms (or CB for non-Gly)."""
    cols = {}
    coords = {"A": [], "B": []}
    in_loop = False

    with open(cif_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("_atom_site."):
                field = line.split(".")[-1].strip()
                cols[field] = len(cols)
                in_loop = True
                continue
            if not in_loop:
                continue
            if line.startswith("ATOM"):
                parts = line.split()
                try:
                    atom_name = parts[cols["label_atom_id"]]
                    chain     = parts[cols["label_asym_id"]]
                    resnum    = int(parts[cols["label_seq_id"]])
                    x         = float(parts[cols["Cartn_x"]])
                    y         = float(parts[cols["Cartn_y"]])
                    z         = float(parts[cols["Cartn_z"]])
                    res_name  = parts[cols["label_comp_id"]]
                except (KeyError, ValueError, IndexError):
                    continue

                # Use CB (Cβ) for all residues except Gly; fall back to CA
                if chain in coords:
                    if atom_name == "CB":
                        # Replace any existing CA entry for this residue
                        coords[chain] = [(r, *xyz) for r, *xyz in coords[chain]
                                         if r != resnum]
                        coords[chain].append((resnum, x, y, z))
                    elif atom_name == "CA":
                        # Only add if CB not already present
                        if not any(r == resnum for r, *_ in coords[chain]):
                            coords[chain].append((resnum, x, y, z))
    return coords


def compute_interface_metrics(coords):
    """Compute contact count, density, patch spread, and sc_proxy."""
    a = np.array([[x, y, z] for _, x, y, z in sorted(coords["A"])])
    b = np.array([[x, y, z] for _, x, y, z in sorted(coords["B"])])

    if len(a) == 0 or len(b) == 0:
        return None

    # Pairwise distances (broadcasting)
    diff   = a[:, None, :] - b[None, :, :]          # (nA, nB, 3)
    dists  = np.sqrt((diff ** 2).sum(axis=-1))       # (nA, nB)

    contact_mask = dists < CB_CUTOFF
    n_contacts   = int(contact_mask.sum())

    if n_contacts == 0:
        return {"n_contacts": 0, "contact_density": 0.0,
                "patch_spread": 0.0, "sc_proxy": 0.0}

    nA, nB       = len(a), len(b)
    contact_density = n_contacts / np.sqrt(nA * nB)

    # Patch spread: PCA of contact midpoints
    ai, bi = np.where(contact_mask)
    midpoints = (a[ai] + b[bi]) / 2.0               # (n_contacts, 3)

    if len(midpoints) >= 3:
        centred  = midpoints - midpoints.mean(axis=0)
        cov      = np.cov(centred.T)
        eigvals  = np.linalg.eigvalsh(cov)
        eigvals  = np.sort(eigvals)[::-1]
        # Spread = geometric mean of top 2 principal components (Å)
        patch_spread = float(np.sqrt(max(eigvals[0], 0) * max(eigvals[1], 0)))
    else:
        patch_spread = 0.0

    # Normalise patch_spread to 0–1 using a sigmoid-like mapping (centre ~8 Å)
    patch_spread_norm = 1 / (1 + np.exp(-(patch_spread - 8) / 4))

    sc_proxy = contact_density * patch_spread_norm

    return {
        "n_contacts":      n_contacts,
        "contact_density": float(contact_density),
        "patch_spread":    float(patch_spread),
        "sc_proxy":        float(sc_proxy),
    }


# ── Main loop ─────────────────────────────────────────────────────────────────
print("Scanning RFdiffusion complex structures...")
results = []

for pred_dir in sorted(glob.glob(f"{COMP_DIR}/RFD_*_complex")):
    seq_id = os.path.basename(pred_dir).replace("_complex", "")

    # Average metrics across available models (seeds)
    cif_files = sorted(glob.glob(f"{pred_dir}/*.cif"))
    if not cif_files:
        continue

    seed_metrics = []
    for cif in cif_files:
        try:
            coords = parse_cif_coords(cif)
            m = compute_interface_metrics(coords)
            if m:
                seed_metrics.append(m)
        except Exception:
            continue

    if not seed_metrics:
        continue

    # Average across seeds
    avg = {k: float(np.mean([m[k] for m in seed_metrics])) for k in seed_metrics[0]}

    # Classify
    if avg["n_contacts"] >= 10 and avg["contact_density"] >= 0.12 and avg["sc_proxy"] >= 0.10:
        flag = "PASS"
    elif avg["n_contacts"] >= 6 and avg["sc_proxy"] >= 0.06:
        flag = "BORDERLINE"
    else:
        flag = "FAIL"

    results.append({"seq_id": seq_id, "flag": flag, **avg})
    print(f"  {seq_id:<30} contacts={avg['n_contacts']:5.1f}  "
          f"density={avg['contact_density']:.3f}  "
          f"spread={avg['patch_spread']:5.1f}Å  "
          f"sc={avg['sc_proxy']:.3f}  [{flag}]")

# ── Summary ───────────────────────────────────────────────────────────────────
results.sort(key=lambda r: -r["sc_proxy"])
pass_     = [r for r in results if r["flag"] == "PASS"]
border    = [r for r in results if r["flag"] == "BORDERLINE"]
fail      = [r for r in results if r["flag"] == "FAIL"]

print(f"\n{'='*60}")
print(f"Total screened : {len(results)}")
print(f"PASS           : {len(pass_)}   (good geometry → more ProteinMPNN)")
print(f"BORDERLINE     : {len(border)}  (marginal — include if slots available)")
print(f"FAIL           : {len(fail)}   (broken geometry → skip)")

print(f"\nTop 20 by sc_proxy:")
print(f"{'seq_id':<30} {'contacts':>8} {'density':>8} {'spread':>8} {'sc_proxy':>9}")
print("-" * 68)
for r in results[:20]:
    print(f"{r['seq_id']:<30} {r['n_contacts']:>8.1f} {r['contact_density']:>8.3f} "
          f"{r['patch_spread']:>8.1f} {r['sc_proxy']:>9.3f}  [{r['flag']}]")

print(f"\nFAIL sequences (broken geometry):")
for r in fail:
    print(f"  {r['seq_id']:<30} contacts={r['n_contacts']:.1f}  sc={r['sc_proxy']:.3f}")

# ── Save outputs ──────────────────────────────────────────────────────────────
with open("shape_screen_results.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["seq_id", "flag", "n_contacts",
                                       "contact_density", "patch_spread", "sc_proxy"])
    w.writeheader()
    w.writerows(results)
print(f"\nSaved shape_screen_results.csv")

# Candidates for more ProteinMPNN — PASS + BORDERLINE, sorted by sc_proxy
rerun = pass_ + border
with open("rerun_mpnn_candidates.txt", "w") as f:
    for r in rerun:
        f.write(f"{r['seq_id']}\t{r['flag']}\t{r['sc_proxy']:.4f}\n")
print(f"Saved rerun_mpnn_candidates.txt ({len(rerun)} sequences)")

# Backbone CIF paths for re-running ProteinMPNN
with open("rerun_mpnn_backbones.txt", "w") as f:
    f.write("# seq_id\tsc_proxy\tflag\tbest_model_cif\n")
    for r in rerun:
        pred_dir = f"{COMP_DIR}/{r['seq_id']}_complex"
        cif_files = sorted(glob.glob(f"{pred_dir}/*.cif"))
        if cif_files:
            f.write(f"{r['seq_id']}\t{r['sc_proxy']:.4f}\t{r['flag']}\t{cif_files[0]}\n")
print(f"Saved rerun_mpnn_backbones.txt")
