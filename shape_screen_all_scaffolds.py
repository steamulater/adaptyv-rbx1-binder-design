#!/usr/bin/env python3
"""
Shape complementarity screen across ALL scaffolds:
  - GLMN redesigns (Batch 1, PDB, 247 AA binder)
  - CUL1_WHB designs (Batch 1, PDB, 72-80 AA binder)
  - RFdiffusion de novo (Batch 2, CIF, 65-95 AA binder)

Sanity check: GLMN (evolved scaffold) should score highest,
RFdiffusion middle, CUL1_WHB lowest.
"""

import os, glob, csv
import numpy as np

B1_COMP  = "boltz_DESIGN_results/boltz_results_complexes/predictions"
B2_COMP  = "boltz_rfdiffusion/complex_results/boltz_results_complex/predictions"
CB_CUT   = 8.0

# ── Parsers ───────────────────────────────────────────────────────────────────
def parse_pdb_coords(pdb_path):
    """Parse PDB: return {chain: [(resnum, x, y, z)]} for CB (or CA for Gly)."""
    coords = {}
    seen_cb = set()
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom = line[12:16].strip()
            chain = line[21]
            try:
                resnum = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            if chain not in coords:
                coords[chain] = []
            key = (chain, resnum)
            if atom == "CB":
                coords[chain] = [r for r in coords[chain] if not (r[0] == resnum)]
                coords[chain].append((resnum, x, y, z))
                seen_cb.add(key)
            elif atom == "CA" and key not in seen_cb:
                coords[chain].append((resnum, x, y, z))
    return coords


def parse_cif_coords(cif_path):
    """Parse mmCIF: return {chain: [(resnum, x, y, z)]} using label_asym_id."""
    cols = {}
    coords = {}
    seen_cb = set()
    in_loop = False
    with open(cif_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("_atom_site."):
                field = line.split(".")[-1].strip()
                cols[field] = len(cols)
                in_loop = True
                continue
            if not in_loop or not line.startswith("ATOM"):
                continue
            parts = line.split()
            try:
                atom  = parts[cols["label_atom_id"]]
                chain = parts[cols["label_asym_id"]]
                resnum= int(parts[cols["label_seq_id"]])
                x     = float(parts[cols["Cartn_x"]])
                y     = float(parts[cols["Cartn_y"]])
                z     = float(parts[cols["Cartn_z"]])
            except (KeyError, ValueError, IndexError):
                continue
            if chain not in coords:
                coords[chain] = []
            key = (chain, resnum)
            if atom == "CB":
                coords[chain] = [r for r in coords[chain] if r[0] != resnum]
                coords[chain].append((resnum, x, y, z))
                seen_cb.add(key)
            elif atom == "CA" and key not in seen_cb:
                coords[chain].append((resnum, x, y, z))
    return coords


def get_ab(coords):
    """Return chain A and B arrays, sorted by resnum."""
    def arr(ch):
        return np.array([[x, y, z] for _, x, y, z in sorted(coords.get(ch, []))])
    return arr("A"), arr("B")


def interface_metrics(a, b):
    if len(a) == 0 or len(b) == 0:
        return None
    dists = np.sqrt(((a[:, None] - b[None]) ** 2).sum(-1))
    mask  = dists < CB_CUT
    n     = int(mask.sum())
    if n == 0:
        return dict(n_contacts=0, contact_density=0.0, patch_spread=0.0, sc_proxy=0.0)

    density = n / np.sqrt(len(a) * len(b))
    ai, bi  = np.where(mask)
    mids    = (a[ai] + b[bi]) / 2.0
    spread  = 0.0
    if len(mids) >= 3:
        c   = mids - mids.mean(0)
        ev  = np.linalg.eigvalsh(np.cov(c.T))[::-1]
        spread = float(np.sqrt(max(ev[0], 0) * max(ev[1], 0)))
    spread_norm = 1 / (1 + np.exp(-(spread - 8) / 4))
    return dict(n_contacts=n, contact_density=float(density),
                patch_spread=float(spread), sc_proxy=float(density * spread_norm))


def process_dir(pred_dir, fmt):
    """Average metrics across all models in a prediction directory."""
    pattern = "*.cif" if fmt == "cif" else "*.pdb"
    files   = sorted(glob.glob(f"{pred_dir}/{pattern}"))
    # exclude confidence json etc — only structure files
    files   = [f for f in files if not os.path.basename(f).startswith("confidence")]
    if not files:
        return None
    seed_metrics = []
    for fp in files:
        try:
            coords = parse_cif_coords(fp) if fmt == "cif" else parse_pdb_coords(fp)
            a, b   = get_ab(coords)
            m      = interface_metrics(a, b)
            if m:
                seed_metrics.append(m)
        except Exception:
            continue
    if not seed_metrics:
        return None
    return {k: float(np.mean([m[k] for m in seed_metrics])) for k in seed_metrics[0]}


# ── Run all scaffolds ─────────────────────────────────────────────────────────
results = []

print("Processing GLMN redesigns (Batch 1)...")
for d in sorted(glob.glob(f"{B1_COMP}/GLMN_*_complex")):
    sid = os.path.basename(d).replace("_complex", "")
    m   = process_dir(d, "pdb")
    if m:
        results.append({"seq_id": sid, "scaffold": "GLMN", **m})

print("Processing CUL1_WHB designs (Batch 1)...")
for d in sorted(glob.glob(f"{B1_COMP}/CUL1_*_complex")):
    sid = os.path.basename(d).replace("_complex", "")
    m   = process_dir(d, "pdb")
    if m:
        results.append({"seq_id": sid, "scaffold": "CUL1_WHB", **m})

print("Processing RFdiffusion de novo (Batch 2)...")
for d in sorted(glob.glob(f"{B2_COMP}/RFD_*_complex")):
    sid = os.path.basename(d).replace("_complex", "")
    m   = process_dir(d, "cif")
    if m:
        results.append({"seq_id": sid, "scaffold": "RFdiffusion", **m})

# ── Summary stats per scaffold ────────────────────────────────────────────────
print("\n" + "="*65)
for scaffold in ["GLMN", "CUL1_WHB", "RFdiffusion"]:
    grp = [r for r in results if r["scaffold"] == scaffold]
    sc  = [r["sc_proxy"]       for r in grp]
    ct  = [r["n_contacts"]     for r in grp]
    cd  = [r["contact_density"]for r in grp]
    ps  = [r["patch_spread"]   for r in grp]
    print(f"\n{scaffold} (n={len(grp)}):")
    print(f"  sc_proxy      : mean={np.mean(sc):.3f}  median={np.median(sc):.3f}  "
          f"p25={np.percentile(sc,25):.3f}  p75={np.percentile(sc,75):.3f}  "
          f"min={min(sc):.3f}  max={max(sc):.3f}")
    print(f"  n_contacts    : mean={np.mean(ct):.1f}  median={np.median(ct):.1f}  "
          f"min={min(ct):.1f}  max={max(ct):.1f}")
    print(f"  contact_density: mean={np.mean(cd):.3f}  median={np.median(cd):.3f}")
    print(f"  patch_spread  : mean={np.mean(ps):.1f}Å  median={np.median(ps):.1f}Å")

# ── Cross-scaffold comparison table ──────────────────────────────────────────
print("\n\nPer-sequence comparison (all scaffolds, sorted by sc_proxy):")
print(f"{'seq_id':<35} {'scaffold':<13} {'contacts':>9} {'density':>8} {'spread':>8} {'sc_proxy':>9}")
print("-"*87)
for r in sorted(results, key=lambda x: -x["sc_proxy"]):
    print(f"{r['seq_id']:<35} {r['scaffold']:<13} {r['n_contacts']:>9.1f} "
          f"{r['contact_density']:>8.3f} {r['patch_spread']:>8.1f} {r['sc_proxy']:>9.3f}")

# ── Save ──────────────────────────────────────────────────────────────────────
with open("shape_screen_all.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["seq_id","scaffold","n_contacts",
                                       "contact_density","patch_spread","sc_proxy"])
    w.writeheader()
    w.writerows(sorted(results, key=lambda x: -x["sc_proxy"]))
print(f"\nSaved shape_screen_all.csv ({len(results)} sequences)")
