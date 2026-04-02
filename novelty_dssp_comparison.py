#!/usr/bin/env python3
"""
Novelty + DSSP secondary structure comparison across all binder design methods.
Methods: GLMN, CUL1-WHB, RFdiffusion, LatentX

Outputs:
  - novelty_dssp_summary.csv   — per-method aggregated stats
  - novelty_dssp_per_design.csv — per-design stats
  - novelty_dssp_comparison.png — summary figure
"""

import os, glob, io, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
import pydssp

BASE = "/Users/tamukamartin/Desktop/adaptyv_competiton"

# ─────────────────────────────────────────────
# 1. Load sequence + binding metric data
# ─────────────────────────────────────────────

res = pd.read_csv(f"{BASE}/rescored_v2.csv")
res['seq_id_clean'] = res['seq_id'].str.replace('_complex', '')

nov1 = pd.read_csv(f"{BASE}/novelty_screen_results.csv")  # GLMN + CUL1-WHB
nov2 = pd.read_csv(f"{BASE}/batch2_novelty_results.csv")  # RFdiffusion

# Merge novelty into res
nov_all = pd.concat([
    nov1[['seq_id','best_pident','novel']].rename(columns={'seq_id':'seq_id_clean'}),
    nov2[['seq_id','best_pident','novel']].rename(columns={'seq_id':'seq_id_clean'}),
])
# Normalise novel column to bool (nov1 has True/False, nov2 has "YES"/"NO")
nov_all['novel'] = nov_all['novel'].map(lambda x: True if x in (True, 'YES', 'yes') else (False if x in (False, 'NO', 'no') else np.nan))
res = res.merge(nov_all, on='seq_id_clean', how='left')

lt = pd.read_csv(f"{BASE}/LatentX/latentx_scored.csv")
lt['scaffold'] = 'LatentX'
lt['seq_id_clean'] = lt['name']
lt['best_pident'] = np.nan   # no SwissProt screen for LatentX locally
lt['novel'] = np.nan
# rename cols for compatibility
lt = lt.rename(columns={'iptm':'avg_iptm', 'iplddt':'avg_iplddt'})

print(f"Prior designs: {len(res)}  (GLMN={len(res[res.scaffold=='GLMN'])}, "
      f"CUL1_WHB={len(res[res.scaffold=='CUL1_WHB'])}, "
      f"RFdiffusion={len(res[res.scaffold=='RFdiffusion'])})")
print(f"LatentX designs: {len(lt)}")

# ─────────────────────────────────────────────
# 2. DSSP helper functions
# ─────────────────────────────────────────────

def cif_to_pdb_text(cif_path, chain_id=None):
    """Convert CIF → PDB text, optionally selecting one chain."""
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", cif_path)

    class ChainSelect(Select):
        def accept_chain(self, chain):
            return chain_id is None or chain.id == chain_id

    pdbio = PDBIO()
    pdbio.set_structure(struct)
    buf = io.StringIO()
    pdbio.save(buf, ChainSelect())
    return buf.getvalue()


def pdb_text_to_dssp(pdb_text):
    """Run pydssp on PDB text string → (ss_array, helix%, sheet%, coil%)."""
    try:
        coord = pydssp.read_pdbtext(pdb_text)
        if coord.shape[0] < 5:
            return None
        ss = pydssp.assign(coord)  # array of 'H','E','-'
        n = len(ss)
        h = np.sum(ss == 'H') / n
        e = np.sum(ss == 'E') / n
        c = np.sum(ss == '-') / n
        return {'ss': ss, 'helix': h, 'sheet': e, 'coil': c, 'length': n}
    except Exception as ex:
        warnings.warn(f"DSSP failed: {ex}")
        return None


def pdb_file_to_dssp(pdb_path, chain_id=None):
    """PDB file → DSSP stats, selecting chain if specified."""
    if chain_id is not None:
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("s", pdb_path)
        class ChainSel(Select):
            def accept_chain(self, c): return c.id == chain_id
        pdbio = PDBIO(); pdbio.set_structure(struct)
        buf = io.StringIO(); pdbio.save(buf, ChainSel())
        pdb_text = buf.getvalue()
    else:
        pdb_text = open(pdb_path).read()
    return pdb_text_to_dssp(pdb_text)


# ─────────────────────────────────────────────
# 3. Compute DSSP per design
# ─────────────────────────────────────────────

records = []

# ── 3a. GLMN + CUL1-WHB: PDB complexes, binder = chain A ──
preds_dir = f"{BASE}/boltz_DESIGN_results/boltz_results_complexes/predictions"

for _, row in res[res.scaffold.isin(['GLMN','CUL1_WHB'])].iterrows():
    sid = row['seq_id']  # e.g. "GLMN_T0.1_s4_complex" or "CUL1_WHB_T0.1_s1_complex"
    design_dir = f"{preds_dir}/{sid}"
    if not os.path.isdir(design_dir):
        # try without _complex
        sid_bare = sid.replace('_complex','') + '_complex'
        design_dir = f"{preds_dir}/{sid_bare}"
    # take model_0
    pdbs = sorted(glob.glob(f"{design_dir}/*_model_0.pdb"))
    if not pdbs:
        continue
    ds = pdb_file_to_dssp(pdbs[0], chain_id='A')
    if ds is None:
        continue
    records.append({
        'seq_id': row['seq_id_clean'],
        'scaffold': row['scaffold'],
        'sequence': row.get('sequence',''),
        'best_pident': row.get('best_pident', np.nan),
        'novel': row.get('novel', np.nan),
        'avg_iptm': row.get('avg_iptm', np.nan),
        'avg_iplddt': row.get('avg_iplddt', np.nan),
        'seq_len': len(str(row.get('sequence',''))),
        'helix': ds['helix'],
        'sheet': ds['sheet'],
        'coil': ds['coil'],
        'struct_len': ds['length'],
    })

print(f"GLMN+CUL1_WHB processed: {len(records)}")

# ── 3b. RFdiffusion: CIF complexes, binder = chain A ──
rfd_preds = f"{BASE}/boltz_rfdiffusion/complex_results/boltz_results_complex/predictions"
rfd_designs = res[res.scaffold == 'RFdiffusion']
rfd_count_before = len(records)

for _, row in rfd_designs.iterrows():
    sid = row['seq_id']  # e.g. "RFD_160_best"
    design_dir = f"{rfd_preds}/{sid}_complex"
    if not os.path.isdir(design_dir):
        continue
    cifs = sorted(glob.glob(f"{design_dir}/*_model_0.cif"))
    if not cifs:
        continue
    pdb_text = cif_to_pdb_text(cifs[0], chain_id='A')
    ds = pdb_text_to_dssp(pdb_text)
    if ds is None:
        continue
    records.append({
        'seq_id': row['seq_id_clean'],
        'scaffold': 'RFdiffusion',
        'sequence': row.get('sequence',''),
        'best_pident': row.get('best_pident', np.nan),
        'novel': row.get('novel', np.nan),
        'avg_iptm': row.get('avg_iptm', np.nan),
        'avg_iplddt': row.get('avg_iplddt', np.nan),
        'seq_len': len(str(row.get('sequence',''))),
        'helix': ds['helix'],
        'sheet': ds['sheet'],
        'coil': ds['coil'],
        'struct_len': ds['length'],
    })

print(f"RFdiffusion processed: {len(records)-rfd_count_before}")

# ── 3c. LatentX: monomer CIFs (chain A = binder) ──
ltx_cifs = sorted(glob.glob(f"{BASE}/LatentX/spellbound_carter/*.cif") +
                  glob.glob(f"{BASE}/LatentX/supreme_korsakov/*.cif"))
ltx_scored = lt.set_index('name')
ltx_count_before = len(records)

for cif_path in ltx_cifs:
    name = os.path.basename(cif_path).replace('.cif','')
    if name not in ltx_scored.index:
        continue
    row = ltx_scored.loc[name]
    pdb_text = cif_to_pdb_text(cif_path, chain_id='A')
    ds = pdb_text_to_dssp(pdb_text)
    if ds is None:
        continue
    seq = str(row.get('seq',''))
    records.append({
        'seq_id': name,
        'scaffold': 'LatentX',
        'sequence': seq,
        'best_pident': np.nan,
        'novel': np.nan,
        'avg_iptm': row.get('avg_iptm', np.nan),
        'avg_iplddt': row.get('avg_iplddt', np.nan),
        'seq_len': len(seq),
        'helix': ds['helix'],
        'sheet': ds['sheet'],
        'coil': ds['coil'],
        'struct_len': ds['length'],
    })

print(f"LatentX processed: {len(records)-ltx_count_before}")

df = pd.DataFrame(records)
print(f"\nTotal designs with DSSP: {len(df)}")
print(df.groupby('scaffold').size())

# ─────────────────────────────────────────────
# 4. Sequence novelty: pairwise self-similarity
# ─────────────────────────────────────────────
# Compute within-method pairwise identity (fraction identical / min_length)
# as a proxy for within-method diversity

from itertools import combinations

def pairwise_identity(seqs):
    """Mean pairwise identity (Hamming over aligned positions) for a list of seqs."""
    if len(seqs) < 2:
        return np.nan
    ids = []
    for s1, s2 in combinations(seqs, 2):
        l = min(len(s1), len(s2))
        if l == 0: continue
        matches = sum(a == b for a, b in zip(s1[:l], s2[:l]))
        ids.append(matches / l)
    return np.mean(ids) * 100  # return as percent

# ─────────────────────────────────────────────
# 5. Amino acid composition
# ─────────────────────────────────────────────

def aa_composition(seqs):
    all_seq = ''.join(seqs)
    n = len(all_seq)
    if n == 0: return {}
    charged_pos = sum(all_seq.count(aa) for aa in 'KRH')
    charged_neg = sum(all_seq.count(aa) for aa in 'DE')
    hydrophobic  = sum(all_seq.count(aa) for aa in 'VILMFYW')
    polar        = sum(all_seq.count(aa) for aa in 'STNQ')
    return {
        'pct_charged_pos': charged_pos/n*100,
        'pct_charged_neg': charged_neg/n*100,
        'pct_hydrophobic': hydrophobic/n*100,
        'pct_polar': polar/n*100,
    }

# ─────────────────────────────────────────────
# 6. Per-method summary
# ─────────────────────────────────────────────

method_order = ['GLMN','CUL1_WHB','RFdiffusion','LatentX']
method_colors = {'GLMN':'#2196F3','CUL1_WHB':'#FF9800','RFdiffusion':'#4CAF50','LatentX':'#9C27B0'}
method_labels = {'GLMN':'GLMN','CUL1_WHB':'CUL1-WHB','RFdiffusion':'RFdiffusion','LatentX':'LatentX'}

summary_rows = []
for method in method_order:
    mdf = df[df.scaffold == method]
    if len(mdf) == 0: continue
    seqs = [s for s in mdf['sequence'].tolist() if isinstance(s, str) and len(s) > 0]
    aac = aa_composition(seqs)
    row = {
        'Method': method_labels[method],
        'N_designs': len(mdf),
        # Sequence novelty
        'SwissProt_pident_mean': mdf['best_pident'].mean(),
        'SwissProt_pident_std': mdf['best_pident'].std(),
        'Novel_pct': (mdf['novel'] == True).mean() * 100 if mdf['novel'].notna().any() else np.nan,
        'Within_method_identity_pct': pairwise_identity(seqs),
        # Length
        'seq_len_mean': mdf['seq_len'].mean(),
        'seq_len_std': mdf['seq_len'].std(),
        # DSSP
        'helix_pct': mdf['helix'].mean() * 100,
        'sheet_pct': mdf['sheet'].mean() * 100,
        'coil_pct':  mdf['coil'].mean() * 100,
        # Binding
        'avg_iptm': mdf['avg_iptm'].mean(),
        'avg_iplddt': mdf['avg_iplddt'].mean(),
        # AA composition
        **{k: v for k,v in aac.items()},
    }
    summary_rows.append(row)

summary = pd.DataFrame(summary_rows)
print("\n=== Method Comparison Summary ===")
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
pd.set_option('display.float_format', '{:.1f}'.format)
print(summary.to_string(index=False))

df.to_csv(f"{BASE}/novelty_dssp_per_design.csv", index=False)
summary.to_csv(f"{BASE}/novelty_dssp_summary.csv", index=False)
print(f"\nSaved: novelty_dssp_per_design.csv ({len(df)} rows)")
print(f"Saved: novelty_dssp_summary.csv ({len(summary)} rows)")

# ─────────────────────────────────────────────
# 7. Figure
# ─────────────────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle('Binder Design Method Comparison: Sequence Novelty, Structure & Composition',
             fontsize=14, fontweight='bold', y=1.01)

methods = [r['Method'] for r in summary_rows]
colors  = [method_colors[m] for m in [r for r in method_order if r in df.scaffold.unique()]]

# ── 7a. SwissProt % identity (novelty) ──
ax = axes[0, 0]
vals  = summary['SwissProt_pident_mean'].values
errs  = summary['SwissProt_pident_std'].values
valid = ~np.isnan(vals)
bars = ax.bar([m for m,v in zip(methods, valid) if v],
              vals[valid], yerr=errs[valid],
              color=[c for c,v in zip(colors, valid) if v],
              capsize=4, alpha=0.85, edgecolor='black', linewidth=0.5)
ax.axhline(30, color='red', linestyle='--', linewidth=1, label='<30% = novel')
ax.set_title('SwissProt Best-Hit %Identity\n(lower = more novel)', fontsize=10)
ax.set_ylabel('% Identity to SwissProt')
ax.set_ylim(0, 60)
ax.legend(fontsize=8)
ax.bar_label(bars, fmt='%.1f%%', fontsize=8, padding=2)
for i, (m, v) in enumerate(zip(methods, valid)):
    if not v:
        ax.text(i-0.15*len(methods), 5, 'N/A', fontsize=9, color='gray')

# ── 7b. Within-method pairwise sequence diversity ──
ax = axes[0, 1]
vals2 = summary['Within_method_identity_pct'].values
valid2 = ~np.isnan(vals2)
bars2 = ax.bar(methods, vals2,
               color=colors, alpha=0.85, edgecolor='black', linewidth=0.5)
ax.set_title('Within-Method Pairwise Seq Identity\n(lower = more diverse)', fontsize=10)
ax.set_ylabel('Mean pairwise %identity')
ax.bar_label(bars2, fmt='%.1f%%', fontsize=8, padding=2)
ax.set_ylim(0, 100)

# ── 7c. Sequence length distribution ──
ax = axes[0, 2]
for method, color in zip(method_order, colors):
    mdf = df[df.scaffold == method]
    if len(mdf) == 0: continue
    ax.scatter([method_labels[method]]*len(mdf), mdf['seq_len'],
               color=color, alpha=0.4, s=20)
vals3 = summary['seq_len_mean'].values
ax.scatter(methods, vals3, color=colors, s=80, zorder=5, edgecolors='black', linewidths=0.8)
ax.set_title('Sequence Length Distribution', fontsize=10)
ax.set_ylabel('Sequence length (AA)')

# ── 7d. DSSP stacked bar ──
ax = axes[1, 0]
helix = summary['helix_pct'].values
sheet = summary['sheet_pct'].values
coil  = summary['coil_pct'].values
x = np.arange(len(methods))
w = 0.55
b1 = ax.bar(x, helix, w, label='Helix (H)', color='#E53935', alpha=0.85)
b2 = ax.bar(x, sheet, w, bottom=helix, label='Sheet (E)', color='#1E88E5', alpha=0.85)
b3 = ax.bar(x, coil,  w, bottom=helix+sheet, label='Coil (-)', color='#43A047', alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(methods)
ax.set_ylim(0, 110)
ax.set_title('Secondary Structure Composition\n(Boltz-2 predicted binder structures)', fontsize=10)
ax.set_ylabel('% of residues')
ax.legend(fontsize=8, loc='upper right')
for xi, (h,s,c) in enumerate(zip(helix, sheet, coil)):
    ax.text(xi, h/2, f'{h:.0f}%', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    ax.text(xi, h+s/2, f'{s:.0f}%', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    ax.text(xi, h+s+c/2, f'{c:.0f}%', ha='center', va='center', fontsize=8, color='white', fontweight='bold')

# ── 7e. AA class composition ──
ax = axes[1, 1]
aa_classes = ['pct_charged_pos','pct_charged_neg','pct_hydrophobic','pct_polar']
aa_labels  = ['Pos charged\n(K,R,H)', 'Neg charged\n(D,E)', 'Hydrophobic\n(V,I,L,M,F,Y,W)', 'Polar\n(S,T,N,Q)']
aa_colors  = ['#E91E63','#FF5722','#795548','#00BCD4']
x = np.arange(len(methods))
w = 0.18
for i, (cls, lbl, col) in enumerate(zip(aa_classes, aa_labels, aa_colors)):
    vals_aa = [r.get(cls, np.nan) for r in summary_rows]
    ax.bar(x + i*w - 1.5*w, vals_aa, w, label=lbl, color=col, alpha=0.85,
           edgecolor='black', linewidth=0.3)
ax.set_xticks(x); ax.set_xticklabels(methods)
ax.set_title('Amino Acid Class Composition', fontsize=10)
ax.set_ylabel('% of residues')
ax.legend(fontsize=7, loc='upper right')

# ── 7f. Binding quality vs novelty scatter ──
ax = axes[1, 2]
for method, color in zip(method_order, colors):
    mdf = df[df.scaffold == method]
    if len(mdf) == 0: continue
    ax.scatter(mdf['best_pident'], mdf['avg_iptm'],
               color=color, alpha=0.6, s=25, label=method_labels[method],
               edgecolors='none')
ax.axvline(30, color='red', linestyle='--', linewidth=1, alpha=0.5)
ax.set_xlabel('SwissProt best-hit %identity (lower = novel)')
ax.set_ylabel('Boltz-2 ipTM')
ax.set_title('Novelty vs Binding Quality', fontsize=10)
ax.legend(fontsize=8)

plt.tight_layout()
out_png = f"{BASE}/novelty_dssp_comparison.png"
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\nSaved figure: {out_png}")
plt.show()

# ─────────────────────────────────────────────
# 8. Print clean comparison table
# ─────────────────────────────────────────────
print("\n" + "="*100)
print("COMPLETE METHOD COMPARISON TABLE")
print("="*100)
header = (f"{'Method':>14}  {'N':>3}  {'SP%id':>6}  {'%Novel':>7}  "
          f"{'IntraDiv%':>9}  {'SeqLen':>6}  "
          f"{'%Helix':>6}  {'%Sheet':>6}  {'%Coil':>6}  "
          f"{'ipTM':>5}  {'ipLDDT':>7}")
print(header)
print("-"*100)
for r in summary_rows:
    sp_str = f"{r['SwissProt_pident_mean']:6.1f}" if not np.isnan(r.get('SwissProt_pident_mean', np.nan)) else "   N/A"
    nv_str = f"{r['Novel_pct']:7.1f}" if not np.isnan(r.get('Novel_pct', np.nan)) else "    N/A"
    print(f"{r['Method']:>14}  {int(r['N_designs']):>3}  {sp_str}  {nv_str}  "
          f"{r['Within_method_identity_pct']:9.1f}  {r['seq_len_mean']:6.1f}  "
          f"{r['helix_pct']:6.1f}  {r['sheet_pct']:6.1f}  {r['coil_pct']:6.1f}  "
          f"{r['avg_iptm']:5.3f}  {r['avg_iplddt']:7.3f}")
print("="*100)
