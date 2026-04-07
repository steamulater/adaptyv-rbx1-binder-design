#!/usr/bin/env python3
"""
Sequence composition audit for final_submission_v2.fasta
Flags: homopolymer runs, dipeptide repeats, charge clusters, low complexity,
       extreme amino acid bias, and other expression/production risk factors.
"""

import re, math
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

FASTA = "/Users/tamukamartin/Desktop/adaptyv_competiton/final_submission_v2.fasta"
BASE  = "/Users/tamukamartin/Desktop/adaptyv_competiton"

# ─────────────────────────────────────────────
# 1. Parse FASTA
# ─────────────────────────────────────────────

def parse_fasta(path):
    records = []
    hdr, seq = None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if hdr: records.append((hdr, ''.join(seq)))
            hdr, seq = line[1:], []
        else:
            seq.append(line.upper())
    if hdr: records.append((hdr, ''.join(seq)))
    return records

records = parse_fasta(FASTA)
print(f"Loaded {len(records)} sequences")

def get_id(hdr):   return hdr.split('|')[0]
def get_scaf(hdr):
    for p in hdr.split('|'):
        if p.startswith('scaffold='): return p.split('=')[1]
    return 'unknown'
def get_rank(hdr):
    for p in hdr.split('|'):
        if p.startswith('rank'): return int(p.replace('rank',''))
    return 999

# ─────────────────────────────────────────────
# 2. Flag functions
# ─────────────────────────────────────────────

def find_homopolymers(seq, min_run=4):
    """Find runs of the same AA >= min_run."""
    hits = []
    for m in re.finditer(r'(.)\1{%d,}' % (min_run-1), seq):
        hits.append({'aa': m.group(1), 'start': m.start(), 'end': m.end(),
                     'len': len(m.group()), 'text': m.group()})
    return hits

def find_dipeptide_repeats(seq, min_reps=3):
    """Find (XY)n repeats for any dipeptide XY, n >= min_reps."""
    hits = []
    for m in re.finditer(r'(..)(\1{%d,})' % (min_reps-1), seq):
        full = m.group(1) + m.group(2)
        hits.append({'motif': m.group(1), 'start': m.start(), 'end': m.end(),
                     'reps': len(full)//2, 'text': full})
    return hits

def charge_cluster(seq, window=10, min_frac=0.7):
    """Windows with >= min_frac charged (K,R,H,D,E) residues."""
    CHARGED = set('KRHDE')
    hits = []
    for i in range(len(seq)-window+1):
        w = seq[i:i+window]
        frac = sum(c in CHARGED for c in w) / window
        charge = sum((1 if c in 'KRH' else -1) for c in w if c in CHARGED)
        if frac >= min_frac:
            hits.append({'start': i, 'end': i+window, 'frac': frac,
                         'net_charge': charge, 'text': w})
    # merge overlapping
    merged = []
    for h in hits:
        if merged and h['start'] <= merged[-1]['end']:
            merged[-1]['end'] = max(merged[-1]['end'], h['end'])
            merged[-1]['text'] = seq[merged[-1]['start']:merged[-1]['end']]
        else:
            merged.append(h.copy())
    return merged

def polyelectrostatic_stretch(seq, window=10, min_net=7):
    """Runs that are almost entirely K/R (pos) or D/E (neg)."""
    hits = []
    for i in range(len(seq)-window+1):
        w = seq[i:i+window]
        pos = sum(c in 'KR' for c in w)
        neg = sum(c in 'DE' for c in w)
        if pos >= min_net:
            hits.append({'type': 'poly-K/R', 'start': i, 'end': i+window, 'text': w, 'count': pos})
        if neg >= min_net:
            hits.append({'type': 'poly-D/E', 'start': i, 'end': i+window, 'text': w, 'count': neg})
    # deduplicate / merge
    merged = {}
    for h in hits:
        key = (h['type'], h['start']//5)
        if key not in merged or h['count'] > merged[key]['count']:
            merged[key] = h
    return list(merged.values())

def low_complexity_score(seq, window=20):
    """Shannon entropy of a window. Low entropy = low complexity."""
    min_entropy = 100
    worst_window = ''
    for i in range(len(seq)-window+1):
        w = seq[i:i+window]
        cnt = Counter(w)
        h = -sum((v/window)*math.log2(v/window) for v in cnt.values())
        if h < min_entropy:
            min_entropy = h
            worst_window = w
    return min_entropy, worst_window

def aa_bias(seq):
    """Return AAs that deviate most from expected proteome frequencies."""
    # Rough SwissProt avg frequencies
    EXPECTED = {'A':8.25,'R':5.53,'N':4.06,'D':5.45,'C':1.37,'Q':3.93,'E':6.75,
                'G':7.07,'H':2.27,'I':5.96,'L':9.66,'K':5.84,'M':2.42,'F':3.86,
                'P':4.70,'S':6.56,'T':5.34,'W':1.08,'Y':2.92,'V':6.87}
    n = len(seq)
    obs = Counter(seq)
    biases = {}
    for aa, exp_pct in EXPECTED.items():
        obs_pct = obs.get(aa, 0) / n * 100
        fold = obs_pct / exp_pct if exp_pct > 0 else 0
        biases[aa] = {'obs': obs_pct, 'exp': exp_pct, 'fold': fold}
    return biases

def net_charge_at_pH7(seq):
    pos = seq.count('K') + seq.count('R') + seq.count('H')*0.1
    neg = seq.count('D') + seq.count('E')
    return pos - neg

def isoelectric_point_approx(seq):
    """Very rough pI estimate."""
    # Simplified: find pH where net charge ≈ 0
    # Just return dominant charge class
    pos = seq.count('K') + seq.count('R')
    neg = seq.count('D') + seq.count('E')
    ratio = (pos+1)/(neg+1)
    if ratio > 2: return 'basic (pI > 8)'
    if ratio < 0.5: return 'acidic (pI < 6)'
    return 'near-neutral'

def has_rare_aas(seq):
    return {'W': seq.count('W'), 'C': seq.count('C'), 'M': seq.count('M')}

# ─────────────────────────────────────────────
# 3. Run audit on all sequences
# ─────────────────────────────────────────────

FLAG_LEVELS = {
    'homopolymer_6+':  'CRITICAL',
    'homopolymer_4-5': 'WARNING',
    'dipeptide_4+':    'WARNING',
    'charge_cluster':  'INFO',
    'low_complexity':  'WARNING',
    'extreme_A_bias':  'CRITICAL',
    'extreme_KE_bias': 'WARNING',
    'no_hydrophobic':  'WARNING',
    'extreme_length':  'WARNING',
}

results = []

for hdr, seq in records:
    sid   = get_id(hdr)
    scaf  = get_scaf(hdr)
    rank  = get_rank(hdr)
    flags = []

    # Homopolymers
    hp = find_homopolymers(seq, min_run=4)
    for h in hp:
        level = 'CRITICAL' if h['len'] >= 6 else 'WARNING'
        flags.append({'type': f'homopolymer_{h["aa"]}x{h["len"]}',
                      'level': level,
                      'detail': f'{h["aa"]}×{h["len"]} at pos {h["start"]+1}-{h["end"]}',
                      'text': h['text']})

    # Dipeptide repeats
    dp = find_dipeptide_repeats(seq, min_reps=3)
    for d in dp:
        flags.append({'type': f'dipeptide_repeat_{d["motif"]}',
                      'level': 'WARNING',
                      'detail': f'({d["motif"]})×{d["reps"]} at pos {d["start"]+1}-{d["end"]}',
                      'text': d['text']})

    # Poly-electrostatic stretches
    pe = polyelectrostatic_stretch(seq, window=10, min_net=7)
    for p in pe:
        flags.append({'type': f'poly_charged_{p["type"]}',
                      'level': 'WARNING',
                      'detail': f'{p["type"]} stretch at pos {p["start"]+1}: {p["text"]}',
                      'text': p['text']})

    # Low complexity
    ent, ww = low_complexity_score(seq, window=20)
    if ent < 2.5:
        level = 'CRITICAL' if ent < 1.8 else 'WARNING'
        flags.append({'type': 'low_complexity',
                      'level': level,
                      'detail': f'entropy={ent:.2f} (worst window: {ww})',
                      'text': ww})

    # AA biases
    bias = aa_bias(seq)
    n = len(seq)
    pct_A = bias['A']['obs']
    pct_E = bias['E']['obs']
    pct_K = bias['K']['obs']
    pct_EK = pct_E + pct_K
    pct_hydrophob = sum(seq.count(aa) for aa in 'VILMFYW') / n * 100

    if pct_A > 25:
        flags.append({'type': 'extreme_A_bias',
                      'level': 'CRITICAL',
                      'detail': f'Alanine = {pct_A:.1f}% (expected ~8.3%)',
                      'text': ''})
    if pct_EK > 40:
        flags.append({'type': 'extreme_EK_bias',
                      'level': 'WARNING',
                      'detail': f'E+K = {pct_EK:.1f}% (expected ~12.6%)',
                      'text': ''})
    if pct_hydrophob < 10:
        flags.append({'type': 'no_hydrophobic_core',
                      'level': 'WARNING',
                      'detail': f'Hydrophobic AA = {pct_hydrophob:.1f}% (expected ~38%)',
                      'text': ''})

    # Length extremes
    if len(seq) > 200:
        flags.append({'type': 'very_long',
                      'level': 'INFO',
                      'detail': f'Length = {len(seq)} AA (>200)',
                      'text': ''})

    # Net charge
    nc = net_charge_at_pH7(seq)
    pI = isoelectric_point_approx(seq)
    rare = has_rare_aas(seq)

    results.append({
        'seq_id': sid, 'scaffold': scaf, 'rank': rank,
        'length': len(seq), 'seq': seq,
        'flags': flags,
        'n_critical': sum(1 for f in flags if f['level']=='CRITICAL'),
        'n_warning':  sum(1 for f in flags if f['level']=='WARNING'),
        'n_info':     sum(1 for f in flags if f['level']=='INFO'),
        'net_charge': nc, 'pI_class': pI,
        'pct_A': pct_A, 'pct_EK': pct_EK, 'pct_hydrophob': pct_hydrophob,
        'entropy_min': ent, 'worst_window': ww,
        'rare_W': rare['W'], 'rare_C': rare['C'],
    })

# ─────────────────────────────────────────────
# 4. Print flagged sequences
# ─────────────────────────────────────────────

critical = [r for r in results if r['n_critical'] > 0]
warned   = [r for r in results if r['n_warning'] > 0 and r['n_critical'] == 0]
clean    = [r for r in results if r['n_critical'] == 0 and r['n_warning'] == 0]

print(f"\n{'='*70}")
print(f"SEQUENCE COMPOSITION AUDIT — {len(results)} submitted designs")
print(f"{'='*70}")
print(f"  CRITICAL issues:  {len(critical)} designs")
print(f"  WARNING only:     {len(warned)} designs")
print(f"  Clean:            {len(clean)} designs")

print(f"\n{'─'*70}")
print("CRITICAL FLAGS")
print(f"{'─'*70}")
for r in sorted(critical, key=lambda x: x['rank']):
    print(f"\n[rank {r['rank']:>3}] {r['seq_id']}  ({r['scaffold']}, {r['length']} AA)")
    print(f"  Sequence: {r['seq']}")
    for f in r['flags']:
        marker = '🚨' if f['level']=='CRITICAL' else '⚠️ '
        print(f"  {marker} {f['type']}: {f['detail']}")

print(f"\n{'─'*70}")
print("WARNING FLAGS (no critical)")
print(f"{'─'*70}")
for r in sorted(warned, key=lambda x: x['rank']):
    print(f"\n[rank {r['rank']:>3}] {r['seq_id']}  ({r['scaffold']}, {r['length']} AA)")
    for f in r['flags']:
        print(f"  ⚠️  {f['type']}: {f['detail']}")

print(f"\n{'─'*70}")
print(f"CLEAN ({len(clean)} designs — no critical or warning flags)")
print(f"{'─'*70}")
for r in sorted(clean, key=lambda x: x['rank']):
    print(f"  [rank {r['rank']:>3}] {r['seq_id']}  ({r['scaffold']}, {r['length']} AA, "
          f"E+K={r['pct_EK']:.0f}%, A={r['pct_A']:.0f}%, hydrophob={r['pct_hydrophob']:.0f}%)")

# ─────────────────────────────────────────────
# 5. Per-scaffold breakdown
# ─────────────────────────────────────────────
print(f"\n{'='*70}")
print("PER-SCAFFOLD ISSUE BREAKDOWN")
print(f"{'='*70}")
from collections import defaultdict
scaf_stats = defaultdict(lambda: {'n':0,'critical':0,'warning':0,'clean':0})
for r in results:
    s = r['scaffold']
    scaf_stats[s]['n'] += 1
    if r['n_critical']: scaf_stats[s]['critical'] += 1
    elif r['n_warning']: scaf_stats[s]['warning'] += 1
    else: scaf_stats[s]['clean'] += 1

for s, st in sorted(scaf_stats.items()):
    print(f"  {s:>15}: {st['n']} total | "
          f"CRITICAL={st['critical']} | WARNING={st['warning']} | CLEAN={st['clean']}")

# ─────────────────────────────────────────────
# 6. Figure
# ─────────────────────────────────────────────

fig = plt.figure(figsize=(18, 14))
gs  = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)
fig.suptitle('Sequence Composition Audit — 100 Submitted RBX1 Binders',
             fontsize=14, fontweight='bold')

scaffolds   = ['GLMN','CUL1_WHB','RFdiffusion']
scaf_colors = {'GLMN':'#2196F3','CUL1_WHB':'#FF9800','RFdiffusion':'#4CAF50'}

# rank→color for scatter
def rank_color(r):
    if r <= 20: return '#E53935'  # top 20 = red
    if r <= 50: return '#FF9800'  # 21-50 = orange
    return '#9E9E9E'              # 51-100 = grey

all_seqs  = [r['seq'] for r in results]
all_ranks = [r['rank'] for r in results]

# ── panel A: %A vs rank ──
ax = fig.add_subplot(gs[0, 0])
ax.scatter([r['rank'] for r in results],
           [r['pct_A'] for r in results],
           c=[rank_color(r['rank']) for r in results], s=25, alpha=0.8)
ax.axhline(25, color='red', linestyle='--', lw=1.2, label='25% threshold (CRITICAL)')
ax.axhline(8.3, color='grey', linestyle=':', lw=1, label='SwissProt avg')
ax.set_xlabel('Submission rank'); ax.set_ylabel('% Alanine')
ax.set_title('Alanine content by rank', fontsize=10)
ax.legend(fontsize=7)
# label worst offenders
for r in results:
    if r['pct_A'] > 22:
        ax.annotate(r['seq_id'].split('_')[0]+'_'+r['seq_id'].split('_')[1],
                    (r['rank'], r['pct_A']), fontsize=6, ha='left', va='bottom')

# ── panel B: %E+K vs rank ──
ax2 = fig.add_subplot(gs[0, 1])
ax2.scatter([r['rank'] for r in results],
            [r['pct_EK'] for r in results],
            c=[rank_color(r['rank']) for r in results], s=25, alpha=0.8)
ax2.axhline(40, color='red', linestyle='--', lw=1.2, label='40% threshold')
ax2.axhline(12.6, color='grey', linestyle=':', lw=1, label='SwissProt avg')
ax2.set_xlabel('Submission rank'); ax2.set_ylabel('% E + K')
ax2.set_title('Glutamate+Lysine content by rank', fontsize=10)
ax2.legend(fontsize=7)

# ── panel C: minimum entropy (low complexity) ──
ax3 = fig.add_subplot(gs[0, 2])
ax3.scatter([r['rank'] for r in results],
            [r['entropy_min'] for r in results],
            c=[rank_color(r['rank']) for r in results], s=25, alpha=0.8)
ax3.axhline(1.8, color='red', linestyle='--', lw=1.2, label='CRITICAL <1.8 bits')
ax3.axhline(2.5, color='orange', linestyle='--', lw=1.2, label='WARNING <2.5 bits')
ax3.set_xlabel('Submission rank'); ax3.set_ylabel('Min window entropy (bits)')
ax3.set_title('Sequence complexity (Shannon entropy)', fontsize=10)
ax3.set_ylim(0, 5)
ax3.legend(fontsize=7)

# ── panel D: %hydrophobic ──
ax4 = fig.add_subplot(gs[1, 0])
ax4.scatter([r['rank'] for r in results],
            [r['pct_hydrophob'] for r in results],
            c=[rank_color(r['rank']) for r in results], s=25, alpha=0.8)
ax4.axhline(10, color='red', linestyle='--', lw=1.2, label='<10% WARNING')
ax4.axhline(38, color='grey', linestyle=':', lw=1, label='SwissProt avg')
ax4.set_xlabel('Submission rank'); ax4.set_ylabel('% hydrophobic (V,I,L,M,F,Y,W)')
ax4.set_title('Hydrophobic residue content', fontsize=10)
ax4.legend(fontsize=7)

# ── panel E: net charge ──
ax5 = fig.add_subplot(gs[1, 1])
ax5.scatter([r['rank'] for r in results],
            [r['net_charge'] for r in results],
            c=[rank_color(r['rank']) for r in results], s=25, alpha=0.8)
ax5.axhline(0, color='black', lw=0.8)
ax5.set_xlabel('Submission rank'); ax5.set_ylabel('Net charge at pH 7')
ax5.set_title('Net charge (pos=K/R/H, neg=D/E)', fontsize=10)
# color band
ax5.axhline(20, color='red', linestyle=':', lw=0.8, label='|charge|>20 risky')
ax5.axhline(-20, color='red', linestyle=':', lw=0.8)
ax5.legend(fontsize=7)

# ── panel F: flag count per design (stacked bar, top 50) ──
ax6 = fig.add_subplot(gs[1, 2])
top50 = sorted(results, key=lambda x: x['rank'])[:50]
xs = [r['rank'] for r in top50]
crits   = [r['n_critical'] for r in top50]
warns   = [r['n_warning']  for r in top50]
infos   = [r['n_info']     for r in top50]
ax6.bar(xs, crits, color='#E53935', label='Critical', width=0.8)
ax6.bar(xs, warns, bottom=crits, color='#FF9800', label='Warning', width=0.8)
ax6.bar(xs, infos, bottom=[c+w for c,w in zip(crits,warns)], color='#90A4AE', label='Info', width=0.8)
ax6.set_xlabel('Submission rank (top 50)'); ax6.set_ylabel('Flag count')
ax6.set_title('Issue flags — top 50 designs', fontsize=10)
ax6.legend(fontsize=7)

# ── panel G: homopolymer run length heatmap (residues, top 50) ──
ax7 = fig.add_subplot(gs[2, 0])
# for each top-50 design, compute max homopolymer run
max_runs = []
run_aas  = []
for r in sorted(results, key=lambda x: x['rank'])[:50]:
    hp = find_homopolymers(r['seq'], min_run=1)
    if hp:
        best = max(hp, key=lambda h: h['len'])
        max_runs.append(best['len'])
        run_aas.append(best['aa'])
    else:
        max_runs.append(1)
        run_aas.append('-')
colors_hp = ['#E53935' if l>=6 else '#FF9800' if l>=4 else '#4CAF50' for l in max_runs]
xs50 = list(range(1,51))
ax7.bar(xs50, max_runs, color=colors_hp, width=0.8)
ax7.axhline(6, color='red', linestyle='--', lw=1.2, label='≥6: CRITICAL')
ax7.axhline(4, color='orange', linestyle='--', lw=1.2, label='≥4: WARNING')
for i, (x, l, aa) in enumerate(zip(xs50, max_runs, run_aas)):
    if l >= 4:
        ax7.text(x, l+0.1, aa, ha='center', va='bottom', fontsize=5)
ax7.set_xlabel('Submission rank (top 50)'); ax7.set_ylabel('Max homopolymer run (AA)')
ax7.set_title('Longest homopolymer run', fontsize=10)
ax7.legend(fontsize=7)

# ── panel H: AA composition heatmap (top 20 designs) ──
ax8 = fig.add_subplot(gs[2, 1:])
AAS = list('ACDEFGHIKLMNPQRSTVWY')
top20 = sorted(results, key=lambda x: x['rank'])[:20]
mat = np.zeros((len(top20), len(AAS)))
for i, r in enumerate(top20):
    n = len(r['seq'])
    for j, aa in enumerate(AAS):
        mat[i, j] = r['seq'].count(aa) / n * 100
labels_y = [f"#{r['rank']} {r['seq_id'][:12]}" for r in top20]
im = ax8.imshow(mat, aspect='auto', cmap='YlOrRd', vmin=0, vmax=35)
ax8.set_xticks(range(len(AAS))); ax8.set_xticklabels(list(AAS), fontsize=8)
ax8.set_yticks(range(len(top20))); ax8.set_yticklabels(labels_y, fontsize=7)
ax8.set_title('AA composition heatmap — top 20 designs (% per residue)', fontsize=10)
plt.colorbar(im, ax=ax8, fraction=0.02, pad=0.02, label='%')

# legend for rank colors
patches = [mpatches.Patch(color='#E53935', label='Rank 1-20'),
           mpatches.Patch(color='#FF9800', label='Rank 21-50'),
           mpatches.Patch(color='#9E9E9E', label='Rank 51-100')]
fig.legend(handles=patches, loc='lower left', fontsize=8, ncol=3,
           bbox_to_anchor=(0.01, 0.0))

out_png = f"{BASE}/sequence_composition_audit.png"
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\nSaved figure: {out_png}")

# ─────────────────────────────────────────────
# 7. Summary stats table
# ─────────────────────────────────────────────
print(f"\n{'='*70}")
print("OVERALL STATISTICS")
print(f"{'='*70}")
pct_As       = [r['pct_A'] for r in results]
pct_EKs      = [r['pct_EK'] for r in results]
entropies    = [r['entropy_min'] for r in results]
hydrophobs   = [r['pct_hydrophob'] for r in results]
net_charges  = [r['net_charge'] for r in results]

# count any design with homopolymer >=4
n_hp4  = sum(1 for r in results if any('homopolymer' in f['type'] for f in r['flags']))
n_hp6  = sum(1 for r in results if any('homopolymer' in f['type'] and int(re.search(r'x(\d+)',f['type']).group(1)) >= 6
                                        for f in r['flags'] if re.search(r'x(\d+)',f['type'])))
n_dp   = sum(1 for r in results if any('dipeptide' in f['type'] for f in r['flags']))
n_lc   = sum(1 for r in results if any('low_complexity' in f['type'] for f in r['flags']))
n_hA   = sum(1 for r in results if r['pct_A'] > 25)
n_hEK  = sum(1 for r in results if r['pct_EK'] > 40)
n_nohyd = sum(1 for r in results if r['pct_hydrophob'] < 10)

print(f"  Designs with homopolymer run ≥ 4:  {n_hp4}/100")
print(f"  Designs with homopolymer run ≥ 6:  {n_hp6}/100  (CRITICAL)")
print(f"  Designs with dipeptide repeats:    {n_dp}/100")
print(f"  Designs with low complexity:       {n_lc}/100")
print(f"  Designs with >25% Alanine:         {n_hA}/100  (CRITICAL)")
print(f"  Designs with >40% E+K:             {n_hEK}/100")
print(f"  Designs with <10% hydrophobic:     {n_nohyd}/100")
print()
print(f"  %A     — mean: {np.mean(pct_As):.1f}%, max: {np.max(pct_As):.1f}%  "
      f"({results[np.argmax(pct_As)]['seq_id']})")
print(f"  %E+K   — mean: {np.mean(pct_EKs):.1f}%, max: {np.max(pct_EKs):.1f}%  "
      f"({results[np.argmax(pct_EKs)]['seq_id']})")
print(f"  Entropy — mean: {np.mean(entropies):.2f}, min: {np.min(entropies):.2f}  "
      f"({results[np.argmin(entropies)]['seq_id']})")
print(f"  Hydrophob — mean: {np.mean(hydrophobs):.1f}%, min: {np.min(hydrophobs):.1f}%  "
      f"({results[np.argmin(hydrophobs)]['seq_id']})")
print(f"  Net charge — range: {min(net_charges):.0f} to {max(net_charges):.0f}")
print()
print("Top 5 most problematic designs:")
results_sorted_by_issues = sorted(results, key=lambda x: (x['n_critical']*10 + x['n_warning']), reverse=True)
for r in results_sorted_by_issues[:5]:
    print(f"  [rank {r['rank']:>3}] {r['seq_id']} — {r['n_critical']} CRITICAL, {r['n_warning']} WARNINGS")
    for f in r['flags']:
        print(f"           {f['level']}: {f['type']}: {f['detail']}")

plt.show()
