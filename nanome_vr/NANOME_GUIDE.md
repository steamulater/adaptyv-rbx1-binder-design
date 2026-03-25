# Nanome VR — RBX1 Binder Design Session Guide

## Files in this folder

| File | What it is | Why it matters |
|---|---|---|
| `00_RBX1_target_RING_domain.pdb` | RBX1 RING domain (chain A, 77 AA) | The target — show this first |
| `01_RFD167_top_denovo_rank1.cif` | RFD_167 + RBX1 complex, best model (ipTM 0.905) | #1 ranked design overall |
| `02_RFD114_top_denovo_rank2.cif` | RFD_114 + RBX1 complex (ipTM 0.883) | #2 ranked, 85 AA helix bundle |
| `03_RFD160_highest_iptm_0910.cif` | RFD_160 + RBX1 complex (ipTM 0.913) | Highest raw ipTM in entire submission |
| `04_RFD97_denovo_helix.cif` | RFD_97 + RBX1 complex (ipTM 0.887) | Clean compact helical binder |
| `05_RFD38_denovo_compact.cif` | RFD_38 + RBX1 complex (ipTM 0.911) | High ipLDDT (0.816) — well-defined interface |
| `06_RFD52_denovo.cif` | RFD_52 + RBX1 complex (ipTM 0.857) | |
| `07_GLMN_s4_scaffold_top.pdb` | GLMN redesign s4 + RBX1 (ipTM 0.918) | Best GLMN sequence — natural binding mode |
| `08_GLMN_s11_scaffold.pdb` | GLMN redesign s11 + RBX1 (ipTM 0.901) | Temperature diversity example |
| `09_CUL1WHB_induced_fit_example.pdb` | CUL1_WHB + RBX1 (ipTM 0.761) | Induced-fit failure mode — RBX1 deforms |

Chain convention in all complexes:
- **Chain A** = binder (designed protein)
- **Chain B** = RBX1 RING domain

---

## Nanome Setup Checklist

1. Load files via **File > Open** or drag into workspace
2. For side-by-side comparisons: load two structures, use **Align** (Structure panel) to superpose on chain B (RBX1)
3. Hotspot residues to highlight (RBX1 chain B): `4, 6, 12, 14, 15, 23, 24, 26, 28, 52, 56, 60, 64, 66`
4. Recommended colour scheme:
   - RBX1 (chain B): **grey surface** — the target
   - Binder (chain A) RFdiffusion: **green cartoon/ribbon**
   - Binder (chain A) GLMN: **blue cartoon/ribbon**
   - Binder (chain A) CUL1_WHB: **orange cartoon/ribbon**
   - Hotspot residues: **red sticks**
5. Interface view: select residues within 5Å of chain A on chain B → **Show sticks + H-bonds**

---

## Scene Sequence for Recording

### Scene 1 — "Meet the target"
Load `00_RBX1_target_RING_domain.pdb`
- Show as surface, colour by hydrophobicity
- Rotate to show the concave RING-H2 face
- Highlight the 14 hotspot residues in red sticks
- **Narration:** "This is RBX1 — the master E3 ligase subunit involved in ~20% of all protein degradation in your cells. We're designing proteins that block this exact surface."

### Scene 2 — "The natural blueprint (GLMN scaffold)"
Load `07_GLMN_s4_scaffold_top.pdb`
- GLMN (chain A) = blue ribbon; RBX1 = grey surface
- Zoom into interface, show H-bond contacts
- **Narration:** "Strategy 1: we started from a natural protein — Glomulin — that already occupies this surface in cells. ProteinMPNN redesigned its sequence to create 48 novel binders. Every single one predicted to bind."

### Scene 3 — "Teaching AI to design from scratch (RFdiffusion)"
Load `01_RFD167_top_denovo_rank1.cif` alongside `07_GLMN_s4_scaffold_top.pdb`
- Superpose on RBX1 (chain B)
- RFdiffusion binder = green; GLMN = blue; RBX1 = grey
- **Narration:** "Strategy 2: we gave RFdiffusion nothing but the target surface and let it hallucinate a completely new protein. Different shape, different approach angle — but same surface. This is our #1 ranked design."

### Scene 4 — "What the AI confidence score looks like"
Load `03_RFD160_highest_iptm_0910.cif`
- Colour by pLDDT (B-factor column) — blue=confident, red=uncertain
- Zoom into interface region
- **Narration:** "Boltz-2 predicted this complex and gave it an ipTM of 0.91. Colour = confidence. The interface is deep blue — the AI is very sure these two proteins grip each other exactly like this."

### Scene 5 — "When a design fails — induced fit"
Load `09_CUL1WHB_induced_fit_example.pdb`
- Show RBX1 deformation vs `00_RBX1_target_RING_domain.pdb` (superpose)
- **Narration:** "Not everything worked. This compact 72 AA binder predicted to bind — but only by deforming RBX1's structure by 5 Å. That's a red flag. Good binders should grip the target as-is, not reshape it."

### Scene 6 — "De novo vs scaffold — interface comparison"
Load `01_RFD167_top_denovo_rank1.cif` and `07_GLMN_s4_scaffold_top.pdb` superposed
- Split-screen or sequential zoom into each interface
- Measure interface area, count contacts
- **Narration:** "Same target. Completely different solutions. GLMN: 247 amino acids, natural binding mode. RFD_167: 70 amino acids, AI-hallucinated, zero similarity to any known protein. This is what modern protein design looks like."

---

## TikTok / Reels Series Plan

### Series title: "I designed 100 proteins in a week — here's what it looked like"

---

**Episode 1 — Hook**
*"I used AI to design 100 proteins to fight cancer. Here's the target."*
- 15-30s: fly into RBX1 surface in Nanome, reveal hotspot patch in red
- Text overlay: "RBX1 — master switch for 20% of protein degradation in your cells"
- End hook: "We needed something that grips this exact surface. Spoiler: AI did it better than nature expected."

---

**Episode 2 — The scaffold strategy**
*"We started by copying nature's homework"*
- Show Glomulin wrapping around RBX1 in Nanome (file 07)
- "This protein — Glomulin — already does this in cells. We kept the shape, redesigned every amino acid."
- Show ProteinMPNN text on screen: "48 new sequences. 48/48 predicted to bind."
- Tools mentioned: **ProteinMPNN**, **Boltz-2**, **PyMOL**

---

**Episode 3 — RFdiffusion hallucination**
*"Then we asked AI to invent a protein from scratch"*
- Screen recording of RFdiffusion running (or terminal output scrolling)
- Cut to Nanome: fly through the RFD_167 binder docking into RBX1
- "200 backbones hallucinated. 0% similarity to any known protein. This one [zoom in] scores higher than anything nature gave us."
- Tools mentioned: **RFdiffusion**, **Colab A100**

---

**Episode 4 — What ipTM actually means**
*"Here's how we know a design will work (before testing it)"*
- In Nanome: colour complex by pLDDT (confidence)
- Show PAE matrix side by side (screenshot from analysis)
- "Blue = confident. Red = uncertain. When the interface is all blue — the AI thinks these two proteins lock together exactly like this."
- Tools mentioned: **Boltz-2**, **ipTM**, **ipSAE**, **ipLDDT**

---

**Episode 5 — The failure**
*"This one looked great on paper — until we checked the fine print"*
- Load `09_CUL1WHB_induced_fit_example.pdb`
- Show RBX1 RMSD animation (superpose native vs deformed)
- "ipTM: 0.76 ✓. But RBX1 moved 5 Å to accommodate it. That's like a lock that only fits your key if you bend the doorframe."
- "We submitted 7 of these anyway — as a diversity bet."
- Tools mentioned: **RMSD analysis**, **PyMOL**

---

**Episode 6 — The Nipah retrospective**
*"We cheated — we looked at 1,000 proteins that were already tested"*
- Show nipah_analysis.png on screen
- "Before finalising our 100 sequences, we analysed results from a previous competition: 1,030 proteins, real experimental binding data."
- "Turns out the metric we were using (ipTM) is only the 7th best predictor. This one — complex ipLDDT — is 15 points better."
- "We re-scored everything. Recovered 26 sequences we had incorrectly filtered out."
- Tools mentioned: **retrospective analysis**, **AUROC**, **Adaptyv Nipah dataset**

---

**Episode 7 — The final 100**
*"100 sequences. 3 strategies. 1 submission."*
- Nanome: fly through all 6 RFdiffusion structures loaded together, orbiting RBX1
- Scorecard text overlay:
  - "47 de novo (RFdiffusion) — 0% identity to any known protein"
  - "46 scaffold redesign (GLMN) — 100% pass rate"
  - "7 compact binders (CUL1-WHB) — diversity bet"
- "We find out if any of them actually bind in 4-6 weeks."
- "Follow to see the experimental results."

---

## Nanome Tips for Recording

- **Screen record** the VR headset view via Nanome's built-in stream/mirror to your Mac
- Use **Presentation Mode** in Nanome to hide UI chrome for clean shots
- **Slow rotation** (joystick at low speed) films better than fast movement
- For interface zooms: get into the binding pocket, pause, rotate slowly
- **Colour by chain** for establishing shots; **colour by B-factor/pLDDT** for confidence shots
- Record at least 60s of each scene — you'll cut down to 10-15s per clip

## Hashtags
`#proteindesign #AIbiology #nanome #rfddiffusion #boltz2 #proteinnengineering #syntheticbiology #computationalbiology #drugdiscovery #cancer`
