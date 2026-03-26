# Lab Notebook: ProteWorks Studio — BoltGen Design Sprint
**Platform:** ProteWorks Studio (studio.proteworks.com) — BoltGen pipeline
**Competition:** GEM x Adaptyv Bio — RBX1 Binder Design Challenge
**Objective:** Generate novel RBX1 binder candidates using BoltGen; concurrently stress-test the ProteWorks platform and log structured feedback for Brandon (founder)
**Experimenter:** Tamuka Martin Chidyausiku @steamulater
**Notebook started:** March 25, 2026

---

## Context

This is a parallel design sprint running alongside (and after) the primary RFdiffusion + ProteinMPNN pipeline documented in `../lab_notebook.md`. That pipeline produced 100 sequences submitted to the competition on March 25, 2026.

This notebook documents a second wave of designs generated entirely within ProteWorks Studio, using their BoltGen workflow. Two goals run in parallel:

1. **Science:** Produce high-quality RBX1 binder candidates for future submission rounds or follow-up experimental work. Apply everything learned from the Nipah retrospective (composite scoring, ipSAE weighting, RING RMSD screening) to evaluate outputs rigorously.

2. **Platform evaluation:** Stress-test ProteWorks Studio end-to-end — input handling, job management, output quality, UI/UX friction — and compile actionable feedback for Brandon. Treat this like a structured beta test with a real design problem.

---

## Target: RBX1 (same as primary campaign)

**UniProt:** P62877 | **PDB:** 2LGV | **Organism:** *Homo sapiens*

**Full sequence (108 AA):**
```
MAAAMDVDTPSGTNSGAGKKRFEVKKWNAVALWAWDIVVDNCAICRNHIMDLCIECQANQASATSEECTVAWGVCNHAFHFHCISRWLKTRQVCPLDNREWEFQKYGH
```

**RING domain (target region, residues 32–90):**
```
NCAICRNHIMDLCIECQANQASATSEECTVAWGVCNHAFHFHCISRWLKTRQVCPLDNREWEFQKYGH
```

**Hotspot residues (renumbered RING domain):** `4, 6, 12, 14, 15, 23, 24, 26, 28, 52, 56, 60, 64, 66`

**Therapeutic rationale:** RBX1 is the catalytic RING subunit of Cullin-RING E3 ligase complexes, directing ubiquitination of ~20% of all cellular proteins. Blocking the E2-binding surface disrupts oncogenic CRL activity. Natural inhibitor: Glomulin (GLMN), which occupies the RING-H2 face with high affinity.

**Key structural files (in `../nanome_vr/`):**
- `00_RBX1_target_RING_domain.pdb` — isolated RING domain
- `01_RFD167_top_denovo_rank1.cif` — best de novo binder (ipTM 0.905)
- `07_GLMN_s4_scaffold_top.pdb` — best GLMN scaffold design (ipTM 0.918)

---

## Scoring Framework (from Nipah retrospective)

Metrics ranked by AUROC on 1,030 experimentally validated binders:

| Metric | AUROC | Use |
|---|---|---|
| complex_ipLDDT | 0.691 | Primary confidence metric |
| min_ipSAE | 0.638 | Interface tightness |
| ipTM | 0.603 | Structural compatibility |
| mono pTM | 0.501 | Not used (random) |

**Composite score:** `composite = 0.4 × ipTM + 0.3 × ipLDDT + 0.3 × norm_ipSAE`

**Hard filters:**
- ipTM >= 0.70
- No poly-Ala or low-complexity sequences
- RING RMSD < 3.0 Å (RBX1 must not deform to accommodate binder)

**Soft thresholds (relative, not absolute):**
- ipLDDT: rank within cohort, do not apply Nipah's 0.85 cutoff
- ipSAE: prefer lower Å; Nipah's 0-1 scale is not directly comparable

---

## Platform: ProteWorks Studio

**URL:** https://studio.proteworks.com/
**Tool:** BoltGen
**Contact:** Brandon (founder)

### What BoltGen is (working understanding)
BoltGen is a web-based protein design tool built on top of Boltz-2. Unlike the Colab-based pipeline used in the primary campaign (RFdiffusion → ProteinMPNN → Boltz-2), BoltGen appears to integrate structure prediction and sequence design into a single interface. To be confirmed through use.

### Evaluation criteria for platform feedback

Track the following for each session:

| Category | Things to note |
|---|---|
| **Input UX** | How is target specified? PDB upload, UniProt ID, sequence? Hotspot specification workflow |
| **Job submission** | Queue times, job naming, batch submission capability |
| **Output quality** | Metrics returned (ipTM? pLDDT? PAE? ipSAE?), structure download format |
| **Output UX** | Result dashboard, sorting/filtering, export options |
| **Reproducibility** | Same inputs → same outputs? Seed control? |
| **Failure modes** | What breaks, what errors look like, recovery path |
| **Speed** | Wall time per job, throughput |
| **Gaps vs Colab pipeline** | What can't be done that we needed before |

---

## Design Plan

### Round 1 — Baseline replication
**Goal:** Reproduce a known-good result from the primary campaign to calibrate BoltGen output quality.
- Input: RBX1 RING domain + hotspot residues (same as primary campaign)
- Design a small set (5–10 sequences) with default BoltGen parameters
- Score with same composite metric
- Compare ipTM/ipLDDT distribution vs primary campaign baseline

### Round 2 — Parameter exploration
**Goal:** Explore BoltGen-specific settings (if any) to understand design space coverage.
- Vary whatever parameters BoltGen exposes (temperature, number of sequences, etc.)
- Document what is and isn't controllable vs Colab pipeline

### Round 3 — Novel design angles
**Goal:** Use BoltGen for design directions not explored in primary campaign.
- Shorter binders (< 60 AA) — if BoltGen supports constrained length
- Beta-sheet dominated scaffolds — primary campaign was helix-heavy
- Mixed alpha/beta if platform generates diverse topologies

---

## Entries

---

### Entry 001 — Platform First Access
**Date:** 2026-03-25
**Status:** COMPLETE
**Screenshots:** `1.png` → `5.png`

#### Session goals
- [x] Create account / log in
- [x] Explore dashboard and available tools
- [x] Resolve WebGL error and confirm 3D viewer working
- [x] Document billing/compute structure
- [ ] Submit first BoltGen design job (next session)

---

#### Step-by-step session log

**Step 1 — First load (1.png)**
Landed on `studio.proteworks.com`. Platform is currently in **public beta** — a prominent banner reads:
> "PROTEIN STUDIO IS IN BETA: BUGS & TEMPORARY UI ELEMENTS MAY OCCUR"

This is important context for the whole evaluation: we are testing a beta product.

Immediately hit the **WebGL error** in the 3D sequence viewer (centre panel). The rest of the UI loaded correctly — left sidebar (Workflows, File Browser), right panel (Protein Tools, Load Files, Save Files) all functional. Error is isolated to the WebGL canvas.

**Step 2 — Fix applied**
Navigated to `chrome://settings/system` → enabled **"Use hardware acceleration when available"** → relaunched Chrome. WebGL error resolved. See BUG-001 in feedback log for full details.

**Step 3 — Onboarding / subscription (2.png)**
After login, platform redirected through an onboarding subscription flow. Key details observed:

- **Trial is active** — free tier with compute credits
- **No payment method on file** — platform warns that workflows will stop when credits are exhausted and intermediate data may be deleted at trial end
- **Subscription period:** ends April 25, 2026
- **Seats:** 1 purchased, 1 in use
- Orange warning banner is prominent — good UX choice to surface this early

**Step 4 — Compute pricing (3.png)**
Billing page reveals compute cost rates:

| Resource | Unit | Rate (USD/hr) |
|---|---|---|
| CPU | vCPU | $0.05 |
| RAM | GB | $0.0083 |
| NVIDIA T4 | GPU | $0.45 |
| NVIDIA A100 | GPU | $3.72 |

Billing history for March 2026: **$0.00** — trial credits in use, no charges yet.

**Step 5 — Credits overview (4.png)**
Organisation: **STEAMulater**

| Metric | Value |
|---|---|
| Total compute cost | $0.00 |
| Compute credits | $500.00 |
| Outstanding balance | $0.00 |
| Remaining soft budget | $500.00 |
| Remaining hard budget | $500.00 |
| Auto-billing threshold | $50.00 |

$500 in trial credits is generous — enough to run many BoltGen jobs on T4 GPUs (~1,100 GPU hours at $0.45/hr) before needing to add a payment method.

**Step 6 — 3D viewer confirmed working (5.png)**
After the WebGL fix, loaded **2lgv** (RBX1) directly into the Sequence Viewer. Structure rendered correctly as teal cartoon ribbon. Sequence bar visible at top. STRUCTURES panel on the right lists `2lgv`. This confirms:
- PDB ID input works (direct lookup, no upload needed for known structures)
- 3D viewer is functional post-fix
- RBX1 is accessible by its PDB code

---

#### Platform observations (first impressions)

**What works well:**
- Clean, uncluttered UI layout — sequence viewer centre, tools right, workflows left
- Direct PDB ID loading (no manual upload required for known structures)
- Billing transparency — costs shown upfront before any jobs are run
- $500 trial credits gives real runway for genuine design work

**Open questions for next session:**
- Where is BoltGen accessed? Is it under "Run Workflow"?
- What inputs does BoltGen require — PDB upload, chain selection, hotspot residues?
- How many sequences per job? Batch submission?
- What metrics are returned — ipTM, ipLDDT, PAE matrices?
- Can we download structure files (CIF/PDB) for downstream analysis?

#### Feedback for Brandon
*(Log friction points, bugs, missing features, things that worked well)*

| # | Category | Observation | Severity |
|---|---|---|---|
| 1 | Bug | WebGL unavailable on Chrome — blocks 3D viewer entirely. Fixed by enabling hardware acceleration. See BUG-001 in feedback log. | High |
| 2 | UX improvement | WebGL error message should include a direct link to `chrome://settings/system` with one-sentence fix instructions. Current message links to caniuse.com which is not actionable for most users. | Medium |
| 3 | UX | Beta banner is appropriately prominent — good call to set expectations early | Positive |
| 4 | Billing UX | Orange warning about data deletion at trial end is clear and well-placed in the onboarding flow | Positive |
| 5 | Feature | No payment method warning could also suggest "save your data locally before trial ends" — data loss risk is non-obvious | Low |

---

## Platform Feedback Log (running)

This section accumulates all feedback across sessions for eventual handoff to Brandon.

| Session | Category | Observation | Severity (Low/Med/High) | Suggested fix |
|---|---|---|---|---|
| 001 | Browser compatibility | WebGL unavailable error on Chrome: "WebGL does not seem to be available. This can be caused by an outdated browser, graphics card driver issue, or bad weather." Blocks access to the 3D viewer entirely. | High | **RESOLVED:** Enable hardware acceleration — `chrome://settings/system` → "Use hardware acceleration when available" → Relaunch |

#### BUG-001: WebGL unavailable on Chrome

**Error message:**
```
WebGL does not seem to be available.
This can be caused by an outdated browser, graphics card driver issue, or bad weather.
Sometimes, just restarting the browser helps. Also, make sure hardware acceleration
is enabled in your browser.
For a list of supported browsers, refer to http://caniuse.com/#feat=webgl.
```

**Root cause:** Chrome has disabled hardware acceleration or GPU compositing, either due to a driver blocklist entry, a crashed GPU process, or the setting being manually turned off.

**Fix checklist (try in order):**

1. **Enable hardware acceleration in Chrome**
   - Go to `chrome://settings/system`
   - Toggle **"Use hardware acceleration when available"** ON
   - Click **Relaunch**

2. **Force-enable GPU / override blocklist**
   - Go to `chrome://flags`
   - Search `Override software rendering list` → set to **Enabled**
   - Search `WebGL Developer Extensions` → set to **Enabled**
   - Click **Relaunch**

3. **Check GPU status**
   - Go to `chrome://gpu`
   - Look for `WebGL: Hardware accelerated` — if it says "Software only, hardware acceleration unavailable", the driver is the issue
   - Update your graphics card driver (NVIDIA GeForce Experience / AMD Adrenalin / macOS Software Update)

4. **macOS specific**
   - System Settings → Battery → uncheck **"Enable Power Nap"** and ensure the machine is not in low-power mode (throttles GPU)
   - If on a MacBook with discrete + integrated GPU, ensure Chrome is using the discrete GPU (Energy tab in Activity Monitor)

5. **Try a clean Chrome profile**
   - Open `chrome://version` and note profile path
   - Launch Chrome with `--disable-gpu-sandbox` flag as a test:
     ```
     /Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome --disable-gpu-sandbox
     ```

6. **Alternative browsers**
   - Firefox and Safari both support WebGL by default on macOS — try these as a workaround while debugging Chrome
   - Edge (Chromium-based) typically inherits Chrome's WebGL support but with fewer blocklist entries

**Recommended platform-side fix for Brandon:**
- Add a pre-load WebGL detection check that surfaces a clearer error with a direct link to the Chrome hardware acceleration setting (`chrome://settings/system`) and one-sentence instructions: *"Go to Chrome Settings > System > enable Use hardware acceleration > Relaunch"*
- The GPU override flag (`chrome://flags/#override-software-rendering-list`) as a secondary option
- Consider falling back to a non-WebGL viewer (e.g., minimal 2D sequence representation) so users aren't fully blocked from the rest of the platform while they fix their browser
- **Confirmed fix:** enabling hardware acceleration in Chrome settings resolves the issue immediately (tested 2026-03-25)

---

### Entry 002 — BoltGen Job Setup & First Submission
**Date:** 2026-03-25
**Status:** JOB RUNNING
**Screenshots:** `Screenshot 2026-03-25 at 9.18.35 PM.png` → `Screenshot 2026-03-25 at 9.30.29 PM.png`

---

#### Step-by-step session log

**Step 1 — Finding BoltGen (9.18.35)**
Clicked "Run Workflow" in the left sidebar. The Workflows page lists available pipelines. Only one visible: **BoltGen Design Pipeline** — described as "Generate protein designs using BoltGen with parallel GPU processing and filtering." Design Specification panel offers:
- Add PROTEIN
- Add RNA
- Add DNA
- Add LIGAND
- Add Structure File
- Load From File (import a saved spec)

Note: at this point the org dropdown still showed "No Organization" — had to switch to STEAMulater org first.

**Step 2 — Org switch + persistent warning (9.21.18)**
Switched to STEAMulater. Two persistent banners now visible on every page:
1. `"No payment method set remaining compute credits: $500.00. When credits are exhausted, running workflows will be stopped."`
2. `"SUBSCRIPTION WILL NOT RENEW: DOWNLOAD YOUR DATA OR RESUBSCRIBE BEFORE ACCESS ENDS"`

Banner #2 is a new red alert that wasn't present in Entry 001 — likely triggered by completing the billing onboarding flow without adding a payment method.

**Step 3 — Design Specification inputs (9.21.51)**
Scrolled through the full BoltGen form. Complete parameter set documented below.

**Step 4 — Protocol selection (9.22.10)**
Protocol dropdown options:
- Protein-Anything *(selected)*
- Peptide-Anything
- Protein-Small Molecule
- Nanobody-Anything
- Antibody-Anything

Selected **Protein-Anything** — most general; appropriate for de novo binder design with no scaffold constraints.

**Step 5 — Hardware + output config (9.22.57)**
GPU Type dropdown:
- NVIDIA Tesla T4
- **NVIDIA Tesla A100** *(selected — "High performance, faster for large structures")*

Output directory set to: `~/Boltz_Gen_RBX1`

Job created in file browser: **Boltz_Gen_RBX1** (2026-03-25, 9:22:44 PM)

**Step 6 — Advanced Settings (9.23.09 + 9.23.17)**
Full advanced parameter set exposed (all left at defaults for first run):

| Parameter | Default / Value | Notes |
|---|---|---|
| Alpha (Diversity Trade-off) | protocol-dependent | 0.0–1.0 range |
| Remove Composition Outliers | ON (checked) | Filters biased sequences |
| Skip Inverse Folding | OFF | Keep ProteinMPNN step |
| Skip Merge and Filter | OFF | Keep ranking step |
| Inverse Fold Num Sequences | 1 | Sequences per backbone |
| Inverse Fold Avoid | — | Motifs to avoid |
| Metrics Override | — | Custom scoring override |
| Additional Filters | — | Custom filter config |
| Size Buckets | — | Length bucketing |
| Refolding RMSD Threshold | — | RBX1 deformation filter (relevant!) |
| Config | — | Additional config |
| Step Scale | — | Diffusion step scale |
| Noise Scale | — | Diffusion noise scale |

**Refolding RMSD Threshold** is a notable parameter — directly relevant to the induced-fit failure mode we identified in our primary campaign (CUL1-WHB binders deformed RBX1 by 5 Å). Will set this explicitly on the next run.

**Step 7 — Protein input (9.29.28 + 9.29.47)**
Clicked "Add PROTEIN" → opened **Edit Protein** modal:
- Chain: **A**
- Sequence entered: `GGGGTNSGAGKKRFEVKKSNASAQSAWDIVVDNCAICRNHIMDLCIECQANQASATSEECTVAWGVCNHAFHFHCISRWLK`
  - This is the RBX1 RING domain (77 AA, same as primary campaign input)
- Fuse to Chain ID: empty
- Cyclic: unchecked
- Secondary Structure: none defined
- Binding Types: none defined

**Notable gap:** No explicit hotspot residue input visible in this modal. In our primary RFdiffusion pipeline, hotspot residues (`A4, A6, A12, A14, A15, A23, A24, A26, A28, A52, A56, A60, A64, A66`) were the key conditioning signal that directed binder generation to the E2-binding surface. It's unclear whether BoltGen infers binding sites automatically or whether hotspots can be specified via the **Binding Types** field. This needs investigation — if BoltGen has no hotspot conditioning, designs may target arbitrary surfaces of RBX1 rather than the therapeutically relevant RING-H2 face.

**Step 8 — Job submitted + running (9.30.29)**
Navigated to Workflow Status. Job confirmed running:

| Field | Value |
|---|---|
| Job ID | `boltzgen-652e8823` |
| Submitted | 2026-03-25, 9:30:18 PM |
| Status | **Running** |
| GPU | NVIDIA A100 |
| Output dir | ~/Boltz_Gen_RBX1 |

File browser shows `2LGV.pdb` (2415.1 KB) already saved in the output folder — the target structure pulled from PDB. Workflow Status panel is set to **auto-refresh**.

---

#### Full BoltGen parameter summary (Run 1)

| Parameter | Value |
|---|---|
| Protocol | Protein-Anything |
| Target | RBX1 RING domain, chain A, 77 AA |
| Number of Designs | 1 |
| Final Design Budget | 30 |
| GPU | NVIDIA Tesla A100 |
| Alpha | default |
| Hotspots | NOT SET (to investigate) |
| Refolding RMSD Threshold | NOT SET (to set next run) |
| Output | ~/Boltz_Gen_RBX1 |

---

#### BoltGen YAML — expected structure

BoltGen auto-generates a YAML from the UI inputs before submission. Based on Boltz-2's input format, the submitted YAML for Run 1 most likely contains:

```yaml
version: 1

sequences:
  - protein:
      id: A
      sequence: GGGGTNSGAGKKRFEVKKSNASAQSAWDIVVDNCAICRNHIMDLCIECQANQASATSEECTVAWGVCNHAFHFHCISRWLK
      msa: null

  - protein:
      id: B
      sequence: 70-90        # length range — BoltGen fills via diffusion

protocol: protein-anything
num_designs: 1
final_design_budget: 30
gpu: nvidia-a100
remove_composition_outliers: true
inverse_fold_num_sequences: 1
```

**What is almost certainly MISSING from Run 1 YAML:**

```yaml
# This block was NOT specified — hotspot conditioning absent
constraints:
  - pocket:
      binder: B
      contacts:
        - [A, 4]
        - [A, 6]
        - [A, 12]
        - [A, 14]
        - [A, 15]
        - [A, 23]
        - [A, 24]
        - [A, 26]
        - [A, 28]
        - [A, 52]
        - [A, 56]
        - [A, 60]
        - [A, 64]
        - [A, 66]
```

Without the hotspot constraints block, BoltGen will hallucinate a binder that docks *somewhere* on RBX1 — but not necessarily the therapeutically relevant RING-H2 E2-binding face. Run 1 is still valuable as a baseline: we can inspect which surface BoltGen chose autonomously and compare it to our target hotspot patch.

**Run 2 plan:** Add the full hotspot contacts block above, set `refolding_rmsd_threshold: 2.0`, and increase `num_designs`. See Entry 003.

---

#### Open questions going into results
1. Which surface did BoltGen target without hotspot conditioning?
2. Where does hotspot specification go in the UI — Binding Types field? String format? Direct YAML edit?
3. What metrics are returned per design — ipTM, ipLDDT, PAE matrices?
4. Are structure files (CIF/PDB) downloadable for each design?
5. What does "Final Design Budget: 30" mean exactly — 30 sequences total after filtering, or 30 per backbone?
6. What is the A100 wall time for this job?

---

#### Feedback for Brandon (Entry 002 additions)

| # | Category | Observation | Severity | Suggested fix |
|---|---|---|---|---|
| 6 | Feature gap | No visible hotspot/epitope input in the Edit Protein modal. Hotspot conditioning is critical for targeting specific binding surfaces (e.g., RBX1 RING-H2 E2-binding face). Without it, BoltGen may design binders to arbitrary surfaces. | High | Add hotspot residue field to Edit Protein dialog (residue IDs or selection string). Consider a 3D picker if the viewer is active. |
| 7 | UX | "Binding Types" field exists but has no documentation tooltip. Unclear if this is where hotspots/epitopes go. | Medium | Add inline help text or tooltip explaining Binding Types and providing an example |
| 8 | UX | Job IDs are auto-generated hashes (e.g., `boltzgen-652e8823`). No user-defined job name field visible. | Low | Allow optional user-defined job name alongside the hash ID for easier tracking across multiple runs |
| 9 | Feature | "Refolding RMSD Threshold" exists in Advanced Settings but has no default value shown or documentation. This is highly relevant for binder design — it should have a sensible default and tooltip explaining what it controls. | Medium | Add default value (e.g., 2.0 Å) and tooltip: "Maximum allowed RMSD of the target structure between input and predicted complex. Filters designs that require the target to deform to bind." |
| 10 | Feature | "Previous Execution Outputs" section allows accumulating designs across runs — excellent feature for iterative campaigns. Worth highlighting in onboarding. | Positive | — |

---

### Entry 003 — Job Cancellation + Run 2 Setup Investigation
**Date:** 2026-03-26
**Status:** COMPLETE
**Screenshots:** `Screenshot 2026-03-26 at 12.08.04 AM.png` → `Screenshot 2026-03-26 at 12.15.46 AM.png`

---

#### Step-by-step session log

**Step 1 — Job cancelled (12.08.04)**
Stopped `boltzgen-652e8823` via Workflow Status → stop button. Platform warned:
> "Warning: Any intermediate data will be lost if you stop these workflows."

No data lost — job had never actually started computing (no GPU usage detected in 2.5 hours).

**Step 2 — Returned to BoltGen form (12.15.13)**
Protein A still populated with RBX1 sequence. Noticed the **Constraints** section at the bottom of the Design Specification panel. Expanded it.

**Step 3 — Constraints section is insufficient for hotspots (12.15.17)**
Constraints panel only offers:
- **Add Bond Constraint**
- **Add Total Length Constraint**

Neither of these is hotspot/binding site conditioning. The YAML Preview panel below shows the auto-generated input:

```yaml
entities:
  - protein:
      id: A
      sequence: GGGGTNSGAGKKRFEVKKSNASAQSAWDIVVDNCAICRNHIMDLCIEQANQS...
```

**Critical finding: the YAML Preview is read-only.** It cannot be edited in the browser. Hotspot contacts cannot be injected manually via the UI. This is the core limitation blocking hotspot-conditioned design.

**Step 4 — GPU switched to T4 (12.15.24)**
Changed Hardware Configuration from A100 → **NVIDIA Tesla T4** ("Cost-effective, good for most workloads"). This should resolve the provisioning stall from Run 1.

**Step 5 — Refolding RMSD Threshold still not set (12.15.34)**
Field is visible in Advanced Settings but empty. Will set for Run 2 once hotspot issue is resolved. Output directory still pointing to old job path (`~/Boltz_Gen_RBX1/boltzgen-652e8823`) — needs updating.

**Step 6 — Discovered "Add Structure File" modal (`Screenshot 2026-03-26 at 12.15.42 AM.png` + `Screenshot 2026-03-26 at 12.15.46 AM.png`)**
Clicked **Add Structure File** button. This opens a richer **Edit Structure File** modal with fields that sequence-only input doesn't have:

| Field | Purpose |
|---|---|
| Structure File (CIF/PDB) | Upload the target structure |
| Include Chains | Specify which chains to use |
| Exclude Chains | Chains to ignore |
| **Design Regions** | Residues where the binder should be designed — **likely hotspot entry point** |
| **Fixed Positions (Not Designed)** | Target residues locked during design |
| **Binding Types** | Binding site specification |
| Structure Groups | Chain groupings |
| Secondary Structure | SS constraints |
| Design Insertions | Insert design regions |
| **Include Proximity** | Proximity-based site definition — **possible alternative hotspot method** |
| Fuse | Chain fusion |
| Reset Residue Index | Renumber residues |

**This is the correct path for hotspot-conditioned design.** The "Design Regions" and "Include Proximity" fields in the Structure File modal — not the sequence-based "Add PROTEIN" input — are almost certainly where hotspot residues need to be specified. The workflow requires uploading the RBX1 PDB file rather than just pasting the sequence.

---

#### Revised workflow for Run 2

Instead of:
> Add PROTEIN → paste RBX1 sequence → (no hotspot possible)

Use:
> **Add Structure File** → upload `rbx1_ring_renumbered.pdb` → expand **Design Regions** → specify residues `4,6,12,14,15,23,24,26,28,52,56,60,64,66`

Need to confirm:
- Does Design Regions accept a residue list? What format? (comma-separated? range notation?)
- Is "Include Proximity" an alternative — define a centre residue and radius rather than explicit list?
- Can the uploaded PDB be the full RING domain or does it need the full RBX1 sequence?

---

#### Run 2 parameters (revised plan)

| Parameter | Value |
|---|---|
| Input method | Add Structure File (not Add PROTEIN) |
| Structure file | `rbx1_ring_renumbered.pdb` (chain A, 77 AA) |
| Design Regions | `4,6,12,14,15,23,24,26,28,52,56,60,64,66` (RING-H2 hotspots) |
| Protocol | Protein-Anything |
| Number of Designs | 5 (increase from 1) |
| Final Design Budget | 30 |
| GPU | NVIDIA Tesla T4 |
| Refolding RMSD Threshold | 2.0 Å |
| Output directory | `~/Boltz_Gen_RBX1_hotspot` |

---

#### Feedback for Brandon (Entry 003 additions)

| # | Category | Observation | Severity | Suggested fix |
|---|---|---|---|---|
| 11 | Infrastructure | A100 jobs stall on provisioning indefinitely in beta — no compute detected after 2.5 hrs | High | Show estimated queue time; default to T4 for beta users; alert after 15 min wait |
| 12 | UX | 5 "design-pipeline" sub-jobs queued simultaneously with no explanation of job graph | Medium | Show expected step count before submission |
| 13 | Feature gap | YAML Preview is read-only — power users cannot inject custom parameters (e.g. hotspot contacts) that aren't exposed in the UI | High | Either make YAML editable, or add a "Raw YAML input" toggle for advanced users |
| 14 | UX / discoverability | Hotspot specification is only possible via "Add Structure File" → Design Regions, not via "Add PROTEIN" sequence input. This is non-obvious — most binder design users will start with a sequence paste and never find the hotspot field | High | Add a "Hotspot Residues" field directly to the Add PROTEIN modal, or show a tooltip on Add PROTEIN saying "For hotspot-conditioned design, use Add Structure File instead" |
| 15 | UX | "Design Regions", "Include Proximity", "Fixed Positions" in the Structure File modal have no tooltips explaining format or purpose | Medium | Add inline help text with example input format (e.g., `4,6,12-15` or `A:4,A:6`) |
| 16 | Bug | Output directory after cancellation still defaults to old job path (`boltzgen-652e8823`) in the rerun form | Low | Reset output directory to a fresh default on rerun |

---

### Entry 004 — Run 2 Submission: Hotspot-Conditioned Design on T4
**Date:** 2026-03-26
**Status:** RUNNING — awaiting results
**Screenshots:** `Screenshot 2026-03-26 at 12.21.46 AM.png`, `Screenshot 2026-03-26 at 12.23.07 AM.png`, `Screenshot 2026-03-26 at 12.28.08 AM.png`, `Screenshot 2026-03-26 at 12.28.27 AM.png`, `Screenshot 2026-03-26 at 12.30.26 AM.png`

---

#### Step-by-step session log

**Step 1 — Verified RBX1 PDB file (`Screenshot 2026-03-26 at 12.21.46 AM.png`)**
Opened `rbx1_ring_renumbered.pdb` directly in Chrome to verify structure before uploading. ATOM records confirmed clean: chain A, residues numbered from 1, 77 residues total. File already present in the ProteWorks sandbox file browser from the earlier upload attempt.

**Step 2 — Hotspot entry via Design Regions (`Screenshot 2026-03-26 at 12.23.07 AM.png`)**
Used **Add Structure File** → uploaded `rbx1_ring_renumbered.pdb` → expanded **Design Regions (1)**.

Format confirmed: dropdown for **Chain ID** + free text field for **Residue Index** (comma-separated integers).

Entered: `4,6,12,14,15,23,24,26,28,52,56,60,64,66`

This is the full RING-H2 E2-binding surface hotspot set — identical to the conditioning used in the primary RFdiffusion campaign.

**Step 3 — YAML Preview confirmed hotspot block (`Screenshot 2026-03-26 at 12.28.08 AM.png`)**
YAML Preview updated to show hotspot conditioning correctly encoded:

```yaml
entities:
  - file:
      path: "/tmp/inputs/rbx1_ring_renumbered.pdb"
      include: "all"
      design:
        - chain:
            id: A
            res_index: 4,6,12,14,15,23,24,26,28,52,56,60,64,66
constraints:
  - total_len:
      min: 70
      max: 90
```

Key observations:
- Hotspots encoded under `design.chain.res_index` within the structure file entity — not a separate `constraints.pocket` block as initially predicted
- Length constraint `70–90 AA` auto-populated from the Total Length Constraint field ✅
- YAML still read-only but now correctly reflects hotspot input via the UI

**Step 4 — Settings confirmed (`Screenshot 2026-03-26 at 12.28.27 AM.png`)**

| Parameter | Value |
|---|---|
| Protocol | Protein-Anything |
| Number of Designs | **10** (increased from 1) |
| Final Design Budget | 30 |
| GPU | NVIDIA Tesla T4 |

Note: Refolding RMSD Threshold was not set (field left blank). Will add in Run 3.

**Step 5 — Job submitted (`Screenshot 2026-03-26 at 12.30.26 AM.png`)**

| Field | Value |
|---|---|
| Job ID | `boltzgen-628e1038` |
| Submitted | 2026-03-26, 12:30:05 AM |
| Status | **Running** |
| GPU | NVIDIA Tesla T4 |
| Sub-jobs | 1 design-pipeline worker (vs 5 on A100 Run 1) |
| Old job | `boltzgen-652e8823` confirmed **Cancelled** |

T4 provisioned with only 1 sub-job — much cleaner than Run 1's 5 parallel workers that all stalled.

---

#### What to check when results arrive
1. Which surface did BoltGen target? Compare to hotspot patch — did `res_index` conditioning work?
2. How many designs returned vs budget of 30?
3. Metrics per design: ipTM, ipLDDT, PAE available?
4. Are CIF/PDB structure files downloadable per design?
5. Wall time on T4 for 10 designs / budget 30
6. Compare ipTM distribution vs primary campaign (GLMN mean 0.878, RFD mean ~0.82)

---

### Entry 004 — Addendum: Workflow Inputs JSON + Provisioning Failure
**Date:** 2026-03-26
**Screenshots:** `Screenshot 2026-03-26 at 1.00.32 AM.png`, `Screenshot 2026-03-26 at 3.45.11 AM.png`

#### Finding 1 — Actual backend JSON reveals hotspot conditioning was NOT passed

Clicked **Inputs** on the Workflow Steps bar to reveal the raw JSON sent to the BoltGen backend:

```json
{
  "workflow_name": "boltzgen",
  "execution_id": "boltzgen-628e1038",
  "submitted_at": "2026-03-26T04:30:05.400001+00:00",
  "files": {
    "extra_files": {
      "Boltz_Gen_RBX1/boltzgen-652e8823/rbx1_ring_renumbered.pdb"
    }
  },
  "parameters": {
    "protocol": "protein-anything",
    "num-designs": "10",
    "gpu-type": "nvidia-tesla-t4",
    "budget": "30",
    "refolding-rmsd-threshold": "2",
    "output_folder": "Boltz_Gen_RBX1/Boltz_Gen_RBX1_hotspot"
  }
}
```

**What made it through correctly:**
- `refolding-rmsd-threshold: "2"` ✅
- `gpu-type: "nvidia-tesla-t4"` ✅
- `num-designs: "10"` ✅
- `budget: "30"` ✅

**Critical bug — hotspot res_index missing from backend JSON:**
The `design.chain.res_index: 4,6,12,14,15,23,24,26,28,52,56,60,64,66` visible in the YAML Preview was **never serialised into the actual workflow parameters**. Despite the UI accepting the Design Regions input and the YAML Preview showing it, the backend received no hotspot conditioning. This run would produce unconditioned designs — effectively the same as Run 1.

**Secondary issue — stale file path:**
The PDB is referenced as `Boltz_Gen_RBX1/boltzgen-652e8823/rbx1_ring_renumbered.pdb` — the old cancelled job's directory. This is fragile; if that folder is cleaned up the file reference breaks.

#### Finding 2 — T4 also fails to provision after 3+ hours

At 3:45 AM (3h 15m after submission), job `boltzgen-628e1038` still shows:
- Overall Status: Running
- design-pipeline: "Waiting for resources to become available"
- Cost Information: "No compute usage has been detected"

T4 provisioning failure is as severe as the A100 failure in Run 1. This is a platform-wide infrastructure issue on the beta, not a configuration error. Cancelling this job.

#### Feedback for Brandon (Entry 004 additions)

| # | Category | Observation | Severity | Suggested fix |
|---|---|---|---|---|
| 17 | Bug (Critical) | Design Regions `res_index` field appears in YAML Preview but is silently dropped from the actual backend JSON. Hotspot conditioning is completely non-functional despite UI accepting input. | Critical | Validate that all UI fields serialise correctly to backend parameters; add a "Submitted Parameters" confirmation step showing exactly what the backend received |
| 18 | Bug | Uploaded structure file is referenced by old cancelled job's directory path (`boltzgen-652e8823/rbx1_ring_renumbered.pdb`). Fragile — breaks if old job folder is cleaned up. | High | File references should use a stable sandbox path, not a job-specific directory |
| 19 | Infrastructure | T4 GPU also fails to provision after 3+ hours with no error, no ETA, no retry logic. Platform appears resource-constrained on beta infrastructure overnight. | High | Implement a provisioning timeout (e.g. 30 min) with automatic retry or user notification; show queue position if resources are unavailable |
| 20 | UX | "Waiting for resources to become available" message gives users no actionable path — no cancel suggestion, no ETA, no alternative GPU option | Medium | Add: estimated wait time, link to cancel, option to switch GPU type without resubmitting full job |

---

## Results

*(Populated as designs come in)*

| ID | Length | ipTM | ipLDDT | ipSAE (Å) | composite | RING RMSD | Notes |
|---|---|---|---|---|---|---|---|
| | | | | | | | |

---

## References

- Adaptyv competition: https://design.adaptyv.bio
- Primary campaign notebook: `../lab_notebook.md`
- Nipah retrospective analysis: `../nipah_analysis/`
- Composite scoring rationale: `../submission_writeup.md`
- Boltz-2: Wohlwend et al. 2024, bioRxiv
- RFdiffusion (primary campaign): Watson et al. 2023, Nature
- ProteinMPNN (primary campaign): Dauparas et al. 2022, Science
