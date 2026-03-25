#!/usr/bin/env python3
"""
Novelty screen via EBI BLAST REST API against SwissProt.
Flags sequences with >75% identity to any known protein.
For competition: a ProteinMPNN redesign of GLMN/CUL1-WHB is expected to
show high identity to its own scaffold (natural protein). We use this to:
  1. Flag any sequence that is essentially unmodified from its scaffold
  2. Report the scaffold identity so we can decide which sequences are
     genuinely redesigned vs. just copied
"""
import requests, time, json, csv, os, re
from pathlib import Path

EMAIL = "research@example.com"   # required by EBI API
BASE  = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
FASTA = Path("candidates_novelty_screen.fasta")
OUT_CSV = Path("novelty_screen_results.csv")

# ---- parse FASTA --------------------------------------------------------
def parse_fasta(fpath):
    seqs = {}
    cur_id = None
    for line in fpath.read_text().splitlines():
        if line.startswith(">"):
            cur_id = line[1:].strip()
            seqs[cur_id] = ""
        elif cur_id:
            seqs[cur_id] += line.strip()
    return seqs

# ---- EBI BLAST helpers --------------------------------------------------
def submit(seq_id, sequence):
    params = {
        "email": EMAIL,
        "program": "blastp",
        "database": "uniprotkb_swissprot",   # SwissProt first (fast)
        "sequence": sequence,
        "stype": "protein",
        "exp": "1e-3",
        "scores": "5",
        "alignments": "5",
        "format": "json",
    }
    r = requests.post(f"{BASE}/run", data=params, timeout=30)
    r.raise_for_status()
    job_id = r.text.strip()
    print(f"  Submitted {seq_id}: job {job_id}")
    return job_id

def wait_for_result(job_id, max_wait=300):
    for _ in range(max_wait // 10):
        time.sleep(10)
        r = requests.get(f"{BASE}/status/{job_id}", timeout=15)
        status = r.text.strip()
        if status in ("FINISHED", "ERROR", "FAILURE", "NOT_FOUND"):
            return status
    return "TIMEOUT"

def get_results(job_id):
    r = requests.get(f"{BASE}/result/{job_id}/json", timeout=30)
    if r.status_code != 200:
        return None
    return r.json()

def parse_best_hit(result_json):
    """Return (hit_id, pident, evalue) for the top hit, or None."""
    try:
        hits = result_json["hits"]
        if not hits:
            return None
        best = hits[0]
        hsp = best["hit_hsps"][0]
        pident = float(hsp["hsp_identity"]) / float(hsp["hsp_align_len"]) * 100.0
        return best["hit_acc"], round(pident, 1), float(hsp["hsp_expect"])
    except Exception:
        return None

# ---- main ---------------------------------------------------------------
def main():
    seqs = parse_fasta(FASTA)
    print(f"Loaded {len(seqs)} sequences from {FASTA}")

    # Resume if partial results exist
    done = {}
    if OUT_CSV.exists():
        with open(OUT_CSV) as f:
            for row in csv.DictReader(f):
                done[row["seq_id"]] = row
        print(f"Resuming: {len(done)} already done")

    results = list(done.values())

    remaining = {k: v for k, v in seqs.items() if k not in done}
    print(f"Submitting {len(remaining)} sequences in batches...")

    # Submit in batches of 10 (EBI rate limit ~10 concurrent)
    batch_size = 10
    items = list(remaining.items())

    for batch_start in range(0, len(items), batch_size):
        batch = items[batch_start:batch_start + batch_size]
        jobs = {}

        for seq_id, sequence in batch:
            try:
                jid = submit(seq_id, sequence)
                jobs[jid] = (seq_id, sequence)
            except Exception as e:
                print(f"  ERROR submitting {seq_id}: {e}")
                results.append({
                    "seq_id": seq_id, "seq_len": len(sequence),
                    "best_hit_id": "SUBMIT_ERROR", "best_pident": "NA",
                    "evalue": "NA", "novel": "unknown"
                })
            time.sleep(1)

        print(f"  Waiting for {len(jobs)} jobs...")
        for jid, (seq_id, sequence) in jobs.items():
            status = wait_for_result(jid)
            if status == "FINISHED":
                data = get_results(jid)
                hit = parse_best_hit(data) if data else None
                if hit:
                    hit_id, pident, evalue = hit
                    novel = "YES" if pident < 75.0 else "NO"
                    results.append({
                        "seq_id": seq_id, "seq_len": len(sequence),
                        "best_hit_id": hit_id, "best_pident": pident,
                        "evalue": evalue, "novel": novel
                    })
                    flag = "NOVEL" if novel == "YES" else f"  hit={pident:.1f}%"
                    print(f"    {seq_id}: {flag}")
                else:
                    results.append({
                        "seq_id": seq_id, "seq_len": len(sequence),
                        "best_hit_id": "NO_HIT", "best_pident": 0.0,
                        "evalue": 1.0, "novel": "YES"
                    })
                    print(f"    {seq_id}: NO_HIT → NOVEL")
            else:
                print(f"    {seq_id}: job {jid} status={status}")
                results.append({
                    "seq_id": seq_id, "seq_len": len(sequence),
                    "best_hit_id": f"JOB_{status}", "best_pident": "NA",
                    "evalue": "NA", "novel": "unknown"
                })

        # Save after every batch
        with open(OUT_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["seq_id","seq_len","best_hit_id","best_pident","evalue","novel"])
            writer.writeheader()
            writer.writerows(results)
        print(f"  Progress saved ({len(results)}/{len(seqs)} done)")

    print(f"\nDone! Results in {OUT_CSV}")
    novel_seqs = [r for r in results if r["novel"] == "YES"]
    not_novel  = [r for r in results if r["novel"] == "NO"]
    unknown    = [r for r in results if r["novel"] == "unknown"]
    print(f"  Novel (<75% id):   {len(novel_seqs)}")
    print(f"  Not novel (≥75%):  {len(not_novel)}")
    print(f"  Unknown/error:     {len(unknown)}")

if __name__ == "__main__":
    main()
