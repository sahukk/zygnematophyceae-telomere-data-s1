#!/usr/bin/env python3
"""Error-tolerant, motif-free telomere discovery from long reads.

read_end_telomere.py needs perfect consecutive repeat copies, which noisy
reads (PacBio CLR, older ONT) rarely contain. This test instead counts every
k-mer (default k=6, both strands collapsed) in the terminal --win bp of reads
and in an interior window of the same size, then ranks k-mers by end/interior
enrichment. A telomeric array shows up as a family of related k-mers (rotations
of the unit) with high end counts and strong enrichment, even when individual
copies carry sequencing errors.

Usage: kmer_end_enrichment.py LABEL reads.fasta [--k 6] [--max-reads 60000]
"""
import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_telomere_ends import read_fasta, revcomp  # noqa: E402


def count(w, k, c):
    for i in range(len(w) - k + 1):
        km = w[i:i + k]
        if "N" in km:
            continue
        rc = revcomp(km)
        c[km if km < rc else rc] += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("reads")
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--win", type=int, default=300)
    ap.add_argument("--min-read", type=int, default=5000)
    ap.add_argument("--max-reads", type=int, default=60000)
    ap.add_argument("--top", type=int, default=12)
    a = ap.parse_args()

    ends, mids = Counter(), Counter()
    n = 0
    for _, s in read_fasta(a.reads):
        if len(s) < a.min_read:
            continue
        n += 1
        count(s[:a.win], a.k, ends)
        count(s[-a.win:], a.k, ends)
        m = len(s) // 2
        count(s[m:m + a.win], a.k, mids)
        if n >= a.max_reads:
            break
    # interior sampled once per read vs. two ends: scale interior x2
    tot_e = sum(ends.values()) or 1
    tot_m = 2 * sum(mids.values()) or 1
    rows = []
    for km, e in ends.items():
        enr = (e / tot_e) / ((2 * mids[km] + 1) / tot_m)
        rows.append((enr, e, km))
    rows = [r for r in rows if r[1] >= 50 and len(set(r[2])) >= 2]
    rows.sort(reverse=True)
    print(f"# {a.label}: {n} reads, k={a.k}, end window {a.win} bp")
    print("label\tkmer\tend_count\tinterior_count_x2\tenrichment")
    for enr, e, km in rows[:a.top]:
        print(f"{a.label}\t{km}\t{e}\t{2 * mids[km]}\t{enr:.1f}")


if __name__ == "__main__":
    main()
