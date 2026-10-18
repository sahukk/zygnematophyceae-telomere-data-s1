#!/usr/bin/env python3
"""Read-level telomere repeat discovery (assembly-independent).

Long reads (ONT/PacBio) that originate at a chromosome end begin or finish in
the telomeric array, so the true telomere unit should be strongly enriched in
the terminal --win bp of reads relative to read interiors. Assembly trimming
cannot hide it.

Short reads (--mode short): reads that consist almost entirely (>= --short-cov)
of one tandem unit are counted; telomeric arrays are usually the most abundant
such units after rDNA/satellites, so this gives an independent abundance check.

Uses the same unbiased tandem-unit detector as verify_telomere_ends.py.

Usage:
  read_end_telomere.py LABEL reads.fasta --mode long  [--max-reads 150000]
  read_end_telomere.py LABEL reads.fasta --mode short [--max-reads 2000000]
"""
import argparse
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_telomere_ends import read_fasta, tandem_units  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("reads")
    ap.add_argument("--mode", choices=["long", "short"], required=True)
    ap.add_argument("--win", type=int, default=500)
    ap.add_argument("--min-read", type=int, default=5000)
    ap.add_argument("--max-reads", type=int, default=150_000)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--short-cov", type=float, default=0.8)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    random.seed(a.seed)

    end_hits, int_hits = Counter(), Counter()
    n = 0
    for _, s in read_fasta(a.reads):
        if a.mode == "long":
            if len(s) < a.min_read:
                continue
            n += 1
            for w in (s[:a.win], s[-a.win:]):
                end_hits.update(tandem_units(w, min_cov=a.min_cov))
            mid = len(s) // 2
            int_hits.update(tandem_units(s[mid:mid + a.win], min_cov=a.min_cov))
        else:
            n += 1
            end_hits.update(tandem_units(s, min_cov=a.short_cov))
        if n >= a.max_reads:
            break

    print(f"# {a.label}: {n} reads examined (mode={a.mode})")
    if a.mode == "long":
        print("label\tunit\tread_ends_with_array\tper_10k_ends\tper_10k_interiors\tenrichment")
        for u, c in end_hits.most_common(a.top):
            e = c / (2 * n) * 1e4
            i = int_hits[u] / n * 1e4
            enr = f"{e / i:.1f}" if i else "inf"
            print(f"{a.label}\t{u}\t{c}\t{e:.1f}\t{i:.1f}\t{enr}")
    else:
        print("label\tunit\treads_mostly_this_unit\tper_million_reads")
        for u, c in end_hits.most_common(a.top):
            print(f"{a.label}\t{u}\t{c}\t{c / n * 1e6:.1f}")


if __name__ == "__main__":
    main()
