#!/usr/bin/env python3
"""Targeted read-end test for candidate plant-type telomere units.

For each long read: strip ONT adapter remnants from both ends (core
GTTACGTATTGCT / AGCAATACGTAAC within the terminal 250 bp), then ask whether
the terminal --win bp holds an array of >= --copies consecutive copies of each
candidate unit (any rotation, either strand). Interior windows of the same
size (read midpoint) give the background. Units are tested longest-first and a
window is assigned to the first unit that matches, so TTTTTAGGG arrays are not
double-counted as TTTTAGGG/TTTAGGG.

Usage: read_end_motif_test.py LABEL reads.fasta [--win 500] [--copies 4]
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_telomere_ends import read_fasta, revcomp  # noqa: E402

UNITS = ["TTTTTTAGGG", "TTTTTAGGG", "TTTTAGGG", "TTTAGGG", "TTAGGG"]
ADAPT = re.compile("GTTACGTATTGCT|AGCAATACGTAAC")


def rotations(u):
    return {u[i:] + u[:i] for i in range(len(u))}


def patterns(copies):
    pats = []
    for u in UNITS:
        alts = set()
        for r in rotations(u):
            alts.add(r * copies)
            alts.add(revcomp(r) * copies)
        pats.append((u, re.compile("|".join(sorted(alts)))))
    return pats


def trim(s):
    head = [m.end() for m in ADAPT.finditer(s[:250])]
    tail = [m.start() for m in ADAPT.finditer(s[-250:])]
    lo = max(head) if head else 0
    hi = len(s) - 250 + min(tail) if tail else len(s)
    return s[lo:hi], bool(head), bool(tail)


def classify(w, pats):
    for u, p in pats:
        if p.search(w):
            return u
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("reads")
    ap.add_argument("--win", type=int, default=500)
    ap.add_argument("--copies", type=int, default=4)
    a = ap.parse_args()
    pats = patterns(a.copies)
    ends = {u: 0 for u in UNITS}; ints = {u: 0 for u in UNITS}
    adapt_ends = {u: 0 for u in UNITS}
    n_reads = n_ends = n_adapt = 0
    for _, s in read_fasta(a.reads):
        s, h, t = trim(s.upper())
        if len(s) < 3 * a.win:
            continue
        n_reads += 1
        for w, marked in ((s[:a.win], h), (s[-a.win:], t)):
            n_ends += 1; n_adapt += marked
            u = classify(w, pats)
            if u:
                ends[u] += 1; adapt_ends[u] += marked
        mid = len(s) // 2
        for w in (s[mid - a.win:mid], s[mid:mid + a.win]):
            u = classify(w, pats)
            if u:
                ints[u] += 1
    print(f"# {a.label}: {n_reads} reads; {n_ends} ends ({n_adapt} adapter-marked native ends); "
          f"array = >={a.copies} consecutive copies in {a.win} bp")
    print("unit\tends\tper_10k_ends\tinteriors\tper_10k_interiors\tenrichment\tadapter_marked_ends")
    for u in UNITS:
        e, i = ends[u], ints[u]
        enr = f"{(e / n_ends) / (i / n_ends):.1f}" if i else ("inf" if e else "-")
        print(f"{u}\t{e}\t{1e4 * e / n_ends:.1f}\t{i}\t{1e4 * i / n_ends:.1f}\t{enr}\t{adapt_ends[u]}")


if __name__ == "__main__":
    main()
