#!/usr/bin/env python3
"""Candidate-motif density test, tolerant of degenerate telomeric arrays.

Complements verify_telomere_ends.py (which needs perfect consecutive copies).
For each candidate unit, every occurrence of the unit or any rotation of it, on
either strand, is counted in the terminal --win bp of each sequence and in
interior windows. A window is "positive" when those occurrences cover at least
--min-frac of it. Testing the earlier manuscript's motif calls alongside the
newly discovered units shows directly which call the genome supports.

Usage: motif_end_density.py LABEL ASSEMBLY.fasta UNIT [UNIT ...]
Output TSV: label, unit, ends_positive, n_ends, end_frac, interior_frac
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_telomere_ends import read_fasta, revcomp  # noqa: E402


def rotations(u):
    out = set()
    for s in (u, revcomp(u)):
        out.update(s[i:] + s[:i] for i in range(len(s)))
    return out


def covered(w, pattern, k):
    pos = set()
    for m in pattern.finditer(w):
        pos.update(range(m.start(), m.start() + k))
    return len(pos) / max(len(w), 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("fasta")
    ap.add_argument("units", nargs="+")
    ap.add_argument("--win", type=int, default=1000)
    ap.add_argument("--min-len", type=int, default=100_000)
    ap.add_argument("--min-frac", type=float, default=0.2)
    ap.add_argument("--interior", type=int, default=3)
    a = ap.parse_args()

    pats = {u: re.compile("(?=(" + "|".join(sorted(rotations(u))) + "))")
            for u in a.units}
    end_pos = dict.fromkeys(a.units, 0)
    int_pos = dict.fromkeys(a.units, 0)
    n_end = n_int = 0
    for _, s in read_fasta(a.fasta):
        L = len(s)
        if L < a.min_len:
            continue
        ends = [s.lstrip("N")[:a.win], s.rstrip("N")[-a.win:]]
        mids = [s[L * t // (a.interior + 1):][:a.win]
                for t in range(1, a.interior + 1)]
        n_end += 2
        n_int += len(mids)
        for u, p in pats.items():
            k = len(u)
            end_pos[u] += sum(covered(w, p, k) >= a.min_frac for w in ends)
            int_pos[u] += sum(covered(w, p, k) >= a.min_frac for w in mids)
    for u in a.units:
        ef = end_pos[u] / n_end if n_end else 0
        inf = int_pos[u] / n_int if n_int else 0
        print(f"{a.label}\t{u}\t{end_pos[u]}\t{n_end}\t{ef:.3f}\t{inf:.4f}")


if __name__ == "__main__":
    main()
