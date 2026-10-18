#!/usr/bin/env python3
"""Plant-type telomere units in Illumina short reads (no homopolymer errors).

Counts reads that hold an array of >= --copies consecutive units of the form
T{1,9}AGGG (G-strand) or CCCTA{1,9} (C-strand), and records the T-length of
every unit in those arrays. Because Illumina reads rarely mis-call homopolymer
length, the modal T-length shows the telomere unit without the T5/T6/T7 spread
seen in long reads. Telomeres are the most abundant such arrays in a genome,
so the dominant unit in a whole-genome read set is the telomere unit.

Usage: shortread_telomere_units.py LABEL reads.fasta [--copies 4]
"""
import argparse
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_telomere_ends import read_fasta  # noqa: E402

G = re.compile(r"(?:T{1,9}AGGG){4,}")
C = re.compile(r"(?:CCCTA{1,9}){4,}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("reads")
    a = ap.parse_args()
    n = 0
    reads_by_mode = collections.Counter()
    pure_by_mode = collections.Counter()
    units = collections.Counter()
    for _, s in read_fasta(a.reads):
        n += 1
        s = s.upper()
        lens = [len(t) for m in G.finditer(s) for t in re.findall(r"(T+)AGGG", m.group())]
        lens += [len(x) for m in C.finditer(s) for x in re.findall(r"CCCT(A+)", m.group())]
        if not lens:
            continue
        units.update(lens)
        mode = collections.Counter(lens).most_common(1)[0][0]
        reads_by_mode[mode] += 1
        if lens.count(mode) >= 0.9 * len(lens):
            pure_by_mode[mode] += 1
    tot_reads = sum(reads_by_mode.values())
    tot_units = sum(units.values()) or 1
    print(f"# {a.label}: {n:,} reads; {tot_reads:,} with an array of >=4 T(1-9)AGGG units "
          f"({1e6 * tot_reads / max(n, 1):.0f} per million)")
    print("T-length\treads_with_this_mode\tper_million_reads\tof_which_>=90%_pure\t%_of_all_units")
    for k in range(1, 10):
        if reads_by_mode[k] or units[k]:
            print(f"T{k}\t{reads_by_mode[k]}\t{1e6 * reads_by_mode[k] / n:.0f}\t{pure_by_mode[k]}\t"
                  f"{100 * units[k] / tot_units:.1f}")


if __name__ == "__main__":
    main()
