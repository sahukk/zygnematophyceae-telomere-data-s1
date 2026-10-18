#!/usr/bin/env python3
"""T-tract composition of plant-type telomere arrays at read ends.

For every long read (ONT adapter remnants trimmed as in read_end_motif_test.py),
look for arrays of >= --copies consecutive units of the general form T{1,9}AGGG
in the correct telomeric orientation: G-strand (T..TAGGG)n in the last --win bp
of the read, C-strand (CCCTA..A)n in the first --win bp. The same search in two
interior windows of the read gives the background.

Each array is split into its units and the number of T's per unit is counted.
This answers "which unit is the telomere made of, and how many units of other
lengths are mixed in?", e.g. whether any TTTAGGG (T3) units sit inside a
TTTTTTAGGG (T6) telomere.

Usage: telomere_unit_composition.py LABEL reads.fasta [--win 1000] [--copies 4]
Output: a summary on stdout; with --tsv FILE, one line per array.
"""
import argparse
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_telomere_ends import read_fasta  # noqa: E402
from read_end_motif_test import trim  # noqa: E402


def arrays(seq, strand, copies):
    """Yield lists of T-tract lengths, one list per array, G-strand sense."""
    if strand == "G":
        pat = re.compile(r"(?:T{1,9}AGGG){%d,}" % copies)
        for m in pat.finditer(seq):
            yield [len(t) for t in re.findall(r"(T+)AGGG", m.group())]
    else:
        pat = re.compile(r"(?:CCCTA{1,9}){%d,}" % copies)
        for m in pat.finditer(seq):
            yield [len(a) for a in re.findall(r"CCCT(A+)", m.group())]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("reads")
    ap.add_argument("--win", type=int, default=1000)
    ap.add_argument("--copies", type=int, default=4)
    ap.add_argument("--tsv")
    a = ap.parse_args()

    n_reads = 0
    end_units, int_units = collections.Counter(), collections.Counter()
    end_arrays = int_arrays = 0
    modal = collections.Counter()          # modal unit length of each end array
    t3_in_other = 0                        # end arrays of another mode holding >=1 T3 unit
    pure = collections.Counter()           # end arrays where >=90% of units share the mode
    out = open(a.tsv, "w") if a.tsv else None
    if out:
        out.write("read\twhere\tn_units\tmode\tunit_T_lengths\n")

    for name, s in read_fasta(a.reads):
        s, _, _ = trim(s.upper())
        if len(s) < 3 * a.win:
            continue
        n_reads += 1
        windows = [("start", s[:a.win], "C"), ("end", s[-a.win:], "G")]
        mid = len(s) // 2
        interior = [("int_a", s[mid - a.win:mid]), ("int_b", s[mid:mid + a.win])]
        for where, w, strand in windows:
            for units in arrays(w, strand, a.copies):
                end_arrays += 1
                end_units.update(units)
                mode = collections.Counter(units).most_common(1)[0][0]
                modal[mode] += 1
                if units.count(mode) >= 0.9 * len(units):
                    pure[mode] += 1
                if mode != 3 and 3 in units:
                    t3_in_other += 1
                if out:
                    out.write(f"{name}\t{where}\t{len(units)}\t{mode}\t{','.join(map(str, units))}\n")
        for where, w in interior:
            for strand in "GC":
                for units in arrays(w, strand, a.copies):
                    int_arrays += 1
                    int_units.update(units)

    print(f"# {a.label}: {n_reads} reads >= {3 * a.win} bp; arrays = >= {a.copies} units "
          f"T(1-9)AGGG, correct orientation in terminal {a.win} bp")
    print(f"end arrays: {end_arrays}   interior arrays (both strands, 2 windows/read): {int_arrays}")
    print("T-length\tarrays_with_this_mode\tof_which_>=90%_pure\tunits_at_ends\t%_of_end_units\tunits_in_interiors")
    tot = sum(end_units.values()) or 1
    for k in range(1, 10):
        if end_units[k] or int_units[k] or modal[k]:
            print(f"T{k}\t{modal[k]}\t{pure[k]}\t{end_units[k]}\t{100 * end_units[k] / tot:.1f}\t{int_units[k]}")
    print(f"end arrays whose mode is not T3 but contain >=1 T3 unit: {t3_in_other}")


if __name__ == "__main__":
    main()
