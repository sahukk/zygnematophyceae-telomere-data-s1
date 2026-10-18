#!/usr/bin/env python3
"""Uniform, unbiased telomere-repeat verification by end-vs-interior enrichment.

For every sequence >= --min-len in an assembly:
  * the terminal --win bp at each end is searched for tandem repeats with a
    primitive unit of 4-12 bp (no motif is assumed in advance);
  * the same search is run on interior windows of the same size (default 3 per
    sequence, evenly spaced) as a background control.

A unit "occupies" a window if tandem copies of it cover >= --min-cov of the
window. Units are collapsed to one canonical form (minimal rotation over both
strands), so TTTAGGG, AGGGTTT and CCCTAAA are counted as the same repeat.

For each assembly the report gives the top units by number of ends occupied,
the same unit's occupancy rate in interior windows, and an enrichment ratio.
A telomere call needs BOTH high end occupancy and strong end enrichment.

Usage: verify_telomere_ends.py LABEL ASSEMBLY.fasta[.gz] [--win 1000] ...
Output: one TSV block on stdout (LABEL, unit, ends_occupied, n_ends,
        end_frac, interior_frac, enrichment).
"""
import argparse
import gzip
from collections import Counter

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def primitive(u):
    n = len(u)
    for p in range(1, n):
        if n % p == 0 and u[:p] * (n // p) == u:
            return u[:p]
    return u


def canonical(u):
    rots = [u[i:] + u[:i] for i in range(len(u))]
    rc = revcomp(u)
    rots += [rc[i:] + rc[:i] for i in range(len(rc))]
    return min(rots)


def tandem_units(w, kmin=4, kmax=12, min_cov=0.3):
    """Return the set of canonical units whose tandem copies cover >= min_cov of w."""
    L = len(w)
    found = set()
    for k in range(kmin, kmax + 1):
        covered = Counter()
        i = 0
        while i + 2 * k <= L:
            u = w[i:i + k]
            if "N" in u or w[i + k:i + 2 * k] != u:
                i += 1
                continue
            j = i + k
            while j + k <= L and w[j:j + k] == u:
                j += k
            pu = primitive(u)
            if len(pu) >= kmin:
                covered[canonical(pu)] += j - i
            i = j
        for cu, c in covered.items():
            if c / L >= min_cov and len(set(cu)) >= 2:
                found.add(cu)
    return found


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().upper())
    if name:
        yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("fasta")
    ap.add_argument("--win", type=int, default=1000)
    ap.add_argument("--min-len", type=int, default=100_000)
    ap.add_argument("--min-cov", type=float, default=0.3)
    ap.add_argument("--interior", type=int, default=3)
    ap.add_argument("--top", type=int, default=3)
    a = ap.parse_args()

    end_hits, int_hits = Counter(), Counter()
    n_ends = n_int = n_seq = 0
    for _, s in read_fasta(a.fasta):
        L = len(s)
        if L < a.min_len:
            continue
        n_seq += 1
        # strip terminal Ns so a padded end is still tested
        left = s.lstrip("N")[:a.win]
        right = s.rstrip("N")[-a.win:]
        for w in (left, right):
            n_ends += 1
            end_hits.update(tandem_units(w, min_cov=a.min_cov))
        for t in range(1, a.interior + 1):
            st = L * t // (a.interior + 1)
            n_int += 1
            int_hits.update(tandem_units(s[st:st + a.win], min_cov=a.min_cov))

    if n_ends == 0:
        print(f"{a.label}\tNA\t0\t0\t0\t0\tNA\t# no sequences >= {a.min_len} bp")
        return
    for u, c in end_hits.most_common(a.top):
        ef = c / n_ends
        inf = int_hits[u] / n_int if n_int else 0
        enr = ef / inf if inf else float("inf")
        print(f"{a.label}\t{u}\t{c}\t{n_ends}\t{ef:.3f}\t{inf:.4f}\t"
              f"{'inf' if enr == float('inf') else f'{enr:.1f}'}\t# {n_seq} seqs")
    if not end_hits:
        print(f"{a.label}\tnone\t0\t{n_ends}\t0\t0\tNA\t# {n_seq} seqs")


if __name__ == "__main__":
    main()
