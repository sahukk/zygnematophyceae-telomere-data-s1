#!/usr/bin/env python3
"""Combine Phase 2 TERT evidence into one per-assembly table.

Evidence sources (all under data/processed/tert_phase2/):
  A. tert_phase2_results.tsv   - miniprot models, 25 external references (pass 1)
  B. tert_pass2_results.tsv    - miniprot models, + pass-1 survey TERTs of OTHER
                                 species (pass 2, only for pass-1 negatives)
  C. orf_scan/<label>.tbl      - reference-free: all ORFs >= 80 aa, TRBD at GA
     orf_scan/hits/<label>.domtbl (optional) - full domain content of those ORFs
  D. busco/busco_summary.txt   - compleasm eukaryota_odb12 completeness

Call levels:
  Complete  - one locus/ORF carries TRBD + RVT_1 + (TEN_TERT or TERT_C)
  Partial   - TRBD present, but not the full domain set (fragmented model/ORF)
  Not detected - no TRBD by any of A, B, C (interpreted against BUSCO score)
"""
import csv
import glob
import os
import re
import sys

P = sys.argv[1] if len(sys.argv) > 1 else "data/processed/tert_phase2"


def miniprot_calls(path, src):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("#") or line.startswith("label"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "TERT":
            continue
        doms = set(re.findall(r"(\w+)\(", f[3]))
        full = {"Telomerase_RBD", "RVT_1"} <= doms and bool(doms & {"TEN_TERT", "TERT_C"})
        lvl = "Complete" if full else "Partial"
        cur = out.get(f[0])
        if cur is None or (lvl == "Complete" and cur[0] != "Complete"):
            out[f[0]] = (lvl, f"{src}:{f[1]} [{'+'.join(sorted(doms))}]")
        out.setdefault(f[0] + "#n", 0)
        out[f[0] + "#n"] += 1
    return out


def orf_calls():
    out = {}
    for tbl in glob.glob(os.path.join(P, "orf_scan", "*.tbl")):
        lab = os.path.basename(tbl)[:-4]
        ids = [l.split()[0] for l in open(tbl) if not l.startswith("#")]
        if not ids:
            continue
        lvl, note = "Partial", f"ORF-scan: {len(ids)} TRBD ORF(s)"
        dt = os.path.join(P, "orf_scan", "hits", lab + ".domtbl")
        if os.path.exists(dt):
            per = {}
            for l in open(dt):
                if l.startswith("#"):
                    continue
                f = l.split()
                per.setdefault(f[3], set()).add(f[0])
            for orf, doms in per.items():
                if {"Telomerase_RBD", "RVT_1"} <= doms and doms & {"TEN_TERT", "TERT_C"}:
                    lvl = "Complete"
                    note += f"; {orf} [{'+'.join(sorted(doms))}]"
        out[lab] = (lvl, note)
    return out


def busco():
    out = {}
    path = os.path.join(P, "busco", "busco_summary.txt")
    if os.path.exists(path):
        for line in open(path):
            f = line.split()
            if not f:
                continue
            m = re.search(r"S:([\d.]+)%.*?D:([\d.]+)%", line)
            if m:
                out[f[0]] = round(float(m.group(1)) + float(m.group(2)), 1)
    return out


def main():
    labels = [l.split("\t")[0] for l in open(os.path.join(P, "manifest_by_size.tsv"))]
    A = miniprot_calls(os.path.join(P, "tert_phase2_results.tsv"), "pass1")
    B = miniprot_calls(os.path.join(P, "tert_pass2_results.tsv"), "pass2")
    C = orf_calls()
    D = busco()
    rank = {"Complete": 2, "Partial": 1}
    w = csv.writer(sys.stdout, delimiter="\t")
    w.writerow(["assembly", "TERT_call", "n_TERT_loci_models", "busco_complete_pct",
                "miniprot_pass1", "miniprot_pass2", "orf_scan"])
    for lab in sorted(labels):
        ev = [x for x in (A.get(lab), B.get(lab), C.get(lab)) if x]
        call = max((e[0] for e in ev), key=lambda v: rank[v]) if ev else "Not detected"
        n = A.get(lab + "#n", 0) or B.get(lab + "#n", 0)
        w.writerow([lab, call, n, D.get(lab, ""),
                    (A.get(lab) or ("", "-"))[1], (B.get(lab) or ("", "-"))[1],
                    (C.get(lab) or ("", "-"))[1]])


if __name__ == "__main__":
    main()
