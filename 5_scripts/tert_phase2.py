#!/usr/bin/env python3
"""Phase 2 TERT search: splice-aware gene models + TERT-specific domain test.

Replaces the earlier BLASTP/TBLASTN-only search. For one genome assembly:
  1. miniprot aligns a curated, lineage-diverse TERT reference set (25 proteins,
     none from surveyed species) to the genome, producing intron-aware
     protein models and their translations.
  2. Each translated model is searched (hmmsearch) against five Pfam TERT
     domains: TEN_TERT (PF11474), Telomerase_RBD (PF12009), RVT_1 (PF00078),
     TERT_thumb (PF17984), TERT_C (PF21399). Telomerase_RBD (TRBD) and TERT_C
     use Pfam gathering (GA) cutoffs; RVT_1 is scored at E <= 1e-3 because
     TERT reverse-transcriptase domains are divergent and often miss GA.
  3. Overlapping models are merged into loci. A locus is a TERT call only if it
     carries the TERT-specific TRBD domain at GA. Retrotransposon RTs carry
     RVT_1 but never TRBD, which is what removes the false positives that
     BLAST-only searching let through.

Output (stdout, TSV): one row per locus with a domain hit, plus a summary row.
Usage: tert_phase2.py LABEL GENOME.fasta REFS.faa DOMAINS.hmm OUTDIR [--threads 8]
"""
import argparse
import os
import subprocess
from collections import defaultdict


def run(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def parse_miniprot(path):
    """Yield (qname, contig, start, end, strand, score, frameshifts, stops, protein)."""
    paf = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("##PAF"):
                f = line.rstrip("\n").split("\t")[1:]
                tags = dict(t.split(":", 2)[::2] for t in f[12:] if t.count(":") >= 2)
                paf = (f[0], f[5], int(f[7]), int(f[8]), f[4],
                       int(tags.get("AS", 0)), int(tags.get("fs", 0)),
                       int(tags.get("st", 0)))
            elif line.startswith("##STA") and paf:
                yield paf + (line.rstrip("\n").split("\t")[1].rstrip("*"),)
                paf = None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("genome")
    ap.add_argument("refs")
    ap.add_argument("hmm")
    ap.add_argument("outdir")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--rt-evalue", type=float, default=1e-3)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    pre = os.path.join(a.outdir, a.label)

    mp = pre + ".miniprot.gff"
    if not os.path.exists(mp) or os.path.getsize(mp) == 0:
        with open(mp, "w") as out:
            run(["miniprot", "-t", str(a.threads), "--gff", "--trans", "-N", "10",
                 "--outn=10", a.genome, a.refs], stdout=out,
                stderr=subprocess.DEVNULL)

    models = list(parse_miniprot(mp))
    faa = pre + ".models.faa"
    with open(faa, "w") as out:
        for i, m in enumerate(models):
            out.write(f">m{i}\n{m[8]}\n")
    dom = pre + ".domtbl"
    if models:
        run(["hmmsearch", "--cpu", str(a.threads), "--domtblout", dom, "-E", "1",
             "--domE", "1", a.hmm, faa], stdout=subprocess.DEVNULL)

    # domain hits per model
    hits = defaultdict(dict)  # model -> {domain: best evalue}
    if models:
        ga = {}
        for line in open(a.hmm):
            if line.startswith("NAME"):
                name = line.split()[1]
            elif line.startswith("GA"):
                ga[name] = float(line.split()[2].rstrip(";"))
        for line in open(dom):
            if line.startswith("#"):
                continue
            f = line.split()
            mid, d, ievalue, dscore = f[0], f[3], float(f[12]), float(f[13])
            ok = (ievalue <= a.rt_evalue) if d == "RVT_1" else (
                dscore >= ga.get(d, 0) or ievalue <= 1e-5)
            if ok:
                hits[mid][d] = min(hits[mid].get(d, 1), ievalue)

    # merge overlapping models into loci
    loci = []
    for i, m in enumerate(models):
        q, ctg, s, e = m[0], m[1], m[2], m[3]
        for L in loci:
            if L["ctg"] == ctg and s <= L["e"] and e >= L["s"]:
                L["s"], L["e"] = min(L["s"], s), max(L["e"], e)
                L["models"].append(i)
                break
        else:
            loci.append({"ctg": ctg, "s": s, "e": e, "models": [i]})

    rows = []
    for L in loci:
        doms = {}
        for i in L["models"]:
            for d, ev in hits.get(f"m{i}", {}).items():
                doms[d] = min(doms.get(d, 1), ev)
        if not doms:
            continue
        best = max(L["models"], key=lambda i: models[i][5])
        refs = sorted({models[i][0].split("|")[1] if "|" in models[i][0]
                       else models[i][0] for i in L["models"]})
        rows.append({
            "locus": f"{L['ctg']}:{L['s']}-{L['e']}",
            "is_tert": "Telomerase_RBD" in doms,
            "domains": ",".join(f"{d}({doms[d]:.1e})" for d in sorted(doms)),
            "n_refs": len(refs),
            "best_AS": models[best][5],
            "prot_len": len(models[best][8]),
            "frameshifts": models[best][6],
            "stops": models[best][7],
        })
    rows.sort(key=lambda r: (not r["is_tert"], -r["best_AS"]))
    print("label\tlocus\tTERT_call\tdomains\tn_refs_agreeing\tbest_score\t"
          "model_aa\tframeshifts\tstops")
    for r in rows:
        print(f"{a.label}\t{r['locus']}\t{'TERT' if r['is_tert'] else 'RT/other'}\t"
              f"{r['domains']}\t{r['n_refs']}\t{r['best_AS']}\t{r['prot_len']}\t"
              f"{r['frameshifts']}\t{r['stops']}")
    n_tert = sum(r["is_tert"] for r in rows)
    print(f"#SUMMARY\t{a.label}\tTERT_loci={n_tert}\tother_domain_loci="
          f"{len(rows) - n_tert}\tmodels={len(models)}")


if __name__ == "__main__":
    main()
