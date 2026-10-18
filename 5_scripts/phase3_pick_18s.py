"""Phase 3, step 2: choose the best 18S locus per assembly from the BLAST hits
written by phase3_extract_18s.py, extract it with blastdbcmd, and summarise.

Locus = cluster of same-contig, same-strand HSPs within 4 kb. Locus score =
best single-query summed bitscore. The extracted span is the best query's
HSP span extended by the query ends left uncovered (so partial HSPs at the
divergent 5'/3' ends do not truncate the gene).
"""
import csv, subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/processed/phylogeny_phase3"

EXCLUDE = {"EF553458.1", "KU351743.1"}  # carry ITS/28S flanks -> match whole rDNA units


def groups(rows):
    """Collinear HSP groups: same query, contig, strand, and a consistent
    query->subject offset (within 1 kb). Each group is one rDNA copy."""
    out = []
    for r in sorted(rows, key=lambda r: -r["bitscore"]):
        plus = r["sstart"] < r["send"]
        off = (r["sstart"] - r["qstart"]) if plus else (r["sstart"] + r["qstart"])
        for g in out:
            if g["key"] == (r["qseqid"], r["sseqid"], plus) and abs(g["off"] - off) < 1000:
                g["hsps"].append(r); break
        else:
            out.append({"key": (r["qseqid"], r["sseqid"], plus), "off": off, "hsps": [r]})
    for g in out:
        g["score"] = sum(h["bitscore"] for h in g["hsps"])
    return sorted(out, key=lambda g: -g["score"])


def span(g):
    hs, plus = g["hsps"], g["key"][2]
    qlen, slen = hs[0]["qlen"], hs[0]["slen"]
    qs, qe = min(r["qstart"] for r in hs), max(r["qend"] for r in hs)
    ss = min(min(r["sstart"], r["send"]) for r in hs)
    se = max(max(r["sstart"], r["send"]) for r in hs)
    if plus:
        ss -= qs - 1; se += qlen - qe
    else:
        ss -= qlen - qe; se += qs - 1
    cov = (qe - qs + 1) / qlen
    pid = sum(r["pident"] * r["length"] for r in hs) / sum(r["length"] for r in hs)
    return max(1, ss), min(slen, se), round(cov, 3), round(pid, 1)


def main():
    summary, fasta = [], []
    for name, _ in csv.reader(open(ROOT / "data/processed/telomere_reverify/assembly_manifest.tsv"), delimiter="\t"):
        f = OUT / "hits" / f"{name}.tsv"
        rows = [dict(zip("qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen".split(), l.split("\t")))
                for l in f.read_text().splitlines()] if f.exists() else []
        for r in rows:
            for k in ("pident", "evalue", "bitscore"): r[k] = float(r[k])
            for k in ("length", "qstart", "qend", "sstart", "send", "qlen", "slen"): r[k] = int(r[k])
        rows = [r for r in rows if r["qseqid"] not in EXCLUDE]
        gs = groups(rows)
        if not gs:
            summary.append([name, "NO_HIT"] + [""] * 9); continue
        g = gs[0]
        q, sid, plus = g["key"]; st = "+" if plus else "-"
        ss, se, cov, pid = span(g)
        n_copies = sum(1 for x in gs if x["key"][0] == q and x["score"] >= 0.8 * g["score"])
        seq = subprocess.run(["blastdbcmd", "-db", str(OUT / "blastdb" / name), "-entry", sid,
                              "-range", f"{ss}-{se}", "-strand", "plus" if plus else "minus"],
                             capture_output=True, text=True, check=True).stdout.split("\n", 1)[1].replace("\n", "")
        fasta.append(f">{name}\n{seq}\n")
        summary.append([name, "OK", sid, ss, se, st, len(seq), q, cov, pid, n_copies])
    (OUT / "extracted/genome_18s.fasta").write_text("".join(fasta))
    with open(OUT / "extraction_summary.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["assembly", "status", "contig", "start", "end", "strand", "len", "best_query",
                    "query_cov", "pident", "near_best_loci"])
        w.writerows(summary)

if __name__ == "__main__":
    main()
