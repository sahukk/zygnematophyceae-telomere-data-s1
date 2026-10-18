#!/usr/bin/env python3
"""Telomerase RNA (TR) gene candidates, following the USE + template strategy
of Fajkus et al. 2021 (NAR 49:7680, doi:10.1093/nar/gkab545).

Plant/algal TRs are RNA polymerase III genes with a type-3 promoter: an
Upstream Sequence Element (USE) followed by a TATA box. The same promoter
drives U6 snRNA genes, which are highly conserved and easy to find, so:

1. U6 genes are found with blastn (query: a plant U6 snRNA), and the 150 bp
   upstream of each is searched for the most over-represented 8-mers. The top
   8-mer (allowing 1 mismatch) is taken as this genome's USE.
2. Minimal template candidates are every circular permutation of the C-rich
   telomere unit plus one nucleotide (unit length + 1), found on either
   strand, excluding positions inside tandem telomeric arrays.
3. A candidate is kept when a USE match lies 30–200 nt upstream on the same
   strand. Each kept candidate is written with -300/+400 nt of flank, for
   cross-species comparison (true TR orthologs share sequence beyond the
   template; chance hits do not).

Usage: tr_candidates.py LABEL genome.fasta G_STRAND_UNIT u6_query.fasta OUTDIR
"""
import collections
import os
import random
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_telomere_ends import read_fasta, revcomp  # noqa: E402

label, genome, gunit, u6q, outdir = sys.argv[1:6]
os.makedirs(outdir, exist_ok=True)
seqs = {n.split()[0]: s.upper() for n, s in read_fasta(genome)}


def get(chrom, start, end, strand):
    s = seqs[chrom]
    start, end = max(0, start), min(len(s), end)
    sub = s[start:end]
    return sub if strand == "+" else revcomp(sub)


# ---- 1. U6 genes and USE discovery ---------------------------------------
db = os.path.join(outdir, f"{label}.db")
if not os.path.exists(db + ".nsq") and not os.path.exists(db + ".00.nsq"):
    subprocess.run(["makeblastdb", "-in", genome, "-dbtype", "nucl", "-out", db], check=True, capture_output=True)
res = subprocess.run(["blastn", "-query", u6q, "-db", db, "-evalue", "1e-10",
                      "-outfmt", "6 sseqid sstart send length pident", "-max_target_seqs", "500",
                      "-word_size", "11"], capture_output=True, text=True, check=True).stdout
ups = []
for line in res.splitlines():
    sid, ss, se, ln, pid = line.split("\t")
    ss, se = int(ss), int(se)
    if int(ln) < 60:
        continue
    if ss < se:   # gene on + strand; upstream = before ss
        ups.append(get(sid, ss - 1 - 150, ss - 1, "+"))
    else:
        ups.append(get(sid, ss, ss + 150, "-"))
ups = list({u for u in ups if len(u) == 150})

K = 8
pres = collections.Counter()
for u in ups:
    pres.update({u[i:i + K] for i in range(len(u) - K + 1)})
# background: random 150-bp windows
rnd = random.Random(1)
chroms = [c for c in seqs if len(seqs[c]) > 10000]
bg = collections.Counter()
for _ in range(3000):
    c = rnd.choice(chroms)
    p = rnd.randrange(0, len(seqs[c]) - 150)
    w = seqs[c][p:p + 150]
    bg.update({w[i:i + K] for i in range(len(w) - K + 1)})
ranked = sorted(pres, key=lambda k: -(pres[k] / max(len(ups), 1)) / ((bg[k] + 1) / 3000))
ranked = [k for k in ranked if pres[k] >= max(2, 0.4 * len(ups)) and "N" not in k]
use_kmers = ranked[:3]

report = [f"# {label}: {len(ups)} U6 upstream regions; telomere G-unit {gunit}"]
for k in ranked[:8]:
    report.append(f"USE 8-mer {k}: in {pres[k]}/{len(ups)} U6 upstreams; background {bg[k]}/3000 windows")


def mism(a, b):
    return sum(x != y for x, y in zip(a, b))


def use_hit(up):
    """Best USE match (<=1 mismatch to any top 8-mer) in an upstream region."""
    best = None
    for k in use_kmers:
        for i in range(len(up) - K + 1):
            m = mism(up[i:i + K], k)
            if m <= 1 and (best is None or m < best[0]):
                best = (m, i, k)
    return best


# ---- 2. minimal template candidates -------------------------------------
cunit = revcomp(gunit)
L = len(cunit)
templ = {(cunit * 3)[i:i + L + 1] for i in range(L)}
pat = re.compile("(?=(" + "|".join(sorted(templ)) + "))")
two = re.compile("|".join(sorted({(cunit * 3)[i:i + 2 * L] for i in range(L)})))

cands = []
for c, s in seqs.items():
    for strand, seq in (("+", s), ("-", revcomp(s))):
        for m in pat.finditer(seq):
            p = m.start()
            ctx = seq[max(0, p - L):p + 2 * L + 1]
            if two.search(ctx):       # inside a tandem telomeric array
                continue
            up = seq[max(0, p - 200):max(0, p - 30)]
            h = use_hit(up) if use_kmers else None
            cands.append((c, strand, p, m.group(1), h, seq[max(0, p - 300):p + 400]))

kept = [x for x in cands if x[4]]
report.append(f"template sites (not in arrays): {len(cands)}; with USE 30-200 nt upstream: {len(kept)}")
with open(os.path.join(outdir, f"{label}_TR_candidates.fasta"), "w") as fo:
    for c, strand, p, t, h, flank in kept:
        dist = 200 - h[1]
        fo.write(f">{label}|{c}|{strand}|{p}|template={t}|USE={h[2]}mm{h[0]}@-{dist}\n{flank}\n")
open(os.path.join(outdir, f"{label}_report.txt"), "w").write("\n".join(report) + "\n")
print("\n".join(report))
