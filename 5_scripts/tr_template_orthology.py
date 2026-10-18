#!/usr/bin/env python3
"""TR candidates by template + cross-species conservation.

Step 1 (per genome): find every occurrence of the template-length window
(C-rich telomere unit + 3 nt, all circular permutations; e.g. CTAAACCCTA-type
for TTTAGGG) on either strand, outside tandem telomeric arrays, and write it
with -300/+300 nt of flank, oriented so the C-rich template reads 5'->3'.

Step 2 (all genomes together): all-vs-all blastn of those flanks between
species. Telomerase RNA genes are conserved between related species beyond
the template itself, whereas chance template matches are not. Candidates
with a significant hit in another species, outside the template window, are
reported, grouped into ortholog clusters.

Usage: tr_template_orthology.py OUTDIR LABEL:GUNIT:genome.fasta [...]
"""
import collections
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_telomere_ends import read_fasta, revcomp  # noqa: E402

outdir = sys.argv[1]
os.makedirs(outdir, exist_ok=True)
FL = 300
allfa = os.path.join(outdir, "template_flanks_all.fasta")
summary = []
with open(allfa, "w") as fo:
    for spec in sys.argv[2:]:
        label, gunit, genome = spec.split(":", 2)
        cunit = revcomp(gunit)
        L = len(cunit)
        wins = {(cunit * 3)[i:i + L + 3] for i in range(L)}
        pat = re.compile("(?=(" + "|".join(sorted(wins)) + "))")
        two = re.compile("|".join(sorted({(cunit * 3)[i:i + 2 * L] for i in range(L)})))
        n_all = n_kept = 0
        for name, s in read_fasta(genome):
            c = name.split()[0]
            s = s.upper()
            for strand, seq in (("+", s), ("-", revcomp(s))):
                for m in pat.finditer(seq):
                    p = m.start()
                    n_all += 1
                    if two.search(seq[max(0, p - 2 * L):p + 3 * L + 3]):
                        continue
                    flank = seq[max(0, p - FL):p + L + 3 + FL]
                    if len(flank) < 2 * FL or flank.count("N") > 20:
                        continue
                    n_kept += 1
                    fo.write(f">{label}|{c}|{strand}|{p}|{m.group(1)}\n{flank}\n")
        summary.append(f"{label}\t{gunit}\ttemplate windows {n_all}\toutside arrays {n_kept}")

db = os.path.join(outdir, "template_flanks_all")
subprocess.run(["makeblastdb", "-in", allfa, "-dbtype", "nucl", "-out", db], check=True,
               capture_output=True)
hits = subprocess.run(["blastn", "-query", allfa, "-db", db, "-evalue", "1e-8", "-dust", "yes",
                       "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
                       "-max_target_seqs", "50"], capture_output=True, text=True, check=True).stdout

T0, T1 = FL, FL + 15      # template window inside each flank (approx.)
pairs = collections.defaultdict(list)
for line in hits.splitlines():
    q, sid, pid, ln, qs, qe, ss, se, ev, bs = line.split("\t")
    if q.split("|")[0] == sid.split("|")[0]:
        continue                      # same species: paralogs/repeats, skip
    qs, qe = int(qs), int(qe)
    # require the alignment to extend clearly beyond the template window
    outside = max(0, min(qe, T0) - qs) + max(0, qe - max(qs, T1))
    if outside < 40:
        continue
    pairs[q].append((sid, float(pid), int(ln), qs, qe, float(ev)))

with open(os.path.join(outdir, "orthology_report.txt"), "w") as fo:
    fo.write("\n".join(summary) + "\n\n")
    fo.write("candidate\tn_other_species_hits\tbest_hits (subject pident len qstart-qend evalue)\n")
    for q in sorted(pairs, key=lambda k: -len({h[0].split('|')[0] for h in pairs[k]})):
        spp = {h[0].split("|")[0] for h in pairs[q]}
        best = sorted(pairs[q], key=lambda h: h[5])[:4]
        fo.write(f"{q}\t{len(spp)}\t" + "; ".join(f"{h[0]} {h[1]:.0f}% {h[2]}bp {h[3]}-{h[4]} {h[5]:.0e}" for h in best) + "\n")
print("\n".join(summary))
print(f"candidates with cross-species conservation beyond the template: {len(pairs)}")
