"""Phase 3: extract full-length nuclear 18S rRNA from each genome assembly.

For each assembly: makeblastdb, dc-megablast a reference 18S panel, take the
best-scoring locus, extend to the full query span, and export it (plus a
copy-count / identity summary). Output: data/processed/phylogeny_phase3/.
"""
import csv, subprocess, sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/processed/phylogeny_phase3"
PANEL = OUT / "ref/panel.fa"
MANIFEST = ROOT / "data/processed/telomere_reverify/assembly_manifest.tsv"

def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

def process(name, fasta):
    db = OUT / "blastdb" / name
    if not Path(str(db) + ".njs").exists() and not Path(str(db) + ".nal").exists():
        run(["makeblastdb", "-in", str(ROOT / fasta), "-dbtype", "nucl", "-out", str(db),
             "-max_file_sz", "4GB", "-parse_seqids"])
    hits = OUT / "hits" / f"{name}.tsv"
    if not hits.exists():
        run(["blastn", "-task", "dc-megablast", "-query", str(PANEL), "-db", str(db),
             "-evalue", "1e-30", "-num_threads", "2", "-max_target_seqs", "50",
             "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen",
             "-out", str(hits)])
    return name

if __name__ == "__main__":
    jobs = [r for r in csv.reader(open(MANIFEST), delimiter="\t") if r]
    with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as ex:
        for n in ex.map(lambda r: process(*r), jobs):
            print("done", n, flush=True)
