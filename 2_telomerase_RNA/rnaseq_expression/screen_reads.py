"""Stream FASTQ on stdin; keep reads sharing >=2 25-mers with one species' TR region (+/-500 bp).
Writes TSV: read_name, strand(+ = same as TR sense), first_pos, last_pos (region coords), n_kmer_hits, sequence.
Usage: fasterq-dump -Z ... | python3 screen_reads.py tr_regions.fa SPECIES out.tsv"""
import sys
K=25
fa,sp,out=sys.argv[1:4]
recs=open(fa).read().split(">")[1:]
seq=[r.split("\n",1)[1].replace("\n","") for r in recs if r.startswith(sp+"|")][0]
rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
idx={}
for i in range(len(seq)-K+1):
    idx.setdefault(seq[i:i+K],("+",i))
    idx.setdefault(rc(seq[i:i+K]),("-",i))
n=kept=0
with open(out,"w") as o:
    o.write("read\tstrand\tstart\tend\tnk\tseq\n")
    it=iter(sys.stdin)
    for h in it:
        s=next(it).strip(); next(it); next(it); n+=1
        hits=[idx[s[j:j+K]] for j in range(0,len(s)-K+1,6) if s[j:j+K] in idx]
        if len(hits)>=2:
            st=max("+-",key=lambda x:sum(1 for h2 in hits if h2[0]==x))
            ps=[p for x,p in hits if x==st]
            o.write(f"{h.split()[0][1:]}\t{st}\t{min(ps)}\t{max(ps)+K}\t{len(hits)}\t{s}\n"); kept+=1
print(f"{sp}: {n} reads screened, {kept} kept",file=sys.stderr)
