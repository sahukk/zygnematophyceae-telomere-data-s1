# B1: is the TR gene expressed? Download one or two RNA-seq runs per species, keep only reads matching TR +/- 500 bp.
cd /Users/krishnakumarsahu/Desktop/Algal_Telomere_DNA_Study/full_export/algal_telomere_study
O=data/processed/tr_search/rnaseq; R=data/raw/rnaseq_tr; mkdir -p $R/tmp
one(){ sp=$1; run=$2
  prefetch --max-size 20g -O $R $run > $O/$run.prefetch.log 2>&1 || { echo "FAIL prefetch $sp $run"; return; }
  fasterq-dump -Z --split-spot -t $R/tmp -e 2 $R/$run/$run.sra* 2>$O/$run.fqd.log | python3 -I $O/screen_reads.py $O/tr_regions.fa $sp $O/${sp}_$run.tsv 2>>$O/$run.fqd.log
  echo "DONE $sp $run $(tail -1 $O/$run.fqd.log)"
}
one Zcir SRR24286561 &   # total RNA (rRNA-depleted), SAG 698-1b, dark 2h
one Zcir SRR24205694 &   # total RNA (rRNA-depleted), SAG 698-1b, 4C desiccation
one Mend SRR31697631 &   # control, cDNA library
one Zcyl SRR37105910 &   # polyA
wait; echo ALL_DONE
