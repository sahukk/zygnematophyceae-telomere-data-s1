# TR hit (iter2 coordinates) +/- 500 bp flank, oriented to the TR sense strand.
# Output: tr_regions.fa ; header = species|TRstart_in_region|TRend_in_region
G={"Mend":("data/raw/ncbi/mesotaenium_endlicherianum/ncbi_dataset/data/GCA_009602735.1/GCA_009602735.1_ASM960273v1_genomic.fna",("WELA01010214.1",818526,818792,"+")),
"Zcir":("data/raw/ncbi/zygnema_circumcarinatum/ncbi_dataset/data/GCA_037447605.1/GCA_037447605.1_SAG698-1b_1.0_genomic.fna",("CM074816.1",1257292,1257573,"+")),
"Zcyl":("data/raw/ncbi/GCA_037447655.1_dir/ncbi_dataset/data/GCA_037447655.1/GCA_037447655.1_SAG698-1a_XF_1.0_genomic.fna",("JBBHKE010000708.1",21618,21895,"-"))}
F=500
rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
with open("data/processed/tr_search/rnaseq/tr_regions.fa","w") as out:
    for sp,(g,(c,a,b,st)) in G.items():
        seq=[];on=False
        for l in open(g):
            if l[0]==">":
                if on: break
                on=l[1:].split()[0]==c
            elif on: seq.append(l.strip())
        s="".join(seq).upper(); lo=max(1,a-F); hi=min(len(s),b+F)
        r=s[lo-1:hi]
        if st=="-": r=rc(r); t0=hi-b; t1=hi-a+1
        else: t0=a-lo; t1=b-lo+1
        out.write(f">{sp}|{t0}|{t1}|{c}:{lo}-{hi}{st}\n{r}\n")
        print(sp,len(r),"TR at",t0,t1)
