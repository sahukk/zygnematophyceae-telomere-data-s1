set -e
cd /Users/krishnakumarsahu/Desktop/Algal_Telomere_DNA_Study/full_export/algal_telomere_study/data/processed/tr_search/iter2b
/opt/homebrew/Caskroom/miniconda/base/envs/tr_infernal/bin/cmcalibrate --cpu 10 iter2b.cm > cmcalibrate.log 2>&1; echo CALIBRATED > progress.log
for x in "Mend:../../../raw/ncbi/mesotaenium_endlicherianum/ncbi_dataset/data/GCA_009602735.1/GCA_009602735.1_ASM960273v1_genomic.fna" "Zcir:../../../raw/ncbi/zygnema_circumcarinatum/ncbi_dataset/data/GCA_037447605.1/GCA_037447605.1_SAG698-1b_1.0_genomic.fna" "Spra:../../../raw/phycocosm/Spipra1_AssemblyScaffolds_2024-03-24.fasta" "Zcyl:../../../raw/ncbi/GCA_037447655.1_dir/ncbi_dataset/data/GCA_037447655.1/GCA_037447655.1_SAG698-1a_XF_1.0_genomic.fna" "Clos:../../../raw/ncbi/GCA_027943415.1_dir/ncbi_dataset/data/GCA_027943415.1/GCA_027943415.1_CLOP_1.0_genomic.fna" "Spmu:../../../raw/ncbi/spirogloea_muscicola/ncbi_dataset/data/GCA_009602725.1/GCA_009602725.1_ASM960272v1_genomic.fna"; do
  l=${x%%:*}; g=${x#*:}
  [ -f "$g" ] || { echo "MISSING $l $g" >> progress.log; continue; }
  /opt/homebrew/Caskroom/miniconda/base/envs/tr_infernal/bin/cmsearch --cpu 10 -E 10 --tblout $l.tbl iter2b.cm $g > $l.out 2>&1
  echo "DONE $l $(grep -vc '^#' $l.tbl) hits" >> progress.log
done
echo ALL_DONE >> progress.log
