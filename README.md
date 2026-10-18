# Data S1: telomere repeat diversity in Zygnematophyceae

Supporting data and scripts for:

> Sahu, K. K. *Telomere repeat diversity in Zygnematophyceae, the closest algal
> relatives of land plants: evidence from long and short reads.* Manuscript
> submitted to *Journal of Phycology* (2026).

All input data are public; accessions are listed in Table 1, Table S1 and
Table S2 of the article. Scripts contain absolute paths from the author's
computer; change them to local paths before re-running.

## Licence

Data and documentation: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Scripts in `5_scripts/` and the `*.sh` / `*.py` files elsewhere: MIT (see `LICENSE`).

## Contents

```
1_telomere_unit_tables/
  Per-strain tables of telomeric units in telomeric arrays at long-read ends
  (long_*.txt; PacBio HiFi/CLR or ONT) or in array-containing Illumina reads
  (short_*.txt): counts of units by T-tract length (data for Figure 2 and
  Table 1). run_fig2_data.sh = commands used to produce them.
  Abbreviations: Clos = Closterium; Mend = Mesotaenium endlicherianum;
  Spmu = Spirogloea muscicola; Spra = Spirogyra pratensis;
  Zcir = Zygnema circumcarinatum; Zcyl = Zygnema cf. cylindricum.

2_telomerase_RNA/
  Zygnema_TR_sequences.fasta     telomerase RNA genes of Z. circumcarinatum
                                 and Z. cf. cylindricum (genome coordinates
                                 in the headers)
  TR_alignment_streptophyte_plus_Zygnema.stk
                                 structural alignment with published
                                 streptophyte telomerase RNAs
  fajkus2021_streptophyte_TRs.fasta  published seed sequences (Fajkus et al. 2021)
  iter2b.stk / iter2b.cm         seed alignment and Infernal covariance model
                                 (Fajkus et al. 2021 seed + Zygnema hits)
  zygnematophyte_TR_hits.fa      sequences added to the seed
  cmsearch_hits/*.tbl            cmsearch results for the six genomes
  run_iter2b.sh                  commands (cmbuild, cmcalibrate, cmsearch)
  rnaseq_expression/             reads matching each telomerase RNA locus
                                 +/- 500 bp in RNA-seq runs (Zcir total RNA:
                                 SRR24286561, SRR24205694; Mend and Zcyl:
                                 poly(A) libraries), with region file and
                                 screening scripts

3_TERT/
  *.models.faa / *.domtbl        TERT candidate protein models (miniprot) and
                                 their HMMER domain hits for the six taxa
  Zcir_TERT.faa                  Z. circumcarinatum TERT protein
  Zcyl_TRBD_orf.faa, Zcyl_RT_model.faa
                                 Z. cf. cylindricum TERT, split over two
                                 scaffolds (TRBD part and RT/CTE part)
  TERT_calls_all_assemblies.tsv  TERT calls for the assembly set (Table S1)

4_18S_phylogeny/
  18S_alignment_trimmed.fasta    trimmed 18S rRNA alignment (Figure S1)
  tree.* / constrained.treefile  IQ-TREE maximum-likelihood trees
                                 (unconstrained and constrained to the
                                 published topology in constraint.nwk)
  au_test.iqtree                 approximately unbiased (AU) topology test
  run.sh                         commands

5_scripts/
  Python scripts for read-end telomere detection, telomeric unit counting,
  assembly end-versus-interior tests, TERT search, telomerase RNA candidate
  search and template comparison, 18S extraction, Table S1 and figures.
```
