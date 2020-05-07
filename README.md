# Natural selection shapes codon usage in the human genome
 Code repository for calculating synRVIS and synGERP scores

## Calculating synRVIS scores
Calculating synRVIS requires a list of synonymous variants (e.g. from gnomAD). Variants should be listed in a tab-delimited file with the following columns:  

| #CHROM  |   POS  | REF  | ALT |  GENE  | CODON_CHANGE | MAF |
| ------- | ------ | ---- | --- | ------ | ------------ | --- |
|1        | 861327 |   C  |  T  | SAMD11 |    tcC/tcT   |4.08180e-06


The script can then be run using the following command:
>Rscript src/syn_rvis.R \
>    -i \<input_path> \
>    -csc data/csc_scores.csv \
>    -maf \<maf_cutoff> \
>    -o \<output_path>


## Calculating synGERP scores
Calculating synGERP requires a list of GERP scores (http://mendel.stanford.edu/SidowLab/downloads/gerp/) and the human reference genome and GTF, both available from Ensembl. Calculating synGERP requires Python 3.6. This script will output per-gene synGERP scores, as well as the GERP scores for every fourfold degenerate synonymous site in canonical transcripts.

> python src/syn_gerp.py \
>    -e \<input.gtf> \
>    --rscu data/per_gene_rscu.txt \
>    -he data/hgnc_canonical_ens75.txt \
>    -r \<reference_genome.fa> \
>    -g \<gerp_scores_path> \
>    -o \<syn_gerp_output_path> \
>    -bo \<gerp_per_site_output_path>
