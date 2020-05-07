# Natural selection shapes codon usage in the human genome
 Code repository for calculating synRVIS and synGERP scores

## Calculating synRVIS scores
Calculating synRVIS requires a list of synonymous variants (e.g. from gnomAD). Variants should be listed in a tab-delimited file with the following columns:  
>#CHROM	POS	REF	ALT	GENE	CODON_CHANGE	MAF


The script can then be run using the following command:
>Rscript src/syn_rvis.R \
>    -i \<input_path> \
>    -csc data/csc_scores.csv \
>    -maf \<maf_cutoff>
>    -o \<output_path>


## Calculating synGERP scores