#!/bin/bash

##After first run, consider running particular SNPs or genes:
## --multiGene=GENE Restrict analysis to this gene (for analyses that support this feature)

## And summarising by LD with:
##--ldOnly 
	## Compute LD blocks for existing eQTL results. This assumes that previous results can be loaded from the 
	## file implied by '--output'. Implies '--ldBlocks'
##--ldBlocks=TRUE
##--ldPairs
		##Flag indicating whether eQTLs should be grouped based on pairwise measures of LD. 
## And see options --maxDist=MAXDIST --maxSNPs=MAXSNPS	--ldFDR=LDFDR	--ldPvalR2=LDPVALR2 --ldR2=LDR2

geno_file=cut_genotype_data_all_treated_baseline.tsv_matched.tsv
gex_file=cut_GEx_baseline_4000_and_2000.tsv_matched.tsv
snp_pos=snp146Common_MatrixEQTL_snp_pos.txt
gene_pos=biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt

output_file=mePipe_all_treated_baseline
##output_file_trans= #sprintf('trans_%s_%s.txt', threshold, geno_file)
##output_file_cis= #sprintf('cis_%s_%s_%s.txt', cis, cisThreshold, geno_file)

##useModel=linear # default is linear, use -a option to change to anova
##pthreshold=100000 #1e-05
cisThreshold=0.001
cisdist=1000000 #1e+06

EXCLUDECOV=0.001
FILTERTHRESHOLD=0.001
COVTHRESHOLD=0.001
multiPvalue=1000000 #1e-06
SGEOPTIONS=FALSE
##--pthreshold=$pthreshold \

/Library/Frameworks/R.framework/Resources/library/mePipe/exec/runMatrixEQTL.R \
--output=$output_file \
--snpspos=$snp_pos \
--genepos=$gene_pos \
--cisthreshold=$cisThreshold \
--cisdist=$cisdist \
--snpspos=$snp_pos \
--genepos=$gene_pos \
--selectcov \
--filtercov \
--excludecov=EXCLUDECOV \
--filterthreshold=FILTERTHRESHOLD \
--covthreshold=COVTHRESHOLD \
--effectSize \
--multiPeak \
--multiPvalue=$multiPvalue \
--bins=1000 \
--qqplot \
--verbose \
--sgeoptions=SGEOPTIONS \
$gex_file \
$geno_file

