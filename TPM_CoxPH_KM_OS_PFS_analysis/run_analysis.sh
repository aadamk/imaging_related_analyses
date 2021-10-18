#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit


data_dir="../data"
histology_file="${data_dir}/histologies.tsv"
expression_file="${data_dir}/gene-expression-rsem-tpm-collapsed.rds"
short_long_match="${data_dir}/short_long_match.tsv"

ref_dir="../../references"
stat_dir="results"
stat_outfile="${stat_dir}/km_stats_combined.tsv"
stat_outfile_cox="${stat_dir}/coxph_test_stats_combined.tsv"

gene_list="SLC7A5,FOLH1,BRAF,NRAS,PEBP1,MAPK1,MAPK3,MAP2K1,MAP2K2,SLC1A5,SLC3A2,SLC7A11,SLC6A14"
cancer_group_list="LGG,HGG,EPN,Medullo,Cranio,ATRT"

# Plot out the TPM vs. harmonized diagnosis for all genes of interest in all cancer group of interest
Rscript --vanilla 01-TPM_vs_harmonized_diag.R \
--histology $histology_file \
--expression $expression_file \
--cancer_groups $cancer_group_list \
--gene_list $gene_list \
--short_long_match $short_long_match


# Rscript -e "rmarkdown::render('02-lgg_cns_gtex.Rmd', clean = TRUE)"

# Run KM survival analysis for all genes of interest in all cancer group of interest
Rscript --vanilla 03-km_w_logrank_survival.R \
--histology $histology_file \
--expression $expression_file \
--cancer_groups $cancer_group_list \
--gene_list $gene_list \
--short_long_match $short_long_match \
--stat_outfile $stat_outfile

# Run CosPH survival analysis for all genes of interest in all cancer group of interest
Rscript --vanilla 04-coxph_reg_survival.R \
--histology $histology_file \
--expression $expression_file \
--cancer_groups $cancer_group_list \
--gene_list $gene_list \
--short_long_match $short_long_match \
--stat_outfile $stat_outfile_cox


