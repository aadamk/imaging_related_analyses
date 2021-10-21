#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Define directory
data_dir="../data"
ref_dir="../references"
results_dir="results"

# Define files used
histology_file="${data_dir}/histologies.tsv"
expression_file="${data_dir}/gene-expression-rsem-tpm-collapsed.rds"
short_long_match="${data_dir}/short_long_match.tsv"
cg_gene_interest="${data_dir}/cg_gene_interest.tsv"
gtf_file="${ref_dir}/gencode.v27.primary_assembly.annotation.gtf.gz"
gsea_score="${results_dir}/gsea_scores/combined_gsea_scores.tsv"
gsea_pval="${results_dir}/gsea_scores/combined_gsea_pval.tsv"

# Run GSNCA analysis per cancer group of interest
Rscript --vanilla 01-gsnca_analysis_per_gene_disease.R \
--histology $histology_file \
--expression $expression_file \
--cg_gene_interest $cg_gene_interest \
--short_long_match $short_long_match \
--gtf_file $gtf_file

# Calculate GSEA scores per cancer group of interest
Rscript --vanilla 02-ssgsea_analysis_per_gene_disease.R \
--histology $histology_file \
--expression $expression_file \
--cg_gene_interest $cg_gene_interest \
--short_long_match $short_long_match \
--gtf_file $gtf_file \
--out_file $gsea_score

# Calculate statistical significance of GSEA scores for each group
Rscript --vanilla 03-ssgsea_sig_calc_per_gene_disease.R \
--gsea_score_file $gsea_score \
--cg_gene_interest $cg_gene_interest \
--out_file $gsea_pval
