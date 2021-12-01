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
count_file="${data_dir}/gene-counts-rsem-expected_count-collapsed.rds"
short_long_match="${data_dir}/short_long_match.tsv"
cg_gene_interest="${data_dir}/cg_gene_interest_sub.tsv"
gtf_file="${ref_dir}/gencode.v27.primary_assembly.annotation.gtf.gz"

# Run GSNCA analysis per cancer group of interest
Rscript --vanilla geometric_mean_expr.R \
--histology $histology_file \
--count $count_file \
--cg_gene_interest $cg_gene_interest \
--short_long_match $short_long_match \
--gtf_file $gtf_file 

