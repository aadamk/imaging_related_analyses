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
cancer_group_list="LGG,HGG,EPN,Medullo"

# Run CosPH survival analysis for all genes of interest in all cancer group of interest
Rscript --vanilla 01-mboost_survival.R \
--histology $histology_file \
--expression $expression_file \
--cg_interest $cancer_group_list \
--short_long_match $short_long_match 

# Calculate ssGSEA scores per cancer group of interest
Rscript --vanilla 02-ssgsea_analysis.R \
--histology $histology_file \
--expression $expression_file \
--cg_interest $cancer_group_list \
--short_long_match $short_long_match 