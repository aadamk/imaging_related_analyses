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
stat_dir="results"

# Define files used
histology_file="${data_dir}/histologies.tsv"
expression_file="${data_dir}/gene-expression-rsem-tpm-collapsed.rds"
short_long_match="${data_dir}/short_long_match.tsv"
cg_gene_interest="${data_dir}/cg_gene_interest.tsv"
gtf_file="${ref_dir}/gencode.v27.primary_assembly.annotation.gtf.gz"

# Define other information
gene_list="SLC7A5,FOLH1,BRAF,NRAS,PEBP1,MAPK1,MAPK3,MAP2K1,MAP2K2,SLC1A5,SLC3A2,SLC7A11,SLC6A14"
cancer_group_list="LGG,HGG,EPN,Medullo,Cranio,ATRT"

# Run GSNCA analysis per cancer group of interest
Rscript --vanilla GSNCA_analysis_per_disease.R \
--histology $histology_file \
--expression $expression_file \
--cg_gene_interest $cg_gene_interest \
--short_long_match $short_long_match \
--gtf_file $gtf_file


