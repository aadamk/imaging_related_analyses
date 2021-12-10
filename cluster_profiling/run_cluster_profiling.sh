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
input_dir="input"
ref_dir="../references"
results_dir="results"

# Define files used
histology_file="${data_dir}/histologies.tsv"
count_file="${data_dir}/gene-counts-rsem-expected_count-collapsed.rds"
short_long_match="${data_dir}/short_long_match.tsv"
gtf_file="${ref_dir}/gencode.v27.primary_assembly.annotation.gtf.gz"
pathway_file="${data_dir}/all_meta_pathways.tsv"

cc_data_match="${input_dir}/cc_data_match.tsv"


# Run cluster on each cancer group of interest 
Rscript --vanilla 01-consensus_clustering.R \
--histology $histology_file \
--count $count_file \
--cg_interest "HGG,Medullo,EPN" \
--short_long_match $short_long_match \
--gtf_file $gtf_file \
--pathways $pathway_file

Rscript --vanilla 02-network-analysis.R \
--histology $histology_file \
--cc_data_match $cc_data_match \
--short_long_match $short_long_match \
--gtf_file $gtf_file \
--pathways $pathway_file


