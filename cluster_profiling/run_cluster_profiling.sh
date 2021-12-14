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

# Define files used - data files
histology_file="${data_dir}/histologies.tsv"
count_file="${data_dir}/gene-counts-rsem-expected_count-collapsed.rds"

# Results file - this file is generated after manually reviewing results from 01
cc_data_match="${input_dir}/cc_data_match.tsv"

# references file 
short_long_match="${data_dir}/short_long_match.tsv"
gtf_file="${ref_dir}/gencode.v27.primary_assembly.annotation.gtf.gz"
gmt_file="${ref_dir}/c2.cp.kegg.v7.4.symbols.gmt"

# pathways of interest
pathway_file="${data_dir}/all_meta_pathways.tsv"
pathway_file_sub="${data_dir}/pathway_of_interest.tsv"


# Run cluster on each cancer group of interest
Rscript --vanilla 01-consensus_clustering.R \
--histology $histology_file \
--count $count_file \
--cg_interest "HGG,Medullo,EPN" \
--short_long_match $short_long_match \
--gtf_file $gtf_file \
--pathways $pathway_file

# Generate heatmap from results from `01-consensus_clustering.R`
Rscript --vanilla 02-cluster_heatmap.R \
--histology $histology_file \
--cc_data_match $cc_data_match \
--pathways $pathway_file_sub

# Run CEMiTools analysis on each cancer group of interest
Rscript --vanilla 03-network_analysis.R \
--count $count_file \
--cg_interest "HGG,Medullo,EPN" \
--short_long_match $short_long_match \
--gtf_file $gtf_file \
--gmt_file $gmt_file

# Run ssGSEA analysis on each cancer group and clusters
Rscript --vanilla 04-ssgsea_analysis.R \
--cg_interest "HGG,Medullo,EPN" \
--gmt_file $gmt_file

