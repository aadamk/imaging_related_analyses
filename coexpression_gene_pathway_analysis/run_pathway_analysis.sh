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
ssgsea_score="${results_dir}/ssgsea_scores/combined_ssgsea_scores.tsv"
ssgsea_pval="${results_dir}/ssgsea_pval/combined_ssgsea_pval.tsv"
dcga_go_out="${results_dir}/dgca_go_term/combined_dgca_go_term_pval.tsv"
gscna_pval="${results_dir}/gsnca_pval/combined_gscna_pathway_analysis.tsv"

# Run GSNCA analysis per cancer group of interest
Rscript --vanilla 01-gsnca_analysis_per_gene_disease.R \
--histology $histology_file \
--expression $expression_file \
--cg_gene_interest $cg_gene_interest \
--short_long_match $short_long_match \
--gtf_file $gtf_file \
--outfile $gscna_pva

# Use DGCA to calculte co expression changes
Rscript --vanilla 02-dgca_analysis_per_gene_disease.R \
--histology $histology_file \
--expression $expression_file \
--cg_gene_interest $cg_gene_interest \
--short_long_match $short_long_match \
--gtf_file $gtf_file \
--outfile $dcga_go_out


# Calculate ssGSEA scores per cancer group of interest
Rscript --vanilla 03-ssgsea_analysis_per_gene_disease.R \
--histology $histology_file \
--expression $expression_file \
--cg_gene_interest $cg_gene_interest \
--short_long_match $short_long_match \
--gtf_file $gtf_file \
--out_file_score $ssgsea_score \
--out_file_pval $ssgsea_pval


# Generate heatmap for the top 20 most altered pathways
Rscript --vanilla 04-barplot_altered_pathways.R \
--cg_gene_interest $cg_gene_interest \
--ssgsea $ssgsea_pval \
--gsnca $gscna_pval

