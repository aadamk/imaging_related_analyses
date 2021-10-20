# Author: Run Jin
# Adapted from https://github.com/PediatricOpenTargets/OpenPedCan-analysis/blob/
# gsea-pull-update/analyses/gene-set-enrichment-analysis/02-model-gsea.Rmd
  
### Purpose
# The purpose of this analysis is to assess significant differences in GSVA scores for each pathway. 
# Using ANOVA and subsequent Tukey tests, we do:
# - For each unique cancer_group -gene combination in the GSVA scores
#   - For each pathway, are GSVA scores significantly different across `group`?
#     - If so, which cancer_group -gene combination are significantly different?
# We perform this using GSVA scores calculated from `02-ssgsea_analysis_per_gene_disease.R`.

BiocManager::install("BiocParallel")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library(broom))

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-o","--gsva_score_file"),type="character",
              help="GSVA scores file of each pathway in each disease (.tsv)"),
  make_option(c("-l","--cg_gene_interest"),type="character",
              help="file containing gene of interest and matching cancer group (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "coexpression_gene_pathway_analysis")
results_dir <- file.path(analysis_dir, "results", "gsva_anova_tukey")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}

#### Read in files necessary for analyses --------------------------------------
# GSVA scores
gsva_score_file <- readr::read_tsv(opt$gsva_score_file)

# This script contains functions used to modeling GSVA scores
source(file.path("util", "pathway_models.R"))

# Significance testing universal threshold
SIGNIFICANCE_THRESHOLD <- 0.01

########## Define output file for each unique cancer_group - gene combination

cancer_group_anova_outpaths <- lapply(rna_library_list, function(x){
  x<-gsub(" ", "_", x)
  x<-stringr::str_to_lower(gsub("-", "", x))
  file.path(results_dir, paste0("gsva_anova_", x, "_cancer_group.tsv"))
})
cancer_group_tukey_outpaths <- lapply(rna_library_list, function(x){
  x<-gsub(" ", "_", x)
  x<-stringr::str_to_lower(gsub("-", "", x))
  file.path(results_dir, paste0("gsva_tukey_", x, "_cancer_group.tsv"))
})


