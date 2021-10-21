# Author: Run Jin
# Adapted from https://github.com/PediatricOpenTargets/OpenPedCan-analysis/blob/
# gsea-pull-update/analyses/gene-set-enrichment-analysis/02-model-gsea.Rmd
  
### Purpose
# The purpose of this analysis is to assess significant differences in GSEA scores for each pathway. 
# Using ANOVA and subsequent Tukey tests, we do:
# - For each unique cancer_group -gene combination in the GSEA scores
#   - For each pathway, are GSEA scores significantly different across `group`?
#     - If so, which cancer_group -gene combination are significantly different?
# We perform this using GSEA scores calculated from `02-ssgsea_analysis_per_gene_disease.R`.

BiocManager::install("BiocParallel")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library(broom))

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-o","--gsea_score_file"),type="character",
              help="GSEA scores file of each pathway in each disease (.tsv)"),
  make_option(c("-l","--cg_gene_interest"),type="character",
              help="file containing gene of interest and matching cancer group (.tsv)"),
  make_option(c("-o","--out_file"),type="character",
              help="output file with pval of GSEA scores of each pathway in each disease (.tsv)"),
  
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
outfile <- opt$out_file

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "coexpression_gene_pathway_analysis")
results_dir <- file.path(analysis_dir, "results", "gsea_stats")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}

#### Read in files necessary for analyses --------------------------------------
# GSEA scores
gsea_scores <- readr::read_tsv(opt$gsea_score_file)
# unique cancer_group and gene of interest combination
cg_gene_interest <- readr::read_tsv(opt$cg_gene_interest)

########## Define output file for each unique cancer_group - gene combination
combined_results <- data.frame()
for (i in 1:nrow(cg_gene_interest)){
  # find the cancer group of interest
  cg_interest <- cg_gene_interest[i,1] %>% 
    pull(short_name)
  # find the gene of interest
  gene_interest <- cg_gene_interest[i,2] %>% 
    pull(gene_of_interest)
  
  gsea_score_each <- gsea_scores %>% 
    filter(cancer_group == cg_interest && gene_parsed_by == gene_interest) 
  
  # gather list of each pathway 
  pathway_list <- gsea_score_each %>% 
    pull(pathway) %>% unique
  
  gsea_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(gsea_results) <- c("pathway_id", "description","pval")
  for(j in 1:length(pathway_list)){
    pathway_interest <- pathway_list[j]
    
    # get description of pathway
    description <- gsea_score_each %>%  filter(pathway == pathway_interest) %>% 
      pull(description) %>% unique()
    
    # filter to contain only pathway of interest
    gsea_score_each_per_pathway <- gsea_score_each %>% 
      filter(pathway == pathway_interest)
    # group1 scores
    group1 <- gsea_score_each_per_pathway %>% filter(group == "1") %>% 
      pull(gsea_score)
    # group2 scores
    group2 <- gsea_score_each_per_pathway %>% filter(group == "2") %>% 
      pull(gsea_score)
    
    # compute t.test of the two
    res <- t.test(group1, group2)$p.value
    gsea_results[j,1] <- pathway_interest
    gsea_results[j,2] <- description
    gsea_results[j,3] <- res
    
    # write out individual file
    readr::write_tsv(gsea_results, file.path(results_dir, paste0(cg_interest, "_parsed_by_", gene_interest, "_gsea_pval.tsv" )))
  }
  # annotate gene and cancer group to the result
  gsea_results <- gsea_results %>% 
    mutate(cancer_group = cg_interest) %>% 
    mutate(gene_parsed_by = gene_interest) 
  # combine them to a combined tsv file
  combined_results <- rbind(combined_results, gsea_results)
}

# write out combined results
readr::write_tsv(gsea_results, outfile)
