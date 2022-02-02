# Author: Run Jin
# Compare different modeling methods
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("DescTools")
  library("pec")
  library("riskRegression")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"))
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cg_list <-unlist(strsplit(opt$cg_interest,","))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "relative_importance")
mboost_input_dir <- file.path(analysis_dir, "results", "mboost", "model_fit")
rfsrc_input_dir <- file.path(analysis_dir, "results", "rfsrc")
meta_combined_input_dir <- file.path(analysis_dir, "results", "expression_sub")

results_dir <- file.path(analysis_dir, "results", "model_method_comparison")
if(!dir.exists(results_dir)){
  dir.create(results_dir)
}

#### read in relevant files 
for (i in 1:length(cg_list)){
  # find the cancer group of interest 
  x <- cg_list[i]
  
  # meta data directory
  meta_combined <- readr::read_tsv(file.path(meta_combined_input_dir, 
                                             paste0("meta_exp_combined_data_in_", x, ".tsv"))) %>%
    dplyr::select(Kids_First_Biospecimen_ID, OS_days, PFS_days) %>%
    dplyr::filter(!is.na(OS_days))
  
  # relevant results from rfsrc models
  brier_rfsrc_fit <- readRDS(file.path(rfsrc_input_dir, x, "rfsrc_optimal_brier_output_full_data.RDS"))
  lrs_rfsrc_fit <- readRDS(file.path(rfsrc_input_dir, x, "rfsrc_optimal_lrs_output_full_data.RDS"))
  
  # read in mboost results
  if(x == "LGG"){
    coxph_mboost_fit <- readRDS(file.path(mboost_input_dir, paste0("coxph_pfs_model_fit_in_", x, ".rds"))) 
    loglog_mboost_fit <- readRDS(file.path(mboost_input_dir, paste0("loglog_pfs_model_fit_in_", x, ".rds"))) 
  } else {
    coxph_mboost_fit <- readRDS(file.path(mboost_input_dir, paste0("coxph_os_model_fit_in_", x, ".rds"))) 
    loglog_mboost_fit <- readRDS(file.path(mboost_input_dir, paste0("loglog_os_model_fit_in_", x, ".rds"))) 
  }
  
  
}

