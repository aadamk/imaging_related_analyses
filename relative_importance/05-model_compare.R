# Author: Run Jin
# Compare different modeling methods
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("DescTools")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"))
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cg_list <-unlist(strsplit(opt$cg_interest,","))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "relative_importance")
mboost_input_dir <- file.path(analysis_dir, "results", "mboost", "predicted_scores")
rfsrc_input_dir <- file.path(analysis_dir, "results", "rfsrc")
meta_combined_input_dir <- file.path(analysis_dir, "results", "expression_sub")

results_dir <- file.path("results", "somers_delta")
if(!dir.exists(results_dir)){
  dir.create(results_dir)
}

# define a matrix for storing results 
somer_results <- as.data.frame(matrix(nrow = 4, ncol = length(cg_list)))
colnames(somer_results) <- cg_list
rownames(somer_results) <- c("rfsrc_brier", "rfsrc_logrankscore", "mboost_coxph", "mboost_loglog")

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
  brier_predicted <- readRDS(file.path(rfsrc_input_dir, x, "rfsrc_optimal_brier_output_full_data.RDS"))[["predicted"]] %>%
    # get predicted results
    as.data.frame() 
  colnames(brier_predicted) <- "predicted_score"

  # get lrs results
  lrs_predicted <- readRDS(file.path(rfsrc_input_dir, x, "rfsrc_optimal_lrs_output_full_data.RDS"))[["predicted"]] %>%
    as.data.frame()
  colnames(lrs_predicted) <- "predicted_score"
  
  # read in coxph results
  if(x == "LGG"){
    coxph <- readr::read_tsv(file.path(mboost_input_dir, paste0("coxph_pfs_pred_risk_in_", x, ".tsv"))) %>% 
      left_join(meta_combined)
    loglog <- readr::read_tsv(file.path(mboost_input_dir, paste0("loglog_pfs_pred_risk_in_", x, ".tsv"))) %>% 
      left_join(meta_combined)
  } else {
    coxph <- readr::read_tsv(file.path(mboost_input_dir, paste0("coxph_os_pred_risk_in_", x, ".tsv"))) %>% 
      left_join(meta_combined)
    loglog <- readr::read_tsv(file.path(mboost_input_dir, paste0("loglog_os_pred_risk_in_", x, ".tsv"))) %>% 
      left_join(meta_combined)

  }
  
  ####### calculate Somers Delta
  # brier
  rownames(brier_predicted) <- meta_combined$Kids_First_Biospecimen_ID
  brier_predicted <- brier_predicted %>% 
    tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>% 
    left_join(meta_combined)
  
  # lrs
  rownames(lrs_predicted) <- meta_combined$Kids_First_Biospecimen_ID
  lrs_predicted <- lrs_predicted %>% 
    tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>% 
    left_join(meta_combined)
  
  ###### calculate scores
  if(x == "LGG"){
    stat_brier = DescTools::SomersDelta(as.matrix(brier_predicted[,c(2,4)]))
    stat_lrs = DescTools::SomersDelta(as.matrix(lrs_predicted[,c(2,4)]))
    stat_coxph = DescTools::SomersDelta(as.matrix(coxph[,c(2,4)]))
    stat_loglog = DescTools::SomersDelta(as.matrix(loglog[,c(2,4)]))
  } else {
    stat_brier = DescTools::SomersDelta(as.matrix(brier_predicted[,c(2,3)]))
    stat_lrs = DescTools::SomersDelta(as.matrix(lrs_predicted[,c(2,3)]))
    stat_coxph = DescTools::SomersDelta(as.matrix(coxph[,c(2,3)]))
    stat_loglog = DescTools::SomersDelta(as.matrix(loglog[,c(2,4)]))
  }
  
  #### store the results 
  somer_results[1,i] <- stat_brier
  somer_results[2,i] <- stat_lrs
  somer_results[3,i] <- stat_coxph
  somer_results[4,i] <- stat_loglog
}

# output somers results 
somer_results %>%
  tibble::rownames_to_column("model_method") %>%
  readr::write_tsv(file.path(results_dir, "somers_delta_results.tsv"))
