# Author: Run Jin
# Compare different modeling methods
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("DescTools")
  library("pec")
  library("riskRegression")
  library("survival")
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
cforest_input_dir <- file.path(analysis_dir, "results", "rf_party")

meta_combined_input_dir <- file.path(analysis_dir, "results", "expression_sub")

results_dir <- file.path(analysis_dir, "results", "model_method_comparison")
if(!dir.exists(results_dir)){
  dir.create(results_dir)
}

plots_dir <- file.path(analysis_dir, "plots", "model_method_comparison")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir)
}

#### read in relevant files 
for (i in 1:length(cg_list)){
  # find the cancer group of interest 
  x <- cg_list[i]
  
  # meta data directory
  meta_combined <- readr::read_tsv(file.path(meta_combined_input_dir, 
                                             paste0("meta_exp_combined_data_in_", x, ".tsv")))
  
  # get the gene variables 
  gene_variables <- meta_combined[, 7:37] %>% colnames()
  
  ############################### handle rfsrc
  # relevant results from rfsrc models
  rsfrc_grid_optimal <- readr::read_tsv(file.path(rfsrc_input_dir, x, "rfsrc_error_brier_lrs_results.tsv"))
  
  brier_rfsrc_fit <- readRDS(file.path(rfsrc_input_dir, x, "rfsrc_optimal_brier_output_full_data.RDS"))
  lrs_rfsrc_fit <- readRDS(file.path(rfsrc_input_dir, x, "rfsrc_optimal_lrs_output_full_data.RDS"))

  # get the optimal parameters for logrankscore method 
  mtry_min_lrs <- rsfrc_grid_optimal[[1,2]]
  ntree_min_lrs <- rsfrc_grid_optimal[[1,3]]
  node_min_lrs <- rsfrc_grid_optimal[[1,4]]
  nspilt_min_lrs <- rsfrc_grid_optimal[[1,5]]
  
  # get the optimal parameters for brier method 
  mtry_min_brier <- rsfrc_grid_optimal[[1,7]]
  ntree_min_brier <- rsfrc_grid_optimal[[1,8]]
  node_min_brier <- rsfrc_grid_optimal[[1,9]]
  nspilt_min_brier <- rsfrc_grid_optimal[[1,10]]
  
  ############################### read in cforest parameters
  # relevant results from cforest models
  cforest_grid <- readr::read_tsv(file.path(cforest_input_dir, 
                                            paste0("rfparty_hyperparam_grid_results_for_", x, ".tsv")))
  # select the mtry, ntree and mincriterion value with the highest c index
  cindex_max <- max(cforest_grid$c_index)
  mtry_max<- cforest_grid[cforest_grid$c_index == cindex_max, "mtry"] %>%
    pull(mtry) %>% unique() %>% sample(1)
  ntree_max <- cforest_grid[cforest_grid$c_index == cindex_max, "ntree"] %>% 
    pull(ntree) %>% unique() %>% sample(1)
  
  ############################### read in mboost results
  # read in mboost results 
  if(x == "LGG"){
    mboost_coxph <- readRDS(file.path(mboost_input_dir, paste0("coxph_pfs_model_fit_in_", x, ".rds"))) 
    mboost_loglog <- readRDS(file.path(mboost_input_dir, paste0("loglog_pfs_model_fit_in_", x, ".rds"))) 
  } else {
    mboost_coxph <- readRDS(file.path(mboost_input_dir, paste0("coxph_os_model_fit_in_", x, ".rds"))) 
    mboost_loglog <- readRDS(file.path(mboost_input_dir, paste0("loglog_os_model_fit_in_", x, ".rds"))) 
  }
  
  ##################### Define formula and subset relevant dataframes for downstream
  # extract variates for coxph from mboost coxph model
  coxph_variates <- coef(mboost_coxph) %>% 
    as.data.frame() %>%
    rownames()
  # extract only the relevant covariates
  coxph_variates <- coxph_variates[2:length(coxph_variates)]
  
  # extract variates for coxph from mboost loglog model
  loglog_variates <- coef(mboost_loglog) %>% 
    as.data.frame() %>%
    rownames()
  # extract only the relevant covariates
  loglog_variates <- loglog_variates[2:length(loglog_variates)]
  
  # piece together a model
  if(x == "LGG"){
    # model based on PFS days and status
    coxph.formula <- paste0("Surv(PFS_days, PFS_status_recode)", " ~ ",
                               paste0(coxph_variates, collapse = " + "),
                               sep = "")
    loglog.formula <- paste0("Surv(PFS_days, PFS_status_recode)", " ~ ",
                            paste0(loglog_variates, collapse = " + "),
                            sep = "")
    
    # model for all gene varialbes 
    total.formula <- paste0("Surv(PFS_days, PFS_status_recode)", " ~ ",
                            paste0(gene_variables, collapse = " + "),
                            sep = "")
    
    # make the combined data compatible with the analysis 
    combined_data_sub <- meta_combined %>% 
      dplyr::select(-c("Kids_First_Participant_ID", "OS_status",
                       "OS_status_recode", "OS_days", "PFS_status")) %>% 
      dplyr::filter(!is.na(PFS_days)) %>%
      tibble::column_to_rownames("Kids_First_Biospecimen_ID")
    
  } else {
    # model based on OS days and status
    coxph.formula <- paste0("Surv(OS_days, OS_status_recode)", " ~ ",
                               paste0(coxph_variates, collapse = " + "),
                               sep = "")
    loglog.formula <- paste0("Surv(OS_days, OS_status_recode)", " ~ ",
                             paste0(loglog_variates, collapse = " + "),
                             sep = "")
    
    # model for all gene variables 
    total.formula <- paste0("Surv(OS_days, OS_status_recode)", " ~ ",
                            paste0(gene_variables, collapse = " + "),
                            sep = "")
    
    # make the combined data compatible with the analysis 
    combined_data_sub <- meta_combined %>% 
      dplyr::select(-c("Kids_First_Participant_ID", "PFS_status",
                       "PFS_status_recode", "PFS_days", "OS_status")) %>% 
      dplyr::filter(!is.na(OS_days)) %>%
      tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  }
  
  ##################### Run CoxPH, cforest and LASSO cox for final step
  # generate a coxph model based on covaraites of interest
  coxph_mboost_model <- coxph(as.formula(coxph.formula),
                              data = combined_data_sub,
                              x = TRUE,
                              y = TRUE)
  
  loglog_mboost_model <- coxph(as.formula(loglog.formula),
                               data = combined_data_sub,
                               x = TRUE,
                               y = TRUE)
  
  # build a cforest with the best model
  fit.cforest <- pec::pecCforest(as.formula(total.formula), 
                                 data = combined_data_sub, 
                                 control = party::cforest_classical(ntree = ntree_max,
                                                                    mtry = mtry_max))
  
  # generate a glmnet fit for loglog mboost results
  # first a matrix with only covariates of interest
  combined_data_sub_glmnet <- combined_data_sub %>%
    dplyr::select(all_of(loglog_variates)) %>%
    as.matrix

  # generate another matrix with status and time and extract time and status as vector
  if(x == "LGG"){
    time_status_matrix <- combined_data_sub %>% dplyr::select("PFS_days", "PFS_status_recode") %>%
      dplyr::rename(time = PFS_days,
                    status = PFS_status_recode)
  } else {
    time_status_matrix <- combined_data_sub %>% dplyr::select("OS_days", "OS_status_recode") %>%
      dplyr::rename(time = OS_days,
                    status = OS_status_recode)
  }
  # find the time and status
  time <- time_status_matrix$time
  status <- time_status_matrix$status

  # fit the model when there are two or more covariates
  if(ncol(combined_data_sub_glmnet)>=2){
    # run glmnet
    coxlasso_model <- glmnet(x = combined_data_sub_glmnet,
                             y = survival::Surv(time, status),
                             family = 'cox')
    coxlasso_cv_model <- cv.glmnet(x = combined_data_sub_glmnet,
                             y = survival::Surv(time, status),
                             family = 'cox')
    
    # run survfit
    fit_coxlasso <- survfit(coxlasso_model,
                            x = combined_data_sub_glmnet,
                            y = survival::Surv(time, status))
    
    fit_cv_coxlasso <- survfit(coxlasso_cv_model,
                               x = combined_data_sub_glmnet,
                              y = survival::Surv(time, status))
  }
  
  ##################### Generate final output using pec package
  # generate an overall pec model
  pec <- pec::pec(object = list("CoxPH_mboost" = coxph_mboost_model, 
                                "Loglog_mboost" = loglog_mboost_model,
                                # "Loglog_mboost_glmnet" = fit_coxlasso, 
                                # "Loglog_mboost_cv_glmnet" = fit_cv_coxlasso, 
                                "rfsrc_lrs" = lrs_rfsrc_fit, 
                                "rfsrc_brier" = brier_rfsrc_fit,
                                "cforst" = fit.cforest),
                 formula = as.formula(total.formula), 
                 data = combined_data_sub, 
                 exact = TRUE,
                 splitMethod = "BootCv", 
                 B = 500)
  
  # save the results 
  pec %>%
    saveRDS(file.path(results_dir, paste0("pec_pred_error_comp_results_in_", x, ".rds")))
  
  # plot out the figures 
  pdf(file.path(plots_dir, paste0("pec_pred_error_comp_results_in_", x, ".pdf")))
  plot(pec)
  dev.off()
}

