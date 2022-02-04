# Author: Run Jin
# Output one representative tree using the party package
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("party")
  library("Hmisc")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"),
  make_option(c("-e","--path_to_meta"),type="character",
              help="path to the file directory that contains combined meta and gene expression data")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cg_list <-unlist(strsplit(opt$cg_interest,","))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "relative_importance")

plots_dir <- file.path(analysis_dir, "plots", "rf_party")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

results_dir <- file.path(analysis_dir, "results", "rf_party")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

meta_exp_input_dir <- opt$path_to_meta

# run the tree generation function for all cancer groups 
for (i in 1:length(cg_list)){
  x <- cg_list[i]
  
  # read in the combined data with metadata + gene expression 
  meta_exp_combined <- readr::read_tsv(file.path(meta_exp_input_dir, 
                                                 paste0("meta_exp_combined_data_in_", x, ".tsv")))
  
  # get the gene variables 
  gene_variables <- meta_exp_combined[, 7:37] %>% colnames()
  
  ########### prepare the formula and data based on which cancer group we are looking at
  if(x == "LGG"){
    # piece together a model
    survival.formula <- paste0("Surv(PFS_days, PFS_status_recode)", " ~ ",
                                paste0(gene_variables, collapse = " + "),
                                sep = "")
    # make the combined data compatible with the analysis 
    combined_data_sub <- meta_exp_combined %>% 
      dplyr::select(-c("Kids_First_Participant_ID", "OS_status",
                       "OS_status_recode", "OS_days", "PFS_status")) %>% 
      dplyr::filter(!is.na(PFS_days)) %>%
      tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  
  } else{
    # piece together a model
    survival.formula <- paste0("Surv(OS_days, OS_status_recode)", " ~ ",
                                paste0(gene_variables, collapse = " + "),
                                sep = "")
    # make the combined data compatible with the analysis 
    combined_data_sub <- meta_exp_combined %>% 
      dplyr::select(-c("Kids_First_Participant_ID", "PFS_status",
                       "PFS_status_recode", "PFS_days", "OS_status")) %>% 
      dplyr::filter(!is.na(OS_days)) %>%
      tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  }
  
  # run the grid search of mtry and ntree and find the optimal setting
  k=1
  cforest_tune_df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(cforest_tune_df) <- c("c_index", "somer_dxy", "mtry", "ntree", "mincriterion")
  
  # find the best parameters for cforest
  for(m in 1:length(gene_variables)){
    for(n in seq(from = 100, to = 2000, by = 100)){
      for(j in seq(from = 0.5, to = 1, by = 0.05)){
        
        # define the control object using hyperparameter
        cforest_ctrl_obj <- cforest_control(teststat = "max",
                                            testtype = "Teststatistic",
                                            mincriterion = j,
                                            savesplitstats = FALSE,
                                            ntree = n, 
                                            mtry = m, 
                                            replace = TRUE,
                                            trace = FALSE)
        
        # run the cforest fit with the hyperparameters
        cforest_fit <- cforest(as.formula(survival.formula), 
                               data = as.data.frame(combined_data_sub), 
                               subset = NULL, 
                               weights = NULL,
                               controls = cforest_ctrl_obj,
                               xtrafo = ptrafo, 
                               ytrafo = ptrafo, 
                               scores = NULL)
        
        # run cindex
        if(sum(is.finite(predict(cforest_fit, type = 'response', OOB = T)))==nrow(combined_data_sub)){
          if(x == "LGG"){
            cor_index <- rcorr.cens(predict(cforest_fit, type = 'response', OOB = T),
                                    with(combined_data_sub,
                                         Surv(PFS_days,PFS_status_recode)))
          } else {
            cor_index <- rcorr.cens(predict(cforest_fit, type = 'response', OOB = T),
                                    with(combined_data_sub,
                                         Surv(PFS_days,PFS_status_recode)))
          }
        }
        
        # store the results
        cforest_tune_df[k,1] <- cor_index[["C Index"]] %>% as.numeric()
        cforest_tune_df[k,2] <- cor_index[["Dxy"]] %>% as.numeric()
        cforest_tune_df[k,3] <- m
        cforest_tune_df[k,4] <- n
        cforest_tune_df[k,5] <- j
        k <- k+1
      }
    }
  }
  
  # write out the results 
  cforest_tune_df %>% 
    dplyr::arrange(c_index) %>%
    readr::write_tsv(file.path(results_dir, paste0("rfparty_hyperparam_grid_results_for_", x, ".tsv")))
  
  # find the optimal parameter
  # select the mtry, ntree and mincriterion value with the highest c index
  cindex_max <- max(cforest_tune_df$c_index)
  mtry_max<- cforest_tune_df[cforest_tune_df$c_index == cindex_max, "mtry"]
  mincri_max<- cforest_tune_df[cforest_tune_df$c_index == cindex_max, "mincriterion"]
  
  ####### run the analysis with optimal parameter
  ctree_ctrl_obj <- ctree_control(teststat = "max",
                                  testtype = "Teststatistic",
                                  mincriterion = mincri_max,
                                  minbucket = 5,
                                  stump = FALSE,
                                  maxsurrogate = 0,
                                  mtry = mtry_max, 
                                  savesplitstats = TRUE, 
                                  maxdepth = 0, 
                                  remove_weights = FALSE)
  
  # get the results
  ctree_results <- ctree(as.formula(survival.formula), 
                         combined_data_sub, 
                         subset = NULL,
                         weights = NULL,
                         controls = ctree_ctrl_obj, 
                         xtrafo = ptrafo, 
                         ytrafo = ptrafo,
                         scores = NULL)
  
  # plot the tree and save the plot 
  pdf(file.path(plots_dir, paste0("rf_party_plot_in_", x, ".pdf")))
  plot(ctree_results)
  dev.off()

}
