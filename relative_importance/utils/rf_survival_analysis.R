# Function: Random Forest model 
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(survival)
  library(randomForestSRC)
  library(ggRandomForests)
  library(randomForest)
  library(caret)
  library(ggpubr)
  library(Hmisc)
  library(pec)
})

#' Function: Use metadata, and optionally cluster annotation to run rfsrc survival analysis
#' on covariates of interest
#' 
#' @author Adam Kraya and Run Jin
#'
#' @param metadata metadata file containing survival information and covariates of interest
#' @param ind_var independent variables to run survival analysis on - can be continuous or categorical
#' @param parse_perc parsing ratio for training and testing set 
#' E.g., 0.6 means 60% going to training and 40% as testing - default is 0.5 
#' @param metadata_sample_col column name of the column (in metadata) that specifies sample names 
#' corresponding to the sample column in the cluster annotation file
#' @param metadata_indep_col column name of the column (in metadata) that specifies individual participants
#' survival analysis is performed on independent participant level
#' @param os_days_col column name of the column (in metadata) that specifies time of survival
#' @param os_status_col column name of the column (in metadata) that specifies event of survival
#' This column can contain either `LIVING` or `DECEASED` or NA
#' @param results_dir output directory for fit summary
#' @param plots_dir output directory for survival plots, risk tables and diagnostic plots
#' 
#' @return a dataframe containing grid search results for logrank, another with grid seach results for brier,
#' a dataframe containing hyperparameters for minimal err.rate in the grid search
#' 
#'
#' @examples rf_survival_analysis(metadata = metadata,
#'                                ind_var= c("broad_histology", "CNS_region"),
#'                                results_dir = "results_test",
#'                                plots_dir= "plots_test")
#'
#' @export
rf_survival_analysis <- function(metadata,
                                 ind_var,
                                 parse_perc = 0.5,
                                 metadata_sample_col = "Kids_First_Biospecimen_ID",
                                 metadata_indep_col = "Kids_First_Participant_ID",
                                 os_days_col = "OS_days",
                                 os_status_col = "OS_status",
                                 results_dir,
                                 plots_dir){
  
  if (is.null(ind_var)) {
    stop("Error: No variable has been supplied to test with using the `ind_var` argument. ")
  }
  # List the columns we need
  needed_cols <- c(os_status_col, os_days_col, metadata_sample_col, metadata_indep_col) %>% unique()
  
  # Get logical vector indicating which are in metadata
  if(!all(needed_cols %in% colnames(metadata))){
    stop("Error: Provided column names for OS status and OS days are not present in metadata. Please check.")
  }
  
  # make directory
  if (!dir.exists(results_dir)) {
    dir.create(results_dir,recursive = TRUE)
  }
  
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir,recursive = TRUE)
  }
  
  # first rename the OS days and OS status column
  metadata <- metadata %>% 
    dplyr::rename(OS_days = !!os_days_col,
                  OS_status = !!os_status_col)
  if(!all((metadata$OS_status) %in% c("LIVING", "DECEASED", NA))){
    stop("Error: OS status column can only contain `LIVING` or `DECEASED`. Please modify. ")
  }
  
  if(!is.null(metadata_sample_col)){
    # modify the colnames of the metadata for easier transformation later on
    metadata <- metadata %>% 
      dplyr::rename(sample = !!metadata_sample_col,
                    participant = !!metadata_indep_col) %>%
      # remove non-independent samples 
      dplyr::distinct(participant, .keep_all = TRUE)
  } else {
    # modify the colnames of the metadata for easier transformation later on
    metadata <- metadata %>% 
      dplyr::rename(sample = !!metadata_sample_col) 
  }
  
  # Reformat as a numeric variable where 0 = LIVING and 1 = DECEASED.
  metadata <- metadata %>%
    dplyr::mutate(OS_status_recode = case_when(
      OS_status == "LIVING" ~ 0,
      OS_status == "DECEASED" ~ 1
    ))
  
  # format input data to be compatible with the functions below 
  input_data <- metadata %>%
    dplyr::filter(!is.na(OS_days))
  
  # Set formula for rfsrc
  # tune hyperparameters for the rfsrc implementation
  # run the grid search of mtry and ntree and find the optimal hyperparameters
  
  set.seed(2021)
  # decide the partitioning method for the run
  trainIndex <- createDataPartition(input_data$OS_status, p = parse_perc, 
                                    list = FALSE, 
                                    times = 1)
  # separate data into different 
  data.train <- input_data[trainIndex, ] 
  data.test <- input_data[-trainIndex, ] 

  # piece together a model
  rfsrc.formula <- paste0("Surv(OS_days, OS_status_recode)", " ~ ",
                          paste0(ind_var, collapse = " + "),
                          sep = "")
  rfsrc.formula <- as.formula(rfsrc.formula)
  
  # grid search for the brier method
  k=1
  rfsrc_pbc_brier <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(rfsrc_pbc_brier) <- c("err.rate", "mtry", "ntree", "nodesize", "nsplit")
  
  for(i in 1:length(ind_var)){
    for(j in seq(from =100, to = 2000, by = 100)){
      for(p in c(5, 10, 20, 25, 30)){
        for(q in c(0,1)){
          rfsrc_pbc_brier_rf <- randomForestSRC::rfsrc(rfsrc.formula, 
                                                       data = as.data.frame(data.train), 
                                                       mtry=i, 
                                                       ntree=j, 
                                                       nodesize=p,
                                                       nsplit=q,
                                                       splitrule = "bs.gradient",
                                                       na.action = "na.impute", 
                                                       tree.err = TRUE,
                                                       importance = TRUE,
                                                       forest=T)
         
          # write up results5
          rfsrc_pbc_brier[k,1] <- mean(rfsrc_pbc_brier_rf[["err.rate"]], na.rm=TRUE)
          rfsrc_pbc_brier[k,2] <- i
          rfsrc_pbc_brier[k,3] <- j
          rfsrc_pbc_brier[k,4] <- p
          rfsrc_pbc_brier[k,5] <- q
          k <- k+1
        }
      }
    }
  }

  error_brier_min <- min(rfsrc_pbc_brier$err.rate, na.rm=TRUE)
  mtry_brier_min<- rfsrc_pbc_brier[rfsrc_pbc_brier$err.rate == error_brier_min, "mtry"]
  ntree_brier_min<- rfsrc_pbc_brier[rfsrc_pbc_brier$err.rate == error_brier_min, "ntree"]

  nodesize_brier_min<- rfsrc_pbc_brier[rfsrc_pbc_brier$err.rate == error_brier_min, "nodesize"]
  nsplit_brier_min<- rfsrc_pbc_brier[rfsrc_pbc_brier$err.rate == error_brier_min, "nsplit"]


  # grid search for the method logrankscore
  k=1
  rfsrc_pbc_lrs <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(rfsrc_pbc_lrs) <- c("err.rate", "mtry", "ntree", "nodesize", "nsplit")

  # for(i in 1:length(ind_var)){
  for(i in 1:length(ind_var)){
    for(j in seq(from =100, to = 2000, by = 100)){
      for(p in c(5, 10, 20, 25, 30)){
        for(q in c(0,1)){
          rfsrc_pbc_lrs_rf <- randomForestSRC::rfsrc(rfsrc.formula,
                                                    data = as.data.frame(data.train),
                                                    mtry=i,
                                                    ntree=j,
                                                    nodesize=p,
                                                    nsplit=q,
                                                    splitrule = "logrankscore",
                                                    na.action = "na.impute",
                                                    tree.err = TRUE,
                                                    importance = TRUE,
                                                    forest=T)

          # write up results
          rfsrc_pbc_lrs[k,1] <- mean(rfsrc_pbc_lrs_rf[["err.rate"]], na.rm=TRUE)
          rfsrc_pbc_lrs[k,2] <- i
          rfsrc_pbc_lrs[k,3] <- j
          rfsrc_pbc_lrs[k,4] <- p
          rfsrc_pbc_lrs[k,5] <- q
          k <- k+1
        }
      }
    }
  }

  error_lrs_min <- min(rfsrc_pbc_lrs$err.rate)
  mtry_lrs_min<- rfsrc_pbc_lrs[rfsrc_pbc_lrs$err.rate == error_lrs_min, "mtry"]
  ntree_lrs_min<- rfsrc_pbc_lrs[rfsrc_pbc_lrs$err.rate == error_lrs_min, "ntree"]

  nodesize_lrs_min<- rfsrc_pbc_lrs[rfsrc_pbc_lrs$err.rate == error_lrs_min, "nodesize"]
  nsplit_lrs_min<- rfsrc_pbc_lrs[rfsrc_pbc_lrs$err.rate == error_lrs_min, "nsplit"]

  # generate a combined table with all the minimal parameter
  combined <- data.table('error_lrs_min' = error_lrs_min,
                         'mtry_of_lrs_min' = mtry_lrs_min,
                         'ntree_of_lrs_min' = ntree_lrs_min,
                         'nodesize_of_lrs_min' = nodesize_lrs_min,
                         'nsplit_of_lrs_min'= nsplit_lrs_min,
                         'error_brier_min' = error_brier_min,
                         'mtry_of_brier_min' = mtry_brier_min,
                         'ntree_of_brier_min' = ntree_brier_min,
                         'nodesize_of_brier_min' = nodesize_brier_min,
                         'nsplit_of_brier_min' = nsplit_brier_min)

  readr::write_tsv(rfsrc_pbc_brier, file.path(results_dir, "rfsrc_error_mtryntree_ns_split_grid_brier.tsv"))
  readr::write_tsv(rfsrc_pbc_lrs, file.path(results_dir, "rfsrc_error_mtryntree_ns_split_grid_logrankscore.tsv"))
  readr::write_tsv(combined, file.path(results_dir, "rfsrc_error_brier_lrs_results.tsv"))
  
  # plots the results to see how it looks - brier as splitrule
  # fit the model using the optimal mtry and ntree
  rfsrc_pbc_brier_rf <- randomForestSRC::rfsrc(rfsrc.formula,
                                               data = as.data.frame(data.train),
                                               mtry=mtry_brier_min,
                                               ntree=ntree_brier_min,
                                               nodesize=nodesize_brier_min,
                                               nsplit=nsplit_brier_min,
                                               splitrule = "bs.gradient",
                                               na.action = "na.impute",
                                               tree.err = TRUE,
                                               importance = TRUE,
                                               forest=T)

  lapply(ind_var, function(x){
    pdf(file = file.path(plots_dir, paste0("rf_brier_output_by_",x,".pdf")),
        h = height,
        w = width,
        onefile = T)
    print(plot(na.omit(ggRandomForests::gg_error(rfsrc_pbc_brier_rf))))
    print(plot(ggRandomForests::gg_rfsrc(rfsrc_pbc_brier_rf, by = x)) +
            labs(x = "Time (days)"))
    print(plot(ggRandomForests::gg_vimp(rfsrc_pbc_brier_rf)))
    dev.off()
  })
  
  # plots the results to see how it looks - logrankscore as splitrule
  # fit the model using the optimal mtry and ntree
  rfsrc_pbc_lrs_rf <- randomForestSRC::rfsrc(rfsrc.formula,
                                             data = as.data.frame(data.train),
                                             mtry=mtry_lrs_min,
                                             ntree=ntree_lrs_min,
                                             nodesize=nodesize_lrs_min,
                                             nsplit=nsplit_lrs_min,
                                             splitrule = "logrankscore",
                                             na.action = "na.impute",
                                             tree.err = TRUE,
                                             importance = TRUE,
                                             forest=T)
  
  lapply(ind_var, function(x){
    pdf(file = file.path(plots_dir, paste0("rf_logrankscore_output_by_",x,".pdf")),
        h = height,
        w = width,
        onefile = T)
    print(plot(na.omit(ggRandomForests::gg_error(rfsrc_pbc_lrs_rf))))
    print(plot(ggRandomForests::gg_rfsrc(rfsrc_pbc_lrs_rf, by = x)) +
            labs(x = "Time (days)"))
    print(plot(ggRandomForests::gg_vimp(rfsrc_pbc_lrs_rf)))
    dev.off()
  })

  ################ Test the model - and output C index for different time points
  # filter the test set to contain only entries that do not have NA in any covariate
  data.test_na_fixed <- data.test %>% filter_at(vars(ind_var),all_vars(!is.na(.)))
  
  # fit the brier model to the NA fixed testing data
  rfsrc_pbc_brier.test <- randomForestSRC::predict.rfsrc(rfsrc_pbc_brier_rf,
                                                         newdata = as.data.frame(data.test_na_fixed),
                                                         importance = T,
                                                         na.action = 'na.impute',
                                                         proximity = T,
                                                         forest.wt = T,
                                                         outcome = 'train',
                                                         seed = 2021)
  
  # compute the apparent estimate of the C-index at different time points
  ApparrentCindex  <- pec::cindex(object=rfsrc_pbc_brier.test[["survival"]],
                                  formula=rfsrc.formula,
                                  data=as.data.frame(data.test_na_fixed),
                                  eval.times=rfsrc_pbc_brier.test[["time.interest"]])
  
  # get the c index for time of interest
  brier_cindex <- ApparrentCindex[["AppCindex"]][["matrix"]] %>% as.data.frame()
  colnames(brier_cindex) <- "cindex"
  
  # output the results
  brier_cindex %>% 
    dplyr::mutate(day_of_interest = rfsrc_pbc_brier.test[["time.interest"]]) %>%
    readr::write_tsv(file.path(results_dir, "rfsrc_brier_cindex_per_time_interest.tsv"))

  # fit the logrank model to the NA fixed testing data
  rfsrc_pbc_lrs.test <- randomForestSRC::predict.rfsrc(rfsrc_pbc_lrs_rf,
                                                       newdata = as.data.frame(data.test_na_fixed),
                                                       importance = T,
                                                       na.action = 'na.impute',
                                                       proximity = T,
                                                       forest.wt = T,
                                                       outcome = 'train',
                                                       seed = 2021)
  # compute the apparent estimate of the C-index at different time points
  ApparrentCindex  <- pec::cindex(object=rfsrc_pbc_lrs.test[["survival"]],
                                  formula=rfsrc.formula,
                                  data=as.data.frame(data.test_na_fixed),
                                  eval.times=rfsrc_pbc_lrs.test[["time.interest"]])
  
  # get the c index for time of interest
  lrs_cindex <- ApparrentCindex[["AppCindex"]][["matrix"]] %>% as.data.frame()
  colnames(lrs_cindex) <- "cindex"
  
  # output the results
  lrs_cindex %>% 
    dplyr::mutate(day_of_interest = rfsrc_pbc_lrs.test[["time.interest"]]) %>%
    readr::write_tsv(file.path(results_dir, "rfsrc_lrs_cindex_per_time_interest.tsv"))
}
