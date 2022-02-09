# Author: Adam Kraya, Run Jin
# Compare different modeling methods
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("DescTools")
  library("survival")
  library("glmnet")
  library("survminer")
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
meta_combined_input_dir <- file.path(analysis_dir, "results", "expression_sub")
results_dir <- file.path(analysis_dir, "results", "model_method_comparison")

if(!dir.exists(results_dir)){
  dir.create(results_dir)
}

plots_dir <- file.path(analysis_dir, "plots", "model_method_comparison")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir)
}

for (i in 1:length(cg_list)){
  # find the cancer group of interest 
  x <- cg_list[i]
  
  # meta data directory
  meta_combined <- readr::read_tsv(file.path(meta_combined_input_dir, 
                                             paste0("meta_exp_combined_data_in_", x, ".tsv")))
  
  # get the gene variables 
  gene_variables <- meta_combined[, 7:37] %>% colnames()
  if(x == "LGG"){
    
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
  opt.param.lamda1se = list()
  opt.param.lamda_min = list()
  k=1
  
  if(ncol(combined_data_sub)>=2){
    # run elastic net cox regression, searching for optimal combination of alphas and lambdas
    alphas = seq(from = 0, to = 1, by = 0.05)
    opt.param.lamda_min = sapply(alphas, function(x){
      coxlasso_cv_model <- cv.glmnet(x = as.matrix(combined_data_sub[,gene_variables]),
                                     y = survival::Surv(time, status),
                                     family = 'cox',
                                     type.measure = 'deviance',
                                     alpha = x)

      c('lamda.min' = coxlasso_cv_model$lambda.min, 'partial.likelihood.lmin' = max(coxlasso_cv_model$cvm))})
    
    opt.param.lamda1se = sapply(alphas, function(x){
      coxlasso_cv_model <- cv.glmnet(x = as.matrix(combined_data_sub[,gene_variables]),
                                     y = survival::Surv(time, status),
                                     family = 'cox',
                                     type.measure = 'deviance',
                                     alpha = x)

      c( 'lamda.1se' = coxlasso_cv_model$lambda.1se, 'partial.likelihood.1se' = min(coxlasso_cv_model$cvm))})
    
    opt.param.lamda_min = as.data.table(t(opt.param.lamda_min))
    opt.param.lamda1se = as.data.table(t(opt.param.lamda1se))
    
    # Find maximum partial likelihood by lambda.min and lambda.1se methods
    ind.lambda.min = which.max(opt.param.lamda_min$partial.likelihood.lmin)
    ind.lambda1se = which.max(opt.param.lamda1se$partial.likelihood.1se)
    
    alpha.min = alphas[ind.lambda.min]
    alpha.1se = alphas[ind.lambda1se]
    
    
    ## Compute model coefficients with optimal alpha and lambda
    ## IMPORTANT: this can be re-called in 05-model_compare.R using the optimal alpha and lambda
    ## that are written to a json file. These models below will be the input for pec. 
    
    opt.model.min = glmnet(x = as.matrix(combined_data_sub[,gene_variables]),
                           y = survival::Surv(time, status),
                           family = 'cox',
                           type.measure = 'deviance',
                           alpha = alpha.min,
                           lambda = opt.param.lamda_min$lamda.min[ind.lambda.min])
    
    opt.model.1se = glmnet(x = as.matrix(combined_data_sub[,gene_variables]),
                           y = survival::Surv(time, status),
                           family = 'cox',
                           type.measure = 'deviance',
                           alpha = alpha.1se,
                           lambda = opt.param.lamda1se$lamda.1se[ind.lambda1se])
    
    
    ### Model selected variables from elastic net in traditional cox model and plot hazard ratios
    coef.min = coef(opt.model.min)
    coef.min.names = coef.min@Dimnames[[1]]
    coef.min.names = coef.min.names[which(coef(opt.model.min) != 0)]
    opt.model.min = survival::coxph(formula = as.formula(paste0('Surv(time, status)', ' ~ ',
                                                                paste0(coef.min.names, collapse = '+'),
                                                                sep = '')),
                                    data = combined_data_sub,
                                    x = T,
                                    y = T)
    
    plot.min = ggforest(opt.model.min, 
                        data = combined_data_sub)
    
    
    coef.1se = coef(opt.model.1se)
    coef.1se.names = coef.1se@Dimnames[[1]]
    coef.1se.names = coef.1se.names[which(coef(opt.model.1se) != 0)]
    opt.model.1se = survival::coxph(formula = as.formula(paste0('Surv(time, status)', ' ~ ',
                                                                paste0(coef.1se.names, collapse = '+'),
                                                                sep = '')),
                                    data = combined_data_sub,
                                    x = T,
                                    y = T)
    
    plot.1se = ggforest(opt.model.1se,
                        data = combined_data_sub)
    
    ## Make json of optimal params and write to file for use in 05-model_compare.R
    optimal.params = sprintf('{"lambda.min": "%s", 
                          "alpha.min": "%s", 
                          "lambda.1se": "%s",
                          "alpha.1se": %s}', 
                          opt.param.lamda_min$lamda.min[ind.lambda.min], 
                          alpha.min, 
                          opt.param.lamda1se$lamda.1se[ind.lambda1se],
                          alpha.1se)
    
    write(optimal.params, file.path(results_dir, "coxlasso_params.json"))

  }
}