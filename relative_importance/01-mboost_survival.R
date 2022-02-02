# Author: Run Jin
# Perform mboost survival analysis to measure the relative importance of SLC1A5 
# relative to known amino acid transporters with respect to disease prognosis
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("mboost")
  library("survival")
  library("EGSEA")
  library("org.Hs.eg.db")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-e","--expression"),type="character",
              help="gene expression rsem tpm file for OpenPedCan RNA samples (.rds) "),
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"),
  make_option(c("-m","--short_long_match"),type="character",
              help="match between long and short names (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cg_list <-unlist(strsplit(opt$cg_interest,","))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "relative_importance")
results_dir_sum <- file.path(analysis_dir, "results", "mboost", "model_summary")
if(!dir.exists(results_dir_sum)){
  dir.create(results_dir_sum, recursive=TRUE)
}

results_dir_pred <- file.path(analysis_dir, "results", "mboost", "predicted_scores")
if(!dir.exists(results_dir_pred)){
  dir.create(results_dir_pred, recursive=TRUE)
}

results_dir_fit <- file.path(analysis_dir, "results", "mboost", "model_fit")
if(!dir.exists(results_dir_fit)){
  dir.create(results_dir_fit, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "mboost")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

#### Read in files necessary for analyses -----------------------------------
histology_df <- readr::read_tsv(opt$histology, guess_max=100000)
expression_data <- readRDS(opt$expression)
short_long_match <- readr::read_tsv(opt$short_long_match)

#### Calculate PFS status based on the PFS days and OS days -----------------------------------
histology_df$PFS_days <- as.numeric(histology_df$PFS_days)

histology_df <- histology_df %>%
  # keep only initial CNV tumor
  dplyr::filter(sample_type == "Tumor") %>% 
  dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
  # does not include TCGA or GTEx
  dplyr::filter(cohort %in% c("PBTA","GMKF","TARGET")) %>% 
  dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
  # exclude derived cell line
  dplyr::filter(composition != "Derived Cell Line") %>% 
  # keep only unique Kids First Participant 
  dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>% 
  dplyr::filter(!is.na(OS_status)) %>%
  dplyr::mutate(os_status_level = case_when(
    OS_status == "LIVING" ~ 0,
    OS_status == "DECEASED" ~ 1)) %>%
  dplyr::mutate(PFS_status = if_else(PFS_days < OS_days, 1, 0)) 

####### filter gene expression to contain genes of interest 
## Entrez IDs for SLC1A5 and ACLY
entrezs = c('6510', '47')
reac.annots = buildMSigDBIdx(entrezIDs = entrezs, species = "human", geneSets = "c2")

# find the gene ids that are associated with the entrez ID in pathway of interest 
entrez_symbol_of_interest <- reac.annots[["c2"]]@original[["REACTOME_AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE"]] %>%
  as.data.frame() %>% 
  dplyr::mutate(gene_symbol = mapIds(org.Hs.eg.db, ., "SYMBOL", "ENTREZID")) 

# filter to contain only genes of interest
expression_data <- expression_data %>%
  tibble::rownames_to_column("gene_symbol") %>% 
  dplyr::filter(gene_symbol %in% entrez_symbol_of_interest$gene_symbol) %>% 
  tibble::column_to_rownames("gene_symbol") %>% 
  # also filter to contain only samples in all cancer groups of interest 
  dplyr::select(histology_df$Kids_First_Biospecimen_ID)

####### do the analysis on all the cancer group of interest ----------------------------------
for (i in 1:length(cg_list)){
  # find the cancer group of interest 
  x <- cg_list[i]
  
  # match the long name to the short name
  long_name <- short_long_match %>% filter(short_name == x) %>%
    pull(long_name) 
  
  # filter to the cohort of interest
  cohort_df <- histology_df %>% dplyr::filter(cancer_group %in% long_name) %>%
    dplyr::select(Kids_First_Biospecimen_ID, os_status_level, OS_days, PFS_status, PFS_days) 
  
  # get biospecimen ID's for samples 
  cohort_bsid <- cohort_df %>% pull(Kids_First_Biospecimen_ID) %>% unique()
  
  # find the tpm of target genes for these samples 
  expression_data_of_interest <- expression_data %>%
    dplyr::select(all_of(cohort_bsid)) %>% 
    t() %>%
    as.data.frame() %>% 
    tibble::rownames_to_column("Kids_First_Biospecimen_ID")
  
  # combine histology info with expression 
  combined_data <- cohort_df %>% 
    left_join(expression_data_of_interest) %>%
    tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  
  # generate gene variables
  gene_variables <- rownames(expression_data) %>% paste(collapse="+")
  
  if(x == "LGG"){
    # generate T/F vector for PFS status 
    censored_pfs <- combined_data$PFS_status == "1"
    
    ########## perform glmboost --------------------
    glm_fit_coxph <- glmboost(as.formula(paste0("Surv(combined_data$PFS_days, censored_pfs) ~ ", gene_variables)),
                              data = combined_data, 
                              family = CoxPH(), 
                              center = TRUE,
                              control =  boost_control(mstop = 500))
    # determine cross-validated optimal mstop parameter
    cvm.coxph = cvrisk(glm_fit_coxph)
    invisible(glm_fit_coxph[mstop(cvm.coxph)])
    
    # save the fit results
    glm_fit_coxph %>% 
      saveRDS(file.path(results_dir_fit,  paste0("coxph_pfs_model_fit_in_", x, ".rds")))
    
    # write out the predicted results for later analysis 
    pred_result <- glm_fit_coxph$predict() %>% as.data.frame() 
    colnames(pred_result) <- "predicted_score"
    pred_result %>% 
      tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
      readr::write_tsv(file.path(results_dir_pred, paste0("coxph_pfs_pred_risk_in_", x, ".tsv")))
      
    # write out the formula summary
    sink(file = file.path(results_dir_sum, paste0("coxph_pfs_summary_in_", x, ".txt")))
    print(summary(glm_fit_coxph))
    print(paste0('Log-Likelihood: ', glm_fit_coxph$logLik()))
    sink()
    
    # plot out the results
    pdf(file.path(plots_dir, paste0("mstop_diagnostics", x, ".pdf")))
    p <- plot(cvm.coxph)
    print(p)
    dev.off()
    
    # perform glmboost with Loglog
    glm_fit_loglog <- glmboost(as.formula(paste0("Surv(combined_data$PFS_days, censored_pfs) ~ ", gene_variables)),
                               data = combined_data, 
                               family = Loglog(), 
                               center = TRUE,
                               control =  boost_control(mstop = 500))
    # determine cross-validated optimal mstop parameter
    cvm.log = cvrisk(glm_fit_loglog)
    invisible(glm_fit_loglog[mstop(cvm.log)])
    
    # save the fit results
    glm_fit_loglog %>% 
      saveRDS(file.path(results_dir_fit,  paste0("loglog_pfs_model_fit_in_", x, ".rds")))
    
    # write out the predicted results for later analysis 
    pred_result <- glm_fit_loglog$predict() %>% as.data.frame() 
    colnames(pred_result) <- "predicted_score"
    pred_result %>% 
      tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
      readr::write_tsv(file.path(results_dir_pred, paste0("loglog_pfs_pred_risk_in_", x, ".tsv")))
    
    sink(file = file.path(results_dir_sum, paste0("loglog_pfs_summary_in_", x, ".txt")))
    print(summary(glm_fit_loglog))
    print(paste0('Log-Likelihood: ', glm_fit_loglog$logLik()))
    sink()
    
    # plot out the results
    pdf(file.path(plots_dir, paste0("mstop_diagnostics", x, ".pdf")))
    p <- plot(cvm.log)
    print(p)
    dev.off()
    
    ########## perform glmboost --------------------
    # # calculate IPC weights
    # iw <- IPCweights(Surv(combined_data$PFS_days, censored_pfs))
    # 
    # # filter columns that are not used in the model fitting
    # filtered <- combined_data %>% 
    #   dplyr::select(-c("PFS_status", "os_status_level", "OS_days"))
    # 
    # # fit a weighted linear model by boosting with componentwise linear weighted least squares 
    # filtered_surv <- glmboost(log(PFS_days) ~ ., 
    #                           data = filtered,
    #                           weights = iw, 
    #                           center = TRUE, 
    #                           control = boost_control(mstop = 500))
    # # find the optimal mstop 
    # mstop(aic <- AIC(filtered_surv))
    # filtered_surv <- filtered_surv[mstop(aic)]
    # 
    # # run prediction 
    # predicted <- predict(filtered_surv) %>% 
    #   as.data.frame() %>%
    #   dplyr::rename(predict_log_days = V1) %>%
    #   tibble::rownames_to_column("Kids_First_Biospecimen_ID") 
    # 
    # # combine the weight 
    # names(iw) <- rownames(combined_data)
    # weight_df <- as.data.frame(iw) %>% 
    #   tibble::rownames_to_column("Kids_First_Biospecimen_ID") 
    # colnames(weight_df) <- c("Kids_First_Biospecimen_ID", "IPC_weight") 
    # 
    # # compare with original data
    # combined_predicted_ori <- combined_data %>%
    #   tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>% 
    #   dplyr::select(Kids_First_Biospecimen_ID, PFS_days) %>% 
    #   left_join(weight_df) %>%
    #   left_join(predicted) %>%
    #   dplyr::mutate(ori_log_days = log(PFS_days))
    # 
    # # plot out the results
    # pdf(file.path(plots_dir, paste0("correlation_plots_by_pfs_in_", x, ".pdf")))
    # p <- ggplot(combined_predicted_ori, aes(x=ori_log_days, y=predict_log_days)) +
    #       geom_point(aes(size = IPC_weight))
    # print(p)
    # dev.off()
  }
  
  if(x != "LGG"){
    # generate T/F vector for OS status 
    censored_os <- combined_data$os_status_level == "1"
    
    ########## perform glmboost --------------------
    glm_fit_coxph <- glmboost(as.formula(paste0("Surv(combined_data$OS_days, censored_os) ~ ", gene_variables)),
                              data = combined_data, 
                              family = CoxPH(), 
                              center = TRUE,
                              control =  boost_control(mstop = 1000))
    # determine cross-validated optimal mstop parameter
    cvm.coxph = cvrisk(glm_fit_coxph)
    invisible(glm_fit_coxph[mstop(cvm.coxph)])
    
    # save the fit results
    glm_fit_coxph%>% 
      saveRDS(file.path(results_dir_fit,  paste0("coxph_os_model_fit_in_", x, ".rds")))
    
    # write out the predicted results for later analysis 
    pred_result <- glm_fit_coxph$predict() %>% as.data.frame() 
    colnames(pred_result) <- "predicted_score"
    rownames(pred_result) <- glm_fit_coxph[["rownames"]]
    pred_result %>% 
      tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
      readr::write_tsv(file.path(results_dir_pred, paste0("coxph_os_pred_risk_in_", x, ".tsv")))
    
    sink(file = file.path(results_dir_sum, paste0("coxph_os_summary_in_", x, ".txt")))
    print(summary(glm_fit_coxph))
    print(paste0('Log-Likelihood: ', glm_fit_coxph$logLik()))
    sink()
    
    # plot out the results
    pdf(file.path(plots_dir, paste0("mstop_diagnostics", x, ".pdf")))
    p <- plot(cvm.coxph)
    print(p)
    dev.off()
    
    # perform glmboost with Loglog
    glm_fit_loglog <- glmboost(as.formula(paste0("Surv(combined_data$OS_days, censored_os) ~ ", gene_variables)),
                               data = combined_data, 
                               family = Loglog(), 
                               center = TRUE,
                               control =  boost_control(mstop = 500))
    
    
    # determine cross-validated optimal mstop parameter
    cvm.log = cvrisk(glm_fit_loglog)
    invisible(glm_fit_loglog[mstop(cvm.log)])
    
    # save the fit results
    glm_fit_loglog %>% 
      saveRDS(file.path(results_dir_fit,  paste0("loglog_os_model_fit_in_", x, ".rds")))
    
    # write out the predicted results for later analysis 
    pred_result <- glm_fit_loglog$predict() %>% as.data.frame() 
    colnames(pred_result) <- "predicted_score"
    rownames(pred_result) <- glm_fit_loglog[["rownames"]]
    pred_result %>% 
      tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
      readr::write_tsv(file.path(results_dir_pred, paste0("loglog_os_pred_risk_in_", x, ".tsv")))
    
    sink(file = file.path(results_dir_sum, paste0("loglog_os_summary_in_", x, ".txt")))
    print(summary(glm_fit_loglog))
    print(paste0('Log-Likelihood: ', glm_fit_loglog$logLik()))
    sink()
    
    # plot out the results
    pdf(file.path(plots_dir, paste0("mstop_diagnostics", x, ".pdf")))
    p <- plot(cvm.log)
    print(p)
    dev.off()
    
    # ########## perform glmboost --------------------
    # # calculate IPC weights
    # iw <- IPCweights(Surv(combined_data$OS_days, censored_os))
    # 
    # # filter columns that are not used in the model fitting
    # filtered <- combined_data %>% 
    #   dplyr::select(-c("PFS_status", "os_status_level", "PFS_days"))
    # 
    # # fit a weighted linear model by boosting with componentwise linear weighted least squares 
    # filtered_surv <- glmboost(log(OS_days) ~ ., 
    #                           data = filtered,
    #                           weights = iw, 
    #                           center = TRUE, 
    #                           control = boost_control(mstop = 500))
    # # find the optimal mstop 
    # mstop(aic <- AIC(filtered_surv))
    # filtered_surv <- filtered_surv[mstop(aic)]
    # 
    # # run prediction 
    # predicted <- predict(filtered_surv) %>% 
    #   as.data.frame() %>%
    #   dplyr::rename(predict_log_days = V1) %>%
    #   tibble::rownames_to_column("Kids_First_Biospecimen_ID") 
    # 
    # # combine the weight 
    # names(iw) <- rownames(combined_data)
    # weight_df <- as.data.frame(iw) %>% 
    #   tibble::rownames_to_column("Kids_First_Biospecimen_ID") 
    # colnames(weight_df) <- c("Kids_First_Biospecimen_ID", "IPC_weight") 
    # 
    # # compare with original data
    # combined_predicted_ori <- combined_data %>%
    #   tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>% 
    #   dplyr::select(Kids_First_Biospecimen_ID, OS_days) %>% 
    #   left_join(weight_df) %>%
    #   left_join(predicted) %>%
    #   dplyr::mutate(ori_log_days = log(OS_days))
    
    # # plot out the results
    # pdf(file.path(plots_dir, paste0("correlation_plots_by_os_in_", x, ".pdf")))
    # p <- ggplot(combined_predicted_ori, aes(x=ori_log_days, y=predict_log_days)) +
    #   geom_point(aes(size = IPC_weight))
    # print(p)
    # dev.off()
  }
}

