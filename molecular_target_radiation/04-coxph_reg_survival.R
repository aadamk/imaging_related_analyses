# Author: Run Jin
#
# Obtain CoxPH Survival Statistic Results and Plots
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("survival"))
suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("ggpubr"))

#### Parse command line options ------------------------------------------------

option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-e","--expression"),type="character",
              help="gene expression rsem tpm file for OpenPedCan RNA samples (.rds) "),
  make_option(c("-l","--cancer_groups"),type="character",
              help="list of cancer groups that are of interest"),
  make_option(c("-m","--short_long_match"),type="character",
              help="match between long and short names"),
  make_option(c("-g","--gene_list"),type="character",
              help="list of genes that we are interested in analyzing"),
  make_option(c("-s","--stat_outfile"),type="character",
              help="path and file name for the stat output (.TSV)")
)

opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cancer_group_list <-unlist(strsplit(opt$cancer_groups,","))
gene_list <-unlist(strsplit(opt$gene_list,","))

#### Define Directories ----------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "molecular_target_radiation")

results_dir <- file.path(analysis_dir, "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir , recursive = TRUE)
}

#### Read in files necessary for analyses -----------------------------------
histology_df <- readr::read_tsv(opt$histology)
expression_data <- readRDS(opt$expression)
short_long_match <- readr::read_tsv(opt$short_long_match)

#### Calculate PFS status based on the PFS days and OS days -----------------------------------
histology_df$PFS_days <- as.numeric(histology_df$PFS_days)
histology_df <- histology_df %>%
  dplyr::mutate(PFS_status = if_else(PFS_days < OS_days, 1, 0)) %>% 
  dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
  dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
  dplyr::filter(!is.na(OS_status)) %>%
  dplyr::mutate(os_status_level = case_when(
    OS_status == "LIVING" ~ 0,
    OS_status == "DECEASED" ~ 1
  )) 

#### Do the analysis for all the cohort of interest -----------------------------------

train_stats <- data.frame()
test_stats <- data.frame()

for (i in 1:length(cancer_group_list)){
  
  x <- cancer_group_list[i]
  
  # define directory
  cox_survival_plots_dir <- file.path(root_dir, "molecular_target_radiation", "plots",x, "CoxPH_survival")
  
  if (!dir.exists(cox_survival_plots_dir )) {
    dir.create(cox_survival_plots_dir , recursive = TRUE)
  }
  
  # match the long name to the short name
  long_name <- short_long_match %>% filter(short_name == x) %>%
    pull(long_name) 
  
  # filter to the cohort of interest
  cohort_df <- histology_df %>% dplyr::filter(cancer_group %in% long_name) %>%
    dplyr::select(Kids_First_Biospecimen_ID, os_status_level, OS_days, PFS_status,PFS_days, CNS_region, harmonized_diagnosis) 
  
  # get biospecimen ID's for samples 
  cohort_bsid <- cohort_df %>% pull(Kids_First_Biospecimen_ID) %>% unique()
  
  for (j in 1:length(gene_list)){
    y <- gene_list[j]
    # subset to gene and sample of interest          
    expression_of_interest <- expression_data %>% dplyr::select(all_of(cohort_bsid)) %>% 
      tibble::rownames_to_column("GeneSymbol") %>% 
      dplyr::filter(GeneSymbol ==y) %>% 
      tibble::column_to_rownames("GeneSymbol") %>%
      t() %>% as.data.frame() %>%
      tibble::rownames_to_column("Kids_First_Biospecimen_ID")
    
    #annotate relevant clinical information to the expression data frame
    combined_annotated <- expression_of_interest %>% dplyr::left_join(cohort_df) 
    
    ########################################## survival analysis OS
    # separate 
    set.seed(236)
    
    # change the columns to desired type 
    combined_annotated$harmonized_diagnosis <- as.factor(combined_annotated$harmonized_diagnosis)
    combined_annotated$CNS_region <- as.factor(combined_annotated$CNS_region)
    
    trainIndex <- createDataPartition(combined_annotated$os_status_level, p = 0.6, 
                                      list = FALSE, 
                                      times = 1)
    data_train <- combined_annotated[trainIndex, ]
    data_test <- combined_annotated[-trainIndex, ]
    
    # Define models
    model_gene <- as.formula(paste0("survival::Surv(OS_days, os_status_level) ~ harmonized_diagnosis + CNS_region + ",y))
    
    coxph_train_gene <- survival::coxph(model_gene,
                                        data = data_train)
    # generate the stats summary
    coxph_train_gene_sum_os <- as.data.frame(summary(coxph_train_gene)$coefficients, keep.rownames = T) %>% 
      tibble::rownames_to_column() %>% 
      mutate(cancer_group =x,
             gene_list = y,
             survival = "OS")

    # define between risk groups
    train_risk_gene <- predict(coxph_train_gene, type = "risk")
    med_risk_gene<- median(train_risk_gene)
    # assign the groups between risk groups
    test_risk_gene <- predict(coxph_train_gene, newdata = data_test, type = "risk") %>% 
      as.data.frame() %>% cbind(data_test) %>%
      dplyr::rename("RiskScore" = ".") %>% 
      dplyr::mutate(RiskGroup = case_when(
        RiskScore > med_risk_gene ~ "High", 
        RiskScore <= med_risk_gene ~ "Low"
      )) %>%
      dplyr::rename(gene_of_interest = all_of(y))
    
    #calculate n for each risk group
    test_risk_n <- test_risk_gene %>% 
      dplyr::group_by(RiskGroup) %>%
      dplyr::summarize(n=n()) %>%
      dplyr::mutate(risk_group_n = paste0(RiskGroup, " n=", n))
    
    #add n to the risk groups
    test_risk_gene <- test_risk_gene %>% 
      dplyr::left_join(test_risk_n)
    
    # generate boxplots showing the expression of high vs. low risk group
    exp_plot <- test_risk_gene %>%
      ggplot( aes(x=risk_group_n, y=gene_of_interest)) +
      geom_violin(width=1, trim=TRUE, show.legend = F, aes(fill=risk_group_n)) +
      geom_boxplot(width=0.1, color="black", show.legend = F,aes(fill=risk_group_n)) +
      labs(title=paste0(y," Expression of Risk Groups"),x="Risk Group", y = paste0(y," TPM Value")) +
      theme(axis.text.x = element_text(size = 18),
            axis.title=element_text(size=20,face="bold"))  +
      stat_compare_means(size = 6)
    
    # make risk groups into factors
    test_risk_gene$RiskGroup <- as.factor(test_risk_gene$RiskGroup)
    
    # calculate stats for risk group per gene
    coxph_riskgroup_gene <- survival::coxph(formula = Surv(OS_days, os_status_level) ~  RiskGroup,
                                            data = test_risk_gene)
    coxph_riskgroup_gene_sum <- as.data.frame(summary(coxph_riskgroup_gene)$coefficients) %>% 
      tibble::rownames_to_column() 
    
    # calculate stats for risk score per gene
    coxph_riskscore_gene <- survival::coxph(formula = Surv(OS_days, os_status_level) ~  RiskScore,
                                            data = test_risk_gene)
    coxph_riskscore_gene_sum <- as.data.frame(summary(coxph_riskscore_gene)$coefficients) %>% 
      tibble::rownames_to_column() 
    
    combined_coxph_sum_os <- rbind(coxph_riskgroup_gene_sum, coxph_riskscore_gene_sum) %>%
      mutate(cancer_group =x,
             gene_list = y,
             survival = "OS")

    ########################################## plot generation for OS survival
    fit_gene <- survfit(Surv(OS_days, os_status_level) ~ RiskGroup, data = test_risk_gene)
    plot_coxph_riskgroup_gene <- survminer::ggsurvplot(fit_gene,
                                                       data=test_risk_gene,                                
                                                       pval = TRUE, 
                                                       conf.int = TRUE,
                                                       risk.table = TRUE, # Add risk table
                                                       risk.table.col = "strata", # Change risk table color by groups
                                                       linetype = "strata", # Change line type by groups
                                                       surv.median.line = "hv", # Specify median survival
                                                       ggtheme = theme_bw(), # Change ggplot2 theme
                                                       palette = c("#E7B800", "#2E9FDF"))
    
    # Make this plot a combined plot
    surv_plot_gene <- cowplot::plot_grid(plot_coxph_riskgroup_gene[[1]], 
                                         plot_coxph_riskgroup_gene[[2]], 
                                         nrow = 2, 
                                         rel_heights = c(2.5, 1))
    
    cox_survival_plots_os_dir <- file.path(cox_survival_plots_dir, "OS",y)
    if (!dir.exists(cox_survival_plots_os_dir )) {
      dir.create(cox_survival_plots_os_dir , recursive = TRUE)
    }
    
    # Save the plot
    cowplot::save_plot(filename = file.path(cox_survival_plots_os_dir, "coxph_riskgroup_survival.png"), plot = surv_plot_gene)
    # Save expression plot
    ggsave(filename = file.path(cox_survival_plots_os_dir, "coxph_riskgroup_tpm.png"), plot = exp_plot, height = 4, width=6)
    
    ########################################## survival analysis PFS
    
    combined_annotated <- combined_annotated %>% 
      dplyr::filter(!is.na(PFS_status))
    # separate 
    set.seed(236)
    trainIndex <- createDataPartition(combined_annotated$PFS_status, p = 0.6, 
                                      list = FALSE, 
                                      times = 1)
    data_train <- combined_annotated[trainIndex, ]
    data_test <- combined_annotated[-trainIndex, ]
    
    model_gene <- as.formula(paste0("survival::Surv(PFS_days, PFS_status) ~ harmonized_diagnosis + CNS_region + ",y))
    coxph_train_gene <- survival::coxph(model_gene,
                                        data = data_train)
    # generate the stats summary
    coxph_train_gene_sum_pfs <- as.data.frame(summary(coxph_train_gene)$coefficients, keep.rownames = T) %>% 
      tibble::rownames_to_column() %>% 
      mutate(cancer_group =x,
             gene_list = y,
             survival = "PFS")
    
    # define between risk groups
    train_risk_gene <- predict(coxph_train_gene, type = "risk")
    med_risk_gene<- median(train_risk_gene)
    
    # assign the groups between risk groups
    test_risk_gene <- predict(coxph_train_gene, newdata = data_test, type = "risk") %>% 
      as.data.frame() %>% cbind(data_test) %>%
      dplyr::rename("RiskScore" = ".") %>% 
      dplyr::mutate(RiskGroup = case_when(
        RiskScore > med_risk_gene ~ "High", 
        RiskScore <= med_risk_gene ~ "Low"
      )) %>%
      dplyr::rename(gene_of_interest = all_of(y)) 
    
    #calculate n for each risk group
    test_risk_n <- test_risk_gene %>% 
      dplyr::group_by(RiskGroup) %>%
      dplyr::summarize(n=n()) %>%
      dplyr::mutate(risk_group_n = paste0(RiskGroup, " n=", n))
    
    #add n to the risk groups
    test_risk_gene <- test_risk_gene %>% 
      dplyr::left_join(test_risk_n)
    
    # generate boxplots showing the expression of high vs. low risk group
    exp_plot <- test_risk_gene %>%
      ggplot( aes(x=risk_group_n, y=gene_of_interest)) +
      geom_violin(width=1, trim=TRUE, show.legend = F, aes(fill=risk_group_n)) +
      geom_boxplot(width=0.1, color="black", show.legend = F,aes(fill=risk_group_n)) +
      labs(title=paste0(y," Expression of Risk Groups"),x="Risk Group", y = paste0(y," TPM Value")) +
      theme(axis.text.x = element_text(size = 18),
            axis.title=element_text(size=20,face="bold"))  +
      stat_compare_means(size = 6)
    
    # make risk groups into factors
    test_risk_gene$RiskGroup <- as.factor(test_risk_gene$RiskGroup)
    
    # calculate stats for risk group per gene
    coxph_riskgroup_gene <- survival::coxph(formula = Surv(PFS_days, PFS_status) ~  RiskGroup,
                                            data = test_risk_gene)
    coxph_riskgroup_gene_sum <- as.data.frame(summary(coxph_riskgroup_gene)$coefficients) %>% 
      tibble::rownames_to_column() 
    
    # calculate stats for risk score per gene
    coxph_riskscore_gene <- survival::coxph(formula = Surv(PFS_days, PFS_status) ~  RiskScore,
                                            data = test_risk_gene)
    coxph_riskscore_gene_sum <- as.data.frame(summary(coxph_riskscore_gene)$coefficients) %>% 
      tibble::rownames_to_column() 
    
    combined_coxph_sum_pfs <- rbind(coxph_riskgroup_gene_sum, coxph_riskscore_gene_sum) %>%
      mutate(cancer_group =x,
             gene_list = y,
             survival = "PFS")
    
    combined_coxph_sum <- rbind(combined_coxph_sum_os, combined_coxph_sum_pfs)
    test_stats <- rbind(test_stats, combined_coxph_sum)
    
    combined_coxph_train <- rbind(coxph_train_gene_sum_os, coxph_train_gene_sum_pfs)   
    train_stats <- rbind(train_stats, combined_coxph_train)
    
   
    ########################################## plot generation for PFS survival
    fit_gene_pfs <- survfit(Surv(PFS_days, PFS_status) ~ RiskGroup, data = test_risk_gene)
    plot_coxph_riskgroup_gene_pfs <- survminer::ggsurvplot(fit_gene_pfs,
                                                           data=test_risk_gene,                                
                                                           pval = TRUE, 
                                                           conf.int = TRUE,
                                                           risk.table = TRUE, # Add risk table
                                                           risk.table.col = "strata", # Change risk table color by groups
                                                           linetype = "strata", # Change line type by groups
                                                           surv.median.line = "hv", # Specify median survival
                                                           ggtheme = theme_bw(), # Change ggplot2 theme
                                                           palette = c("#E7B800", "#2E9FDF"))
    
    # Make this plot a combined plot
    surv_plot_gene_pfs <- cowplot::plot_grid(plot_coxph_riskgroup_gene_pfs[[1]], 
                                             plot_coxph_riskgroup_gene_pfs[[2]], 
                                             nrow = 2, 
                                             rel_heights = c(2.5, 1))
    
    cox_survival_plots_pfs_dir <- file.path(cox_survival_plots_dir, "PFS",y)
    if (!dir.exists(cox_survival_plots_pfs_dir )) {
      dir.create(cox_survival_plots_pfs_dir , recursive = TRUE)
    }
    
    # Save the plot
    cowplot::save_plot(filename = file.path(cox_survival_plots_pfs_dir, "coxph_riskgroup_survival.png"), plot = surv_plot_gene_pfs)
    # Save expression plot
    ggsave(filename = file.path(cox_survival_plots_pfs_dir, "coxph_riskgroup_tpm.png"), plot = exp_plot, height = 4, width=6)

  }
}
                  
readr::write_tsv(train_stats, file.path(results_dir,"coxph_train_stats_combined.tsv"))
readr::write_tsv(test_stats, file.path(results_dir,"coxph_test_stats_combined.tsv"))


