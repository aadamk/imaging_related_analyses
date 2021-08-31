# Author: Run Jin
#
# Obtain Kaplan-Meier Survival Statistic Results and Plots

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("survminer"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("survival"))

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
stat_list <- lapply(cancer_group_list, function(x){
  
  # match the long name to the short name
  long_name <- short_long_match %>% filter(short_name == x) %>%
    pull(long_name) 
  # make directory for results and plots
  km_survival_plots_dir <- file.path(root_dir, "molecular_target_radiation", "plots",x, "KM_survival")
  
  # filter to the cohort of interest
  cohort_df <- histology_df %>% dplyr::filter(cancer_group %in% long_name) %>%
    dplyr::select(Kids_First_Biospecimen_ID, os_status_level, OS_days, PFS_status,PFS_days) 
  
  # get biospecimen ID's for samples 
  cohort_bsid <- cohort_df %>% pull(Kids_First_Biospecimen_ID) %>% unique()
  
  # subset to gene and sample of interest          
  expression_of_interest <- expression_data %>% dplyr::select(all_of(cohort_bsid)) %>% 
    tibble::rownames_to_column("GeneSymbol") %>% 
    dplyr::filter(GeneSymbol %in% gene_list) %>% 
    tibble::column_to_rownames("GeneSymbol") %>%
    t() %>% as.data.frame() %>%
    tibble::rownames_to_column("Kids_First_Biospecimen_ID")
  
  #annotate relevant clinical information to the expression data frame
  combined <- expression_of_interest %>% 
    dplyr::left_join(cohort_df) 
  
  #for genes of interest, generate grouping based on that gene
  lapply(gene_list,function(y){
    # generate annotated expression matrix for the modeling below
    median_exp <- median(combined[[y]]) %>% as.numeric()

    # assign high or low expression to samples
    combined_annotated <- combined %>% 
      rename("gene_of_interest" = y) %>%
      mutate(expression_level = case_when(
        gene_of_interest>median_exp ~"High",
        gene_of_interest<=median_exp ~ "low"
    ))
    
    ########################################## survival analysis OS
    # generate the log-rank model 
    fit_gene <- survival::survdiff(
      as.formula("survival::Surv(OS_days, os_status_level) ~ expression_level"),
      data = combined_annotated
    )
    # get the p.values and all statistics for the model 
    fit_gene$p.value <- round(pchisq(fit_gene$chisq, df = 1, lower = FALSE), digits=3)
    fit_gene_df <- as.data.frame(fit_gene[c("n", "obs", "exp", "chisq", "p.value")])
    
    fit_gene_df <- fit_gene_df %>% mutate(cancer_group = x,
                           gene_interest = y, 
                           model = "OS")

    # generate the kaplan-meier model 
    kap_fit_gene <- survival::survfit(
      as.formula("survival::Surv(OS_days, os_status_level) ~ expression_level"),
      data = combined_annotated
    )
    
    # generate and save the survival plot 
    surv_plot <- survminer::ggsurvplot(kap_fit_gene,
                                       pval = round(fit_gene$p.value, digits=3), # use the computed pval from log rank analyses 
                                       data = combined_annotated,
                                       risk.table = TRUE,
                                       xlim = c(0, 2000),
                                       break.time.by = 500,
                                       ggtheme = theme_minimal(),
                                       risk.table.y.text.col = TRUE,
                                       risk.table.y.text = FALSE
    )
    # Make this plot a combined plot
    surv_plot <- cowplot::plot_grid(surv_plot[[1]], surv_plot[[2]], nrow = 2, 
                                    rel_heights = c(2.5, 1))
    
    km_survival_plot_sub_dir <- file.path(km_survival_plots_dir, "OS",y)
    if (!dir.exists(km_survival_plot_sub_dir )) {
      dir.create(km_survival_plot_sub_dir , recursive = TRUE)
    }
    cowplot::save_plot(filename = file.path(km_survival_plot_sub_dir,"kaplan_meier_survival.png"), plot = surv_plot)
    
    
    ########################################## survival analysis PFS
    # generate the log-rank model 
    fit_gene_pfs <- survival::survdiff(
      as.formula("survival::Surv(PFS_days, PFS_status) ~ expression_level"),
      data = combined_annotated
    )
    # get the p.values and all statistics for the model 
    fit_gene_pfs$p.value <- round(pchisq(fit_gene_pfs$chisq, df = 1, lower = FALSE), digits=3)
    fit_gene_pfs_df <- as.data.frame(fit_gene_pfs[c("n", "obs", "exp", "chisq", "p.value")])
    
    fit_gene_pfs_df <- fit_gene_pfs_df %>% mutate(cancer_group = x,
                                          gene_interest = y, 
                                          model = "PFS")
    
    # generate the kaplan-meier model 
    kap_fit_gene <- survival::survfit(
      as.formula("survival::Surv(PFS_days, PFS_status) ~ expression_level"),
      data = combined_annotated
    )
    
    # generate and save the survival plot 
    surv_plot <- survminer::ggsurvplot(kap_fit_gene,
                                       pval = round(fit_gene_pfs$p.value,digits=3), # use the computed pval from log rank analyses 
                                       data = combined_annotated,
                                       risk.table = TRUE,
                                       xlim = c(0, 2000),
                                       break.time.by = 500,
                                       ggtheme = theme_minimal(),
                                       risk.table.y.text.col = TRUE,
                                       risk.table.y.text = FALSE
    )
    # Make this plot a combined plot
    surv_plot <- cowplot::plot_grid(surv_plot[[1]], surv_plot[[2]], nrow = 2, 
                                    rel_heights = c(2.5, 1))
    
    km_survival_plot_pfs_dir <- file.path(km_survival_plots_dir, "PFS",y)
    if (!dir.exists(km_survival_plot_pfs_dir )) {
      dir.create(km_survival_plot_pfs_dir , recursive = TRUE)
    }

    cowplot::save_plot(filename = file.path(km_survival_plot_pfs_dir,"kaplan_meier_survival.png"), plot = surv_plot)
    
    combined_df <- rbind(fit_gene_df, fit_gene_pfs_df)
    return(combined_df)
  })
})

cg_combined <- lapply(stat_list, function(x){
  do.call(rbind,x)
})
all_combined <- do.call(rbind, cg_combined)
readr::write_tsv(all_combined, opt$stat_outfile)

  