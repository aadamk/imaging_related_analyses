# Author: Run Jin
#
# Obtain TPM plots for genes of interest with cancer group of interest

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("car"))
suppressPackageStartupMessages(library("multcomp"))
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
  make_option(c("-s","--stat_outpath"),type="character",
              help="path and file name for the stat output (.TSV)")
  
)

opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cancer_group_list <-unlist(strsplit(opt$cancer_groups,","))
gene_list <-unlist(strsplit(opt$gene_list,","))

#### Define Directories ----------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
stats_anova_dir <- file.path(root_dir, "molecular_target_radiation", "stats", "anova")
stats_tukey_dir <- file.path(root_dir, "molecular_target_radiation", "stats", "tukey")

tpm_plots_dir <- file.path(root_dir, "molecular_target_radiation", "plots", "TPM_violin")

if (!dir.exists(stats_anova_dir )) {
  dir.create(stats_anova_dir , recursive = TRUE)
}

if (!dir.exists(stats_tukey_dir )) {
  dir.create(stats_tukey_dir , recursive = TRUE)
}

if (!dir.exists(tpm_plots_dir )) {
  dir.create(tpm_plots_dir , recursive = TRUE)
}

#### Read in files necessary for analyses -----------------------------------
histology_df <- readr::read_tsv(opt$histology)
expression_data <- readRDS(opt$expression)
short_long_match <- readr::read_tsv(opt$short_long_match)

#### Filter to the cancer group of interest -----------------------------------
cohort_df_list <- lapply(cancer_group_list, function(x){
  cancer_group_long <- short_long_match %>% 
    filter(short_name ==x) %>%
    pull(long_name)
    
  cohort_df_each <- histology_df %>% 
    # keep only initial CNV tumor
    dplyr::filter(sample_type == "Tumor") %>% 
    dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
    # does not include TCGA or GTEx
    dplyr::filter(cohort %in% c("PBTA","GMKF","TARGET")) %>% 
    dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
    dplyr::filter(cancer_group %in% cancer_group_long) %>%
    # exclude derived cell line
    dplyr::filter(composition != "Derived Cell Line") %>% 
    # keep only unique Kids First Participant 
    dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>% 
    dplyr::select(Kids_First_Biospecimen_ID, harmonized_diagnosis) %>%
    dplyr::mutate(short_name = x)
})

cohort_df <- do.call(rbind,cohort_df_list)

#### Generate figures and statistics for each gene ----------------------------------
anova_results <- data.frame()
tukey_results <- data.frame()

for(i in 1:length(gene_list)){
  # filter to gene of interest
  x<- gene_list[i]
  cohort_bsid <- cohort_df %>% pull(Kids_First_Biospecimen_ID) %>% unique()
  expression_of_interest <- expression_data %>% dplyr::select(all_of(cohort_bsid)) %>% 
    tibble::rownames_to_column("GeneSymbol") %>% 
    dplyr::filter(GeneSymbol ==x) %>% 
    tibble::column_to_rownames("GeneSymbol") %>%
    t() %>% as.data.frame() %>%
    tibble::rownames_to_column("Kids_First_Biospecimen_ID")
  
  #annotate relevant clinical information to the expression data frame
  combined <- expression_of_interest %>% dplyr::left_join(cohort_df) %>%
    dplyr::rename(gene_of_interest = x)
  
  # get ready for calculating statistics
  combined_fixed <- combined %>% mutate(harmonized_diagnosis = gsub("-","_", harmonized_diagnosis))
  n_sample_hd <- combined_fixed %>% group_by(harmonized_diagnosis) %>% summarize(n=n())
  
  combined_fixed <- combined_fixed %>%
    dplyr::left_join(n_sample_hd) %>%
    dplyr::mutate(harmonized_diagnosis = paste0(harmonized_diagnosis, " ", "n=", n)) %>%
    dplyr::mutate(harmonized_diagnosis = paste0(short_name, "_", harmonized_diagnosis))
  # calculate statistics 
  res_aov <- aov(combined_fixed$gene_of_interest ~ combined_fixed$harmonized_diagnosis,
                 data = combined_fixed)
  summary(res_aov)
  anova_result <- summary(res_aov)[[1]] %>% 
    tibble::rownames_to_column() %>%
    mutate(gene = x)
  anova_results <- rbind(anova_results, anova_result)
  
  tukey.test<-TukeyHSD(res_aov)[1] %>% as.data.frame() %>% tibble::rownames_to_column() 
  colnames(tukey.test) <- c("Pair Compared", "DIFF", "LWR", "UPR", "P.adj")
  
  tukey_test_heatmap <- tukey.test %>% arrange(P.adj, descending = TRUE) %>%
    mutate(gene=x) %>%
    tidyr::separate(col = `Pair Compared`, into = c("group1", "group2"), sep = "-") %>%
    rename(p.adj = P.adj) 
  
  # generate data frame for adding statistics for violin plot
  tukey_test_violin <- tukey_test_heatmap %>%
    mutate(.y. = "harmonized_diagnosis") %>%
    dplyr::select(.y., group1, group2, p.adj) %>% tibble() %>%
    filter(p.adj <=0.01)
  tukey_test_violin <-tukey_test_violin %>%
    mutate(y.position = seq(0.6*max(combined_fixed$gene_of_interest), 1.5*max(combined_fixed$gene_of_interest), length.out = nrow(tukey_test_violin)))

  tukey_test_violin$p.adj <- round(tukey_test_violin$p.adj, digits=6)
  
  p<-combined_fixed %>%
    ggplot( aes(x=harmonized_diagnosis, y=gene_of_interest)) +
    geom_violin(width=1.4, trim=TRUE, show.legend = F, aes(fill=harmonized_diagnosis)) +
    geom_boxplot(width=0.1, color="black", show.legend = F,aes(fill=harmonized_diagnosis)) +
    labs(title=paste0(x," TPM per Harmonized Diagnosis"),x="Harmonized Diagnosis", y = paste0(x," TPM Value")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    stat_pvalue_manual(tukey_test_violin, label = "p.adj")
    
  tpm_plot_each_dir <- file.path(tpm_plots_dir,x)
  if (!dir.exists(tpm_plot_each_dir )) {
    dir.create(tpm_plot_each_dir , recursive = TRUE)
  }
  # ggsave(file.path(tpm_plot_each_dir,"harmonized_diagnosis_violin.png"), p, height = 10, width=12)
  
  violin_plot <- combined_fixed %>%
    group_by(short_name) %>% 
    ggplot( aes(x=harmonized_diagnosis, y=gene_of_interest)) +
    geom_violin(width=1.4, trim=TRUE, show.legend = F, aes(fill=harmonized_diagnosis)) +
    geom_boxplot(width=0.1, color="black", show.legend = F,aes(fill=harmonized_diagnosis)) +
    labs(title=paste0(x," TPM per Harmonized Diagnosis"),x="Harmonized Diagnosis", y = paste0(x," TPM Value")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(~ short_name, scales = "free", space = "free")
  
  ggsave(file.path(tpm_plot_each_dir,"harmonized_diagnosis_violin_no_stats.png"), violin_plot, height = 7, width=12)
  
 
  # Generate output
  tukey_results <- rbind(tukey_results, tukey.test)
  
  # Generate complete dataframe for plotting heatmap
  tukey_test_heatmap_filled <- tukey_test_heatmap %>%
    rename(group2=group1, 
           group1=group2)
  tukey_heatmap_complete <- rbind(tukey_test_heatmap, tukey_test_heatmap_filled) %>%
    mutate(group_1_short = gsub("\\_.*", "", group1)) %>%
    mutate(group_2_short = gsub("\\_.*", "", group2)) %>%
    # treat all p.adj>=.05 as non significant 
    mutate(p.adj = case_when(
      p.adj < 0.05 ~ p.adj,
      p.adj >=0.05 ~ 0.05
    ))
  
  q<-ggplot(tukey_heatmap_complete, aes(x=group1, group2, fill= p.adj)) +
          geom_tile(aes(fill=p.adj)) +
          labs(title=paste0(x," Harmonized Diagnosis Tukey Test for TPM"),x="Harmonized Diagnosis", y = "Harmonized Diagnosis") + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          scale_fill_gradient(low = "red", high = "white") + 
    facet_grid(~group_1_short, scales = "free", space = "free")
  
  ggsave(file.path(tpm_plot_each_dir,"harmonized_diagnosis_heatmap.png"), q, height = 10, width=10)
}

write_tsv(anova_results, file.path(stats_anova_dir, "anova_result_hd.tsv"))
write_tsv(tukey_results, file.path(stats_tukey_dir, "tukey_result_hd.tsv"))


  


