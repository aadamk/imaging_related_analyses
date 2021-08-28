# Author: Run Jin
#
# Obtain TPM plots for genes of interest with cancer group of interest

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("car"))
suppressPackageStartupMessages(library("multcomp"))
suppressPackageStartupMessages(library("kableExtra"))

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
cancer_group_long <- short_long_match %>% 
  filter(short_name %in% cancer_group_list) %>%
  pull(long_name)

cohort_df <- histology_df %>% dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
  dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
  dplyr::filter(cancer_group %in% cancer_group_long) %>%
  dplyr::select(Kids_First_Biospecimen_ID, harmonized_diagnosis)

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
  
  n_sample_hd <- combined %>% group_by(harmonized_diagnosis) %>% summarize(n=n())
  p<-combined %>%
    dplyr::left_join(n_sample_hd) %>%
    dplyr::mutate(myaxis = paste0(harmonized_diagnosis, "\n", "n=", n)) %>%
    ggplot( aes(x=myaxis, y=gene_of_interest, fill=harmonized_diagnosis)) +
    geom_violin(width=1.4, trim=FALSE, show.legend = F) +
    geom_boxplot(width=0.1, color="black", show.legend = F) +
    labs(title=paste0(x," TPM per Harmonized Diagnosis"),x="Harmonized Diagnosis", y = paste0(x," TPM Value")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  tpm_plot_each_dir <- file.path(tpm_plots_dir,x)
  if (!dir.exists(tpm_plot_each_dir )) {
    dir.create(tpm_plot_each_dir , recursive = TRUE)
  }
  ggsave(file.path(tpm_plot_each_dir,"harmonized_diagnosis_violin.png"), p)

  # calculate statistics 
  res_aov <- aov(combined$gene_of_interest ~ combined$harmonized_diagnosis,
                 data = combined)
  summary(res_aov)
  anova_result <- summary(res_aov)[[1]] %>% 
    tibble::rownames_to_column() %>%
    mutate(gene = x)
  anova_results <- rbind(anova_results, anova_result)
  
  tukey.test<-TukeyHSD(res_aov)[1] %>% as.data.frame() %>% rownames_to_column() 
  colnames(tukey.test) <- c("Pair Compared", "DIFF", "LWR", "UPR", "P.adj")
  
  tukey.test <- tukey.test %>% arrange(P.adj, descending = TRUE) %>%
    mutate(gene=x)
  tukey_results <- rbind(tukey_results, tukey.test)
  
}

write_tsv(anova_results, file.path(stats_anova_dir, "anova_result_hd.tsv"))
write_tsv(tukey_results, file.path(stats_tukey_dir, "tukey_result_hd.tsv"))


