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
results_dir <- file.path(analysis_dir, "results", "mboost")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
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
    left_join(expression_data_of_interest)
  
  # generate gene variables
  gene_variables <- rownames(expression_data) %>% paste(collapse="+")

  # perform glmboost 
  glm_fit <- glmboost(data = combined_data, 
                      as.formula(paste0("survival::Surv(OS_days, os_status_level) ~ ", gene_variables)),
                      family = CoxPH(), 
                      control = boost_control())
  
}

