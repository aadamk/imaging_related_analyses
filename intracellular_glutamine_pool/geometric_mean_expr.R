# Author: Run Jin
# Calculate geometric mean of all markers within the manually curated pathways
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-e","--expression"),type="character",
              help="gene expression tpm data from OpenPedCan RNA samples (.rds) "),
  make_option(c("-l","--cg_gene_interest"),type="character",
              help="file containing gene of interest and matching cancer group (.tsv)"),
  make_option(c("-m","--short_long_match"),type="character",
              help="match between long and short names (.tsv)"),
  make_option(c("-p","--pathway"),type="character",
              help="manually curated pathways and genes of interest (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "intracellular_glutamine_pool")
results_dir <- file.path(analysis_dir, "results")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

#### Read in files necessary for analyses --------------------------------------
# histology file
histology_df <- readr::read_tsv(opt$histology, guess_max=100000) %>% 
  # keep only initial CNV tumor
  dplyr::filter(sample_type == "Tumor") %>% 
  dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
  # do not include TCGA or GTEx
  dplyr::filter(cohort %in% c("PBTA","GMKF","TARGET")) %>% 
  dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
  # exclude derived cell line
  dplyr::filter(composition != "Derived Cell Line") %>% 
  # keep only unique Kids First Participant 
  dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE)

# gene expression TPM file
expression_data <- readRDS(opt$expression)

# file matching short and long histology
short_long_match <- readr::read_tsv(opt$short_long_match)
# file containing genes of interest in each cancer group
cg_gene_interest <- readr::read_tsv(opt$cg_gene_interest)

# file containing pathways of interest and their genes
pathway_df <- readr::read_tsv(opt$pathway)

#### Filter to the cancer group of interest ------------------------------------
pval_geometric <- data.frame(matrix(ncol=4))
colnames(pval_geometric) <- c("short_name", "gene", "quantile", "pval")

for(i in 1:nrow(cg_gene_interest)){
  
  # get the information of cohort of interest 
  short_name_interest <- cg_gene_interest[i,]$short_name
  gene_of_interest <- cg_gene_interest[i,]$gene_of_interest
  quantile_interest <- cg_gene_interest[i,]$percentile
  
  # find the matching cancer group name 
  cancer_group_long <- short_long_match %>% 
    filter(short_name == short_name_interest) %>%
    pull(long_name)
  
  # find the histology associated with the cancer group of interest 
  cohort_df_each <- histology_df  %>% 
    dplyr::filter(cancer_group %in% cancer_group_long) %>% 
    dplyr::select(Kids_First_Biospecimen_ID, harmonized_diagnosis) %>%
    dplyr::mutate(short_name = short_name_interest)
  
  #### Handle gene expression matrix -----------------------------------
  # first filter to our cohort of interest 
  expression_of_interest <- expression_data[,cohort_df_each$Kids_First_Biospecimen_ID]
  
  # annotate high vs. low expression group to the samples in the cohort 
  expression_of_gi <- expression_of_interest[rownames(expression_of_interest) == gene_of_interest, ] %>%
    t() %>% as.data.frame()
  colnames(expression_of_gi) <- "gene_of_interest"
  
  # calculate the quantile and assign groups
  upper_quantile <- quantile(expression_of_gi$gene_of_interest, (1-(as.numeric(quantile_interest)/100)))
  lower_quantile <- quantile(expression_of_gi$gene_of_interest, (as.numeric(quantile_interest)/100))
  
  bs_id_quantile_df <- expression_of_gi %>%
    tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>% 
    mutate(group = case_when(
      gene_of_interest >= upper_quantile ~"high", 
      gene_of_interest <= lower_quantile ~"low",
      TRUE ~ "middle"
    )) %>% 
    filter(group != "middle") %>%
    dplyr::select(Kids_First_Biospecimen_ID, group) %>% 
    dplyr::left_join(cohort_df_each) 
  
  # subset to gene and sample of interest          
  expression_of_interest <- expression_data %>% 
    dplyr::select(bs_id_quantile_df$Kids_First_Biospecimen_ID) 
  
  #### Calculate geometric mean for markers of interest --------------------------------
  # get markers of interest
  markers <- pathway_df$genes %>%
    str_split(",") %>% unlist() %>% unique()
  
  # filter to only markers
  expression_markers <- expression_of_interest[rownames(expression_of_interest) %in% markers,]
  
  #calculate geometric mean 
  geometric_mean <- lapply(expression_markers,function(x)exp(mean(log(x)))) %>% 
    unlist() %>% as.data.frame()
  colnames(geometric_mean) <- "geometric_mean"
  
  # annotate the results with low vs. high group
  geometric_mean <- geometric_mean %>%
    tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
    left_join(bs_id_quantile_df) %>% 
    readr::write_tsv(file.path(results_dir, 
                               paste0("geometric_mean_of_markers_parsed_by_", gene_of_interest, "_in_", short_name_interest, ".tsv")))
  
  # get t.test results
  high_exp <- geometric_mean %>% filter(group=="high") %>% dplyr::select(geometric_mean)
  low_exp <- geometric_mean %>% filter(group=="low") %>% dplyr::select(geometric_mean)
  pvalue <- t.test(high_exp, low_exp)$p.value
  
  # store the results
  pval_geometric[i,1] <- short_name_interest
  pval_geometric[i,2] <- gene_of_interest
  pval_geometric[i,3] <- quantile_interest
  pval_geometric[i,4] <- pvalue
  
}

readr::write_tsv(pval_geometric, file.path(results_dir, "geometric_mean_pval_markers.tsv"))



