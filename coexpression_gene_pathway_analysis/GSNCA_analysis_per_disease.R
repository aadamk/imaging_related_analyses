# Author: Run Jin
#
# GSNCA analysis comparing upper and lower quantile of gene expressions in each disease
BiocManager::install("GSAR")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GSAR"))
suppressPackageStartupMessages(library("edgeR"))


#### Parse command line options ------------------------------------------------

option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-e","--expression"),type="character",
              help="gene expression rsem tpm file for OpenPedCan RNA samples (.rds) "),
  make_option(c("-l","--cg_gene_interest"),type="character",
              help="file containing gene of interest and matching cancer group (.tsv)"),
  make_option(c("-m","--short_long_match"),type="character",
              help="match between long and short names (.tsv)"),
  make_option(c("-g","--gtf_file"),type="character",
              help="gtf file for annotation (.gtf.gz)")
)

opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

#### Read in files necessary for analyses --------------------------------------
histology_df <- readr::read_tsv(opt$histology, guess_max=100000)
expression_data <- readRDS(opt$expression)
short_long_match <- readr::read_tsv(opt$short_long_match)
cg_gene_interest <- readr::read_tsv(opt$cg_gene_interest)
gtf <- opt$gtf_file

#### Filter to the cancer group of interest ------------------------------------
cancer_group_list <- cg_gene_interest %>% 
  pull(short_name) %>% unique() %>%
  as.list()

cohort_df_list <- lapply(cancer_group_list, function(x){
  cancer_group_long <- short_long_match %>% 
    filter(short_name ==x) %>%
    pull(long_name)
  
  cohort_df_each <- histology_df %>% 
    # keep only initial CNV tumor
    dplyr::filter(sample_type == "Tumor") %>% 
    dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
    # do not include TCGA or GTEx
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

#### Subset gene count to sample of interest -----------------------------------
# get biospecimen ID's for samples 
cohort_bsid <- cohort_df %>% pull(Kids_First_Biospecimen_ID) %>% unique()

# subset to gene and sample of interest          
expression_of_interest <- expression_data %>% dplyr::select(all_of(cohort_bsid)) 

# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = gtf)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# filter expression count file to contain only protein coding gene
expression_of_interest_coding <- expression_of_interest[rownames(expression_of_interest) %in% gencode_gtf$gene_name,]

#### Run GSNCA for each combination in the CG Gene match file-------------------
for(i in 1:nrow(cg_gene_interest)){
  # find the cancer group of interest
  cg_interest <- cg_gene_interest[i,1] %>% 
    pull(short_name)
  # find the gene of interest
  gene_interest <- cg_gene_interest[i,2] %>% 
    pull(gene_of_interest)
  
  # get BS ID for that particular cancer group
  cg_bsid_each <- cohort_df %>%
    filter(short_name == cg_interest) %>%
    pull(Kids_First_Biospecimen_ID) %>% unique()
  
  # filter expression to that particular cancer group
  expression_of_interest_coding_each <- expression_of_interest_coding %>%
    select(all_of(cg_bsid_each))
  
  ######## filter out genes with sd < 0.001 and mean tpm < 0.1 for next step
  # calculate sd for each gene and get gene names of sd < 0.001
  gene_sd_low <- apply(expression_of_interest_coding_each, 1, sd, na.rm = TRUE) %>%
    as.data.frame() %>% 
    filter(. <0.001) %>% 
    rownames()
  
  # calculate mean for each gene and get gene names ofmean tpm < 0.1
  gene_mean_low <- rowMeans(expression_of_interest_coding_each, na.rm = TRUE) %>% 
    as.data.frame() %>% 
    filter(. <0.1) %>% 
    rownames()
  
  # filter those genes out
  expression_of_interest_coding_each <- expression_of_interest_coding_each %>%
    tibble::rownames_to_column("Gene_symbol") %>% 
    filter(!Gene_symbol %in% gene_sd_low) %>%
    filter(!Gene_symbol %in% gene_mean_low) %>%
    tibble::column_to_rownames("Gene_symbol")
  
  ######## assign 1 and 2 to all samples in the cancer group base on expression of goi
  expression_of_goi <- expression_of_interest_coding_each %>%
    tibble::rownames_to_column("Gene_symbol") %>% 
    filter(Gene_symbol == gene_interest) %>%
    tibble::column_to_rownames("Gene_symbol") %>%
    t() %>% 
    as.data.frame() 
  colnames(expression_of_goi) <-"gene_interest"
  
  # calculate the quantile and assign groups
  upper_quantile <- quantile(expression_of_goi$gene_interest, 0.75)
  lower_quantile <- quantile(expression_of_goi$gene_interest, 0.25)
  
  bs_id_quantile_df <- expression_of_goi %>%
    tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>% 
    mutate(group = case_when(
      gene_interest >= upper_quantile ~"1", 
      gene_interest <= lower_quantile ~"2",
      TRUE ~ "middle"
    )) %>% 
    filter(group != "middle") %>%
    select(Kids_First_Biospecimen_ID, group) %>% 
    tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  
  # filter count matrix to those samples
  expression_of_interest_coding_each_final <- expression_of_interest_coding_each %>% 
    # by selecting, the matrix is sorted in column based on row order o bs_id_quantile_df
    dplyr::select(rownames(bs_id_quantile_df)) %>%
    as.matrix()
  
  ######## finally run GSNCA test
  results<-GSNCAtest(object=expression_of_interest_coding_each_final, 
                     # since the matrix is selected by order of row, the group will match
                     group=bs_id_quantile_df$group, 
                     nperm=1000, 
                     cor.method="pearson", 
                     check.sd=TRUE, 
                     min.sd=1e-3, 
                     max.skip=10)
}

