# Author: Run Jin
# Perform ssGSEA analysis for particular pathways within cancer group of interest
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("BiocParallel")
  library("GSVA")
  library("EGSEA")
  library("org.Hs.eg.db")
  library("GGally")
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
results_dir <- file.path(analysis_dir, "results", "ssgsea")

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "ssgsea")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

#### Read in files necessary for analyses --------------------------------------
# histology file
histology_df <- readr::read_tsv(opt$histology, guess_max=100000) %>%
  # keep only initial CNV tumor
  dplyr::filter(sample_type == "Tumor") %>% 
  dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
  # does not include TCGA or GTEx
  dplyr::filter(cohort %in% c("PBTA","GMKF","TARGET")) %>% 
  dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
  # exclude derived cell line
  dplyr::filter(composition != "Derived Cell Line") 

# gene expression TPM file
expression_data <- readRDS(opt$expression)

# file matching short and long histology
short_long_match <- readr::read_tsv(opt$short_long_match)

####### filter gene expression to contain genes of interest 
## Entrez IDs for SLC1A5 and ACLY
entrezs = c('6510', '47')
reac.annots <- buildMSigDBIdx(entrezIDs = entrezs, species = "human", geneSets = "c2")[["c2"]]@original

# pathways of interest 
pathway_of_interest <- c("REACTOME_TRIGLYCERIDE_BIOSYNTHESIS",
                         "REACTOME_FATTY_ACYL_COA_BIOSYNTHESIS",
                         "REACTOME_AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE",
                         "REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT",
                         "REACTOME_AMINO_ACID_AND_OLIGOPEPTIDE_SLC_TRANSPORTERS",
                         "REACTOME_METABOLISM_OF_LIPIDS_AND_LIPOPROTEINS",
                         "REACTOME_FATTY_ACID_TRIACYLGLYCEROL_AND_KETONE_BODY_METABOLISM")

# find the gene ids that are associated with the entrez ID in pathway of interest 
pathway_entrez <- reac.annots[pathway_of_interest] 

# define a df to store the results 
pathway_entrez_list <- lapply(pathway_entrez, function(x){
  x <- x %>% mapIds(org.Hs.eg.db, ., "SYMBOL", "ENTREZID") })

# all genes in the pathways of interest
gene_symbols <- pathway_entrez_list %>% unlist %>% unique()

# filter to contain only genes of interest
expression_data <- expression_data %>%
  # also filter to contain only samples in all cancer groups of interest 
  dplyr::select(histology_df$Kids_First_Biospecimen_ID) %>% 
  tibble::rownames_to_column("gene_symbol") %>% 
  dplyr::filter(gene_symbol %in% gene_symbols) %>% 
  tibble::column_to_rownames("gene_symbol")

####### do the analysis on all the cancer group of interest ----------------------------------
for (i in 1:length(cg_list)){
  # find the cancer group of interest 
  x <- cg_list[i]
  
  # match the long name to the short name
  long_name <- short_long_match %>% filter(short_name == x) %>%
    pull(long_name) 
  
  # filter to the cohort of interest
  cohort_bsid <- histology_df %>% dplyr::filter(cancer_group %in% long_name) %>%
    # get biospecimen ID's for samples 
    pull(Kids_First_Biospecimen_ID) %>% unique()
 
  # find the tpm of target genes for these samples 
  expression_data_of_interest <- expression_data %>%
    dplyr::select(all_of(cohort_bsid)) 
  
  # calculate z score for next steps
  expression_data_of_interest_log2_matrix <- as.matrix( log2(expression_data_of_interest  + 1) )
  
  # We then calculate the Gaussian-distributed scores
  ssgsea_scores_each <- GSVA::gsva(expression_data_of_interest_log2_matrix,
                                   pathway_entrez_list,
                                   method = "ssgsea",
                                   min.sz=1, ## has to use 1 since this is a manual list 
                                   max.sz=500,## Arguments from OMPARE
                                   parallel.sz = 8, # For the bigger dataset, this ensures this won't crash due to memory problems
                                   mx.diff = TRUE,
                                   BPPARAM=SerialParam(progressbar=T)) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("pathway_name")
  
  # save the results 
  ssgsea_scores_each %>% 
    readr::write_tsv(file.path(results_dir, paste0("ssgsea_scores_per_", x, ".tsv")))
  
  
}
