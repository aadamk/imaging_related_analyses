# Author: Run Jin
# Perform consensus clustering on expected count of disease of interest 
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("DESeq2")
  library("ConsensusClusterPlus")
  })

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-c","--count"),type="character",
              help="gene count data from OpenPedCan RNA samples (.rds) "),
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"),
  make_option(c("-m","--short_long_match"),type="character",
              help="match between long and short names (.tsv)"),
  make_option(c("-g","--gtf_file"),type="character",
              help="GTF file of the Gencode V27 primary assembly file"),
  make_option(c("-p","--pathways"),type="character",
              help="File contianing all the metabolic related pathways and genes (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cg_list <-unlist(strsplit(opt$cg_interest,","))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "cluster_profiling")
results_dir <- file.path(analysis_dir, "results", "clustering")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "clustering")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

source(file.path(analysis_dir, "utils", "calc_functions.R"))

#### Read in files necessary for analyses --------------------------------------
# histology file
histology_df <- readr::read_tsv(opt$histology, guess_max=100000)

# gene expected count file
count_matrix <- readRDS(opt$count)

# file matching short and long histology
short_long_match <- readr::read_tsv(opt$short_long_match)

# Gencode27 GTF file loading
gtf <- opt$gtf_file

# use pathway file to get gene list
pathways <- readr::read_tsv(opt$pathways) 
genes_in_pathways <-  pathways$genes %>%
  str_split(",") %>% unlist() %>% unique()

# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = gtf)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# filter expression count file to contain only protein coding gene
count_matrix_coding <- count_matrix[rownames(count_matrix) %in% gencode_gtf$gene_name,]

# now perform analysis for each disease of interest 
for(i in 1:length(cg_list)){
  cg_of_interest <- cg_list[i]
  
  # define disease specific pathways 
  results_dir_specific <- file.path(results_dir, cg_of_interest)
  if(!dir.exists(results_dir_specific)){
    dir.create(results_dir_specific, recursive=TRUE)
  }
  
  plots_dir_specific <- file.path(plots_dir, cg_of_interest)
  if(!dir.exists(plots_dir_specific)){
    dir.create(plots_dir_specific, recursive=TRUE)
  }
  
  cancer_group_long <- short_long_match %>% 
    filter(short_name ==cg_of_interest ) %>%
    pull(long_name)
  
  cohort_df_each <- histology_df %>% 
    # keep only initial CNV tumor
    dplyr::filter(sample_type == "Tumor", 
                  tumor_descriptor == "Initial CNS Tumor") %>% 
    # do not include TCGA or GTEx
    dplyr::filter(cohort %in% c("PBTA","GMKF","TARGET")) %>% 
    dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
    dplyr::filter(cancer_group %in% cancer_group_long) %>%
    # exclude derived cell line
    dplyr::filter(composition != "Derived Cell Line") %>% 
    # keep only unique Kids First Participant 
    dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>% 
    dplyr::select(Kids_First_Biospecimen_ID, harmonized_diagnosis) %>%
    dplyr::mutate(short_name = cg_of_interest )
  
  # filter gene count matrix to only contain samples from this disease
  count_of_interest <- count_matrix_coding %>%
    dplyr::select(cohort_df_each$Kids_First_Biospecimen_ID)
  
  # perform variance stabilizing transformation 
  count_transformed <- varianceStabilizingTransformation(round(as.matrix(count_of_interest)),
                                                         blind = TRUE,
                                                         fitType="parametric") %>%
    as.data.frame()
  # filter by diptest
  diptest_gene_list <- get_gene_list_by_diptest(count_transformed, 
                                                min_n=1000)
  # add back the genes of interest if not present already
  combined_gene_list <- c(diptest_gene_list, genes_in_pathways) %>% 
    unique()
  
  # select the df with ombined gene list
  count_transformed_combined <- count_transformed[rownames(count_transformed) %in% combined_gene_list,]
  
  # save the counts file for later analysis 
  saveRDS(count_transformed_combined, 
          file.path(results_dir_specific, "transformed_diptest_coding_goi_filtered_counts.rds"))

  # perform consensus clustering on all possible parameters
  for(distance_param in c("spearman", "euclidean", "manhattan")){
    for(cluster_param in c("pam", "hc")){
      result_CC <- ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(count_transformed_combined),
                                                              finalLinkage = "average",
                                                              distance = distance_param,
                                                              clusterAlg = cluster_param,
                                                              plot = "pdf",
                                                              reps = 100, maxK = 10, pItem = 0.8,
                                                              title = file.path(plots_dir_specific, paste0(distance_param, "_", cluster_param)),
                                                              seed = 123)
      # save the plots for manual inspection
      save(result_CC, file = file.path(results_dir_specific, paste0(distance_param, "_", cluster_param, "_CC.Rdata")))
    }
  }
  
  # additionally, run on eucilean and km
  result_CC <- ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(count_transformed_combined),
                                                          finalLinkage = "average",
                                                          distance = "euclidean",
                                                          clusterAlg = "km",
                                                          plot = "pdf",
                                                          reps = 100, maxK = 10, pItem = 0.8,
                                                          title = file.path(plots_dir_specific, "euclidea_km"),
                                                          seed = 123)
  # save the plots for manual inspection
  save(result_CC, file = file.path(results_dir_specific, "euclidean_km_CC.Rdata"))
}
