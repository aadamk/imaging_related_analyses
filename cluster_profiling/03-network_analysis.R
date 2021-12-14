# Author: Run Jin
# Perform CEMiTools analysis for pathways 
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("CEMiTool")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-c","--count"),type="character",
              help="gene count data from OpenPedCan RNA samples (.rds) "),
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"),
  make_option(c("-g","--gtf_file"),type="character",
              help="GTF file of the Gencode V27 primary assembly file"),
  make_option(c("-t","--gmt_file"),type="character",
              help="gmt file containing the pathway of interest")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cg_list <-unlist(strsplit(opt$cg_interest,","))
gmt_file <- opt$gmt_file

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "cluster_profiling")
input_dir <- file.path(analysis_dir, "results", "cluster_anno")

results_dir <- file.path(analysis_dir, "results", "network_analysis")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "network_analysis")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

source(file.path(analysis_dir, "utils", "run_cemitools.R"))

#### Read in files necessary for analyses --------------------------------------
# gene expected count file
count_matrix <- readRDS(opt$count)

# Gencode27 GTF file loading
gtf <- opt$gtf_file

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
  
  # define disease and cluster specific output 
  results_dir_specific <- file.path(results_dir, cg_of_interest)
  if(!dir.exists(results_dir_specific)){
    dir.create(results_dir_specific, recursive=TRUE)
  }
  
  plots_dir_specific <- file.path(plots_dir, cg_of_interest)
  if(!dir.exists(plots_dir_specific)){
    dir.create(plots_dir_specific, recursive=TRUE)
  }
  
  # read in the cluster annotation file 
  cluster_anno_file <- list.files(input_dir, pattern = cg_of_interest, full.names = TRUE)
  cluster_anno <- readr::read_tsv(cluster_anno_file) %>% 
  # also make rownames kids first biospecimen id
    dplyr::mutate(tmp = Kids_First_Biospecimen_ID) %>% 
    tibble::column_to_rownames("tmp")
    
  
  ################# run CEMiTools on each cluster 
  # filter gene count matrix to only contain samples from this disease
  count_of_interest <- count_matrix_coding %>%
    dplyr::select(cluster_anno$Kids_First_Biospecimen_ID)
  
  # run the tool on particular cluster
  run_cemitools_functions(expr_df=count_of_interest, 
                          annot_df=cluster_anno, 
                          n = 100, 
                          direction = c("signed", "unsigned"), 
                          output_dir = results_dir_specific, 
                          plots_dir = plots_dir_specific,
                          gmt_file = gmt_file)
}
