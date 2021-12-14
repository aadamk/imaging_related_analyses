# Author: Run Jin
# Perform GSNCA analysis across different clusters
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("CEMiTool")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-e","--expression"),type="character",
              help="gene expression TPM data from OpenPedCan RNA samples (.rds) "),
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
anno_input_dir <- file.path(analysis_dir, "results", "cluster_anno")

results_dir <- file.path(analysis_dir, "results", "gsnca_analysis")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "gsnca_analysis")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

source(file.path(analysis_dir, "utils", "gsnca_calc.R"))

#### Read in files necessary for analyses --------------------------------------
# gene expected count file
expression_df <- readRDS(opt$expression)

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
expression_df_coding <- expression_df[rownames(expression_df) %in% gencode_gtf$gene_name,]

# read in the GMT file and generate list for analysis 
c2_cp_kegg <- read_gmt(gmt_file)
c2_cp_kegg$term <- as.character(c2_cp_kegg$term)

# now perform analysis for each disease of interest 
for(j in 1:length(cg_list)){
  cg_of_interest <- cg_list[j]
  
  # define disease and cluster specific output 
  plots_dir_specific <- file.path(plots_dir, cg_of_interest)
  if(!dir.exists(plots_dir_specific)){
    dir.create(plots_dir_specific, recursive=TRUE)
  }
  
  results_dir_specific <- file.path(results_dir, cg_of_interest)
  if(!dir.exists(results_dir_specific)){
    dir.create(results_dir_specific, recursive=TRUE)
  }
  
  # read in the cluster annotation file 
  cluster_anno_file <- list.files(anno_input_dir, pattern = cg_of_interest, full.names = TRUE)
  cluster_anno <- readr::read_tsv(cluster_anno_file) 
  
  # filter to samples in cancer group
  expression_coding_cg_filtered <- expression_df_coding %>%
    dplyr::select(cluster_anno$Kids_First_Biospecimen_ID) %>%
    filter_low_expr_df()
  
  # iterate through all cluster pair combination to get results 
  cluster_n <- length(unique(cluster_anno$cluster_assigned))
  
  # get all combinations
  for(p in 1:(cluster_n - 1)){
    for(q in (p+1):cluster_n){
      group1_cluster <- p
      group2_cluster <- q
      
      # generate annotation file with only cluster of interest
      cluster_anno_each <- cluster_anno %>%
        dplyr::filter(cluster_assigned %in% c(p, q)) %>%
        dplyr::mutate(group = case_when(
          cluster_assigned == p ~ "1",
          cluster_assigned == q ~ "2"
        )) %>% 
        dplyr::select(Kids_First_Biospecimen_ID, group)
      
      # filter to these samples 
      expression_coding_cg_filtered_each <- expression_coding_cg_filtered %>%
        dplyr::select(cluster_anno_each$Kids_First_Biospecimen_ID)
      
      # run the analysis 
      gsnca_analysis_plot(tpm_matrix = as.matrix(expression_coding_cg_filtered_each), 
                           cluster_anno = cluster_anno_each, 
                           pathway_df = c2_cp_kegg, 
                           comparison = paste0("cluster", p, "_vs_cluster", q),
                           output_file_dir = results_dir_specific, 
                           output_plot_dir = plots_dir_specific, 
                           top_bar=20, 
                           top_net=5)
      
    }
  }
}

