# Author: Run Jin
# Run ssGSEA analysis and perform differentially 
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("CEMiTool")
  library("BiocParallel")
  library("limma")
  library("edgeR")
  library("DESeq2")
  library("rtracklayer")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"),
  make_option(c("-t","--gmt_file"),type="character",
              help="gmt file containing the pathway of interest"),
  make_option(c("-c","--count"),type="character",
              help="gene count data from OpenPedCan RNA samples (.rds) "),
  make_option(c("-e","--expression"),type="character",
              help="gene expression TPM data from OpenPedCan RNA samples (.rds) "),
  make_option(c("-g","--gtf_file"),type="character",
              help="GTF file of the Gencode V27 primary assembly file")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cg_list <-unlist(strsplit(opt$cg_interest,","))
gmt_file <- opt$gmt_file

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "cluster_profiling")
anno_input_dir <- file.path(analysis_dir, "results", "cluster_anno")
count_input_dir <- file.path(analysis_dir, "results", "clustering")

results_dir <- file.path(analysis_dir, "results", "ssgsea")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "ssgsea_score_heatmap")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

source(file.path(analysis_dir, "utils", "perform_ssgsea.R"))

#### Read in files necessary for analyses --------------------------------------
# gene expected count file
count_matrix <- readRDS(opt$count)
# gene expected count file
tpm_df <- readRDS(opt$expression)

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
# filter expression count file to contain only protein coding gene
tpm_df_coding <- tpm_df[rownames(tpm_df) %in% gencode_gtf$gene_name,]

# read in the GMT file and generate list for analysis 
c2_cp_kegg <- read_gmt(gmt_file)
c2_cp_kegg_list <- base::split(c2_cp_kegg$gene, list(c2_cp_kegg$term))

# now perform analysis for each disease of interest 
for(i in 1:length(cg_list)){
  cg_of_interest <- cg_list[i]
  
  # define disease and cluster specific output 
  results_dir_specific <- file.path(results_dir, cg_of_interest)
  if(!dir.exists(results_dir_specific)){
    dir.create(results_dir_specific, recursive=TRUE)
  }
  
  # read in the cluster annotation file 
  cluster_anno_file <- list.files(anno_input_dir, pattern = cg_of_interest, full.names = TRUE)
  cluster_anno <- readr::read_tsv(cluster_anno_file) 
  # get list of clusters
  cluster_list <- cluster_anno %>% pull(cluster_assigned) %>% unique()
  
  # read in vst count file 
  vst_count_file <- readRDS(file.path(count_input_dir, cg_of_interest, "transformed_all_coding_counts.rds"))
  
  # filter the count matrix to contain only samples of interest
  count_matrix_coding_cg <- count_matrix_coding %>%
    dplyr::select(cluster_anno$Kids_First_Biospecimen_ID)
  # filter the TPM to contain only samples of interest
  tpm_df_coding_cg <- tpm_df_coding %>%
    dplyr::select(cluster_anno$Kids_First_Biospecimen_ID)
  
  ############### Try out different normalization method
  # first is TMM calculated by edgeR
  group <- factor(rep(c(1), times = length(colnames(count_matrix_coding_cg)))) # give them the same group for all genes
  count_matrix_coding_cg_dge <- DGEList(counts = count_matrix_coding_cg, group = group)

  # compute the normalization factors (TMM) for scaling by library size
  count_matrix_coding_cg_dge_tmm <- calcNormFactors(count_matrix_coding_cg_dge, method = "TMM")
  count_matrix_coding_cg_dge_tmm <- cpm(count_matrix_coding_cg_dge_tmm)

  # second is DESeq2 normalized count
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix_coding_cg),
                                colData = cluster_anno,
                                design = ~1 )
  count_matrix_coding_cg_dge_deseq2 <- estimateSizeFactors(dds)
  count_matrix_coding_cg_dge_deseq2 <- counts(count_matrix_coding_cg_dge_deseq2, normalized=T)

  # third is log TPM
  tpm_df_coding_cg_log2 <- log2(tpm_df_coding_cg + 1)
  
  ####### run the analysis for each normalized matrix
  ssgsea_analysis(normalized_count = count_matrix_coding_cg_dge_tmm,
                  normalized_method = "edgeR_tmm")

  ssgsea_analysis(normalized_count = count_matrix_coding_cg_dge_deseq2,
                  normalized_method = "deseq2")

  ssgsea_analysis(normalized_count = tpm_df_coding_cg_log2,
                  normalized_method = "log2_tpm")
  
  ssgsea_analysis(normalized_count = vst_count_file,
                  normalized_method = "vst_count")
  
}
