# Author: Run Jin
# Run ssGSEA analysis and perform differentially 
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("CEMiTool")
  library("BiocParallel")
  library("limma")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"),
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
count_input_dir <- file.path(analysis_dir, "results", "clustering")

results_dir <- file.path(analysis_dir, "results", "ssgsea")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "ssgsea_score_heatmap")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

#### Read in files necessary for analyses --------------------------------------
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
  
  # read in the VST transformed protein coding file 
  count_file <- readRDS(file.path(count_input_dir, cg_of_interest, "transformed_all_coding_counts.rds"))
  # log transform for the analysis 
  count_log2_matrix <- log2(count_file + 1)
  
  # define dataframe for writing out results - rownames is pathway name
  gsea_combined <- c2_cp_kegg$term %>% unique() %>% as.data.frame() 
  colnames(gsea_combined) <- "term"
  gsea_combined <- gsea_combined %>% 
    tibble::column_to_rownames("term")
  
  # run ssGSEA analysis on each cluster 
  for(j in 1:length(cluster_list)){
    # get samples in cluster
    cluster_of_interest <- cluster_list[j]
    samples_in_cluster <- cluster_anno %>% 
      dplyr::filter(cluster_assigned == cluster_of_interest) %>% 
      pull(Kids_First_Biospecimen_ID) 
    
    # filter expression matrix to contain only those samples
    count_log2_matrix_each_cluster <- count_log2_matrix %>% 
      dplyr::select(all_of(samples_in_cluster))
    
    # ssGSEA analysis
    gsea_scores_df <- GSVA::gsva(as.matrix (count_log2_matrix_each_cluster),
                                  c2_cp_kegg_list,
                                  method = "ssgsea",
                                  parallel.sz = 8, # For the bigger dataset, this ensures this won't crash due to memory problems
                                  mx.diff = TRUE,
                                  ssgsea.norm = F,
                                  BPPARAM=SerialParam(progressbar=T)) %>%
      as.data.frame()
    
    # write out the results 
    gsea_scores_df %>%
      tibble::rownames_to_column("pathway_description") %>% 
      readr::write_tsv(file.path(results_dir_specific, paste0("ssgsea_scores_in_", cg_of_interest, "_cluster_", cluster_of_interest, ".tsv")))
    
    # combine results from each cluster 
    gsea_combined <- cbind(gsea_combined, gsea_scores_df)
  }
  
  ######## run differential gene expression using the combined gsea scores
  # build model matrix
  mod <- model.matrix(~ cluster_anno$cluster_assigned)
  fit <- lmFit(as.matrix(gsea_combined), mod)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=2, n=Inf) 
  
  # generate results with directions
  ssgsea_pathway_results <- tt %>% 
    dplyr::mutate(direction = ifelse((logFC>0), "up", "down")) %>%
    tibble::rownames_to_column("pathway") %>%
    dplyr::select(pathway, logFC, P.Value, adj.P.Val, direction) %>%
    dplyr::filter(adj.P.Val < 0.01) %>%
    dplyr::arrange(adj.P.Val) %>% 
    readr::write_tsv(file.path(results_dir_specific, paste0("enriched_pathway_in_", cg_of_interest, ".tsv")))

  ######## generate heatmap using the combined GSEA scores
  # get annotation file to have sample name as rownames
  cluster_anno <- cluster_anno %>%
    arrange(cluster_assigned, harmonized_diagnosis) %>% 
    tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  
  # filter to pathways with sig (adj.p < 0.01)
  gsea_combined <- gsea_combined[rownames(gsea_combined) %in% ssgsea_pathway_results$pathway,]
  
  # get annotation row
  pathway_anno <- ssgsea_pathway_results %>%
    dplyr::select(pathway, adj.P.Val, direction) %>%
    tibble::column_to_rownames("pathway")
  
  # arrange the GSEA score in the same order 
  gsea_combined <- gsea_combined %>%
    dplyr::select(rownames(cluster_anno))
  
  # output heatmap
  pheatmap::pheatmap(as.matrix(gsea_combined),
                     annotation_col = cluster_anno,
                     annotation_row = pathway_anno, 
                     cluster_rows=FALSE,
                     cluster_cols=FALSE,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     width = 12,
                     height = 8,
                     show_colnames = F,
                     fontsize_row = 6,
                     show_rownames = T,
                     filename = file.path(plots_dir, paste0("ssgsea_scores_in_", 
                                                            cg_of_interest, 
                                                            "_by_cluster_heatmap.pdf")))
  
}
