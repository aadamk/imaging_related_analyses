# function to perform ssgsea analysis and output plots
ssgsea_analysis <- function(normalized_count, 
                            normalized_method,
                            pathway_df = c2_cp_kegg, 
                            pathway_list = c2_cp_kegg_list, 
                            anno_file = cluster_anno, 
                            level_list = cluster_list,
                            cg_results_dir = results_dir_specific,
                            cg_plots_dir = plots_dir){
  
  # define directories
  results_dir_each <- file.path(cg_results_dir, normalized_method)
  if(!dir.exists(results_dir_each)){
    dir.create(results_dir_each, recursive=TRUE)
  }
  
  plots_dir_each <- file.path(cg_plots_dir, normalized_method)
  if(!dir.exists(plots_dir_each)){
    dir.create(plots_dir_each, recursive=TRUE)
  }
  
  # define dataframe for writing out results - rownames is pathway name
  gsea_combined <- pathway_df$term %>% unique() %>% as.data.frame() 
  colnames(gsea_combined) <- "term"
  gsea_combined <- gsea_combined %>% 
    tibble::column_to_rownames("term")
  
  # run ssGSEA analysis on each cluster 
  for(j in 1:length(level_list)){
    # get samples in cluster
    cluster_of_interest <- level_list[j]
    samples_in_cluster <- anno_file %>% 
      dplyr::filter(cluster_assigned == cluster_of_interest) %>% 
      pull(Kids_First_Biospecimen_ID) 
    
    # filter expression matrix to contain only those samples
    normalized_count_each_cluster <- normalized_count %>% 
      as.data.frame() %>% 
      dplyr::select(all_of(samples_in_cluster))
    
    # ssGSEA analysis
    gsea_scores_df <- GSVA::gsva(as.matrix (normalized_count_each_cluster),
                                 pathway_list,
                                 method = "ssgsea",
                                 parallel.sz = 8, # For the bigger dataset, this ensures this won't crash due to memory problems
                                 mx.diff = TRUE,
                                 ssgsea.norm = T,
                                 BPPARAM=SerialParam(progressbar=T)) %>%
      as.data.frame()
    
    # write out the results 
    gsea_scores_df %>%
      tibble::rownames_to_column("pathway_description") %>% 
      readr::write_tsv(file.path(results_dir_each, 
                                 paste0("ssgsea_scores_in_", cg_of_interest, "_cluster_", cluster_of_interest, ".tsv")))
    
    # combine results from each cluster 
    gsea_combined <- cbind(gsea_combined, gsea_scores_df)
  }
  
  ######## run differential gene expression using the combined gsea scores
  # build model matrix
  mod <- model.matrix(~ anno_file$cluster_assigned)
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
    readr::write_tsv(file.path(results_dir_each,
                               paste0("enriched_pathway_in_", cg_of_interest, ".tsv")))
  
  ######## generate heatmap using the combined GSEA scores
  # get annotation file to have sample name as rownames
  anno_file <- anno_file %>%
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
    dplyr::select(rownames(anno_file))
  
  # output heatmap
  pheatmap::pheatmap(as.matrix(gsea_combined),
                     annotation_col = anno_file,
                     annotation_row = pathway_anno, 
                     cluster_rows=TRUE,
                     cluster_cols=FALSE,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     width = 12,
                     height = 8,
                     show_colnames = F,
                     fontsize_row = 6,
                     show_rownames = T,
                     filename = file.path(plots_dir_each,
                                          paste0("ssgsea_scores_in_", cg_of_interest, "_by_cluster_heatmap.pdf")))
}
