# function to perform ssgsea analysis and output plots
ssgsea_analysis <- function(normalized_count, 
                            normalized_method, 
                            pathway_df = c2_cp_kegg, 
                            pathway_list = c2_cp_kegg_list, 
                            anno_file = cluster_anno, 
                            level_list = cluster_list,
                            cg_results_dir = results_dir_specific,
                            cg_plots_dir = plots_dir_specific){
  
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
  
  ######## run differential pathway using contrast fit with combined gsea scores
  # define contrasts
  cluster_n <- length(unique(anno_file$cluster_assigned))
  contrasts <- c()
  # get all combinations
  for(p in 1:(cluster_n - 1)){
    for(q in (p+1):cluster_n){
      contrasts <- c(contrasts, paste0("cluster", p, "-cluster", q))
    }
  }
  
  # define column names
  column_names <- c()
  for(m in 1:cluster_n){
    column_names <- c(column_names, paste0("cluster", m))
  }
  
  # build model matrix
  mod <- model.matrix(~ factor(anno_file$cluster_assigned))
  colnames(mod) <- column_names
  fit_contrast <- lmFit(as.matrix(gsea_combined), mod)
  
  contrasts <- makeContrasts(contrasts=contrasts, levels=column_names)
  fit_contrast <- contrasts.fit(fit_contrast, contrasts)
  fit_contrast <- eBayes(fit_contrast)
  
  topTable(fit_contrast, n=Inf) %>% 
    tibble::rownames_to_column("pathway") %>%
    dplyr::select(pathway, P.Value, adj.P.Val) %>%
    dplyr::filter(adj.P.Val < 0.01) %>%
    dplyr::arrange(adj.P.Val) %>% 
    readr::write_tsv(file.path(results_dir_each,
                               paste0("enriched_pathway_in_", cg_of_interest, "_overall_results.tsv")))
  
  
  for(t in 1:ncol(as.data.frame(contrasts))){
    contrast_name <- colnames(as.data.frame(contrasts))[t]
    ssgsea_pathway_results <- topTable(fit_contrast, coef=t, n=Inf) %>% 
      dplyr::mutate(direction = ifelse((logFC>0), "up", "down")) %>%
      tibble::rownames_to_column("pathway") %>%
      dplyr::select(pathway, logFC, P.Value, adj.P.Val, direction) %>%
      dplyr::filter(adj.P.Val < 0.01) %>%
      dplyr::arrange(adj.P.Val) %>% 
      readr::write_tsv(file.path(results_dir_each,
                                 paste0("enriched_pathway_in_", cg_of_interest, "_comparing_", contrast_name, ".tsv")))
    
    ######## generate heatmap using the combined GSEA scores
    # get annotation file to have sample name as rownames
    anno_plot <- anno_file %>%
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
      dplyr::select(rownames(anno_plot))
    
    # output heatmap
    pheatmap::pheatmap(as.matrix(gsea_combined),
                       annotation_col = anno_plot,
                       annotation_row = pathway_anno, 
                       cluster_rows=TRUE,
                       cluster_cols=FALSE,
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       width = 12,
                       height = 8,
                       show_colnames = F,
                       fontsize_row = 4,
                       show_rownames = T,
                       filename = file.path(plots_dir_each,
                                            paste0("ssgsea_scores_in_", cg_of_interest, "_comparing_", contrast_name, "_by_cluster_heatmap.pdf")))
  }
}
