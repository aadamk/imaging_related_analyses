# apply v.test function per gene
vtest_output <- plyr::ddply(.data = filtered_count_cg_coding_pval_cluster, 
                            .variables = "gene_symbol", 
                            .fun = function(x) compute.v.test(x, clustering_col = "cluster_assigned_nb"))

vtest_output <- vtest_output %>%
  dplyr::rename("geneSymbol" = "gene_symbol", "score" = "v_score") %>% 
  saveRDS(vtest_file)