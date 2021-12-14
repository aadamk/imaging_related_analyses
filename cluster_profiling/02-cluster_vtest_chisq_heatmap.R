# Author: Run Jin
# Generate heatmap with consensus clustering results + harmonized diagnosis
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("broom")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-c","--cc_data_match"),type="character",
              help="file matching cancer group to optimal cluster method and n (.tsv) "),
  make_option(c("-p","--pathways"),type="character",
              help="file containing subset of metabolic pathways of interest (.tsv) ")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "cluster_profiling")
input_dir <- file.path(analysis_dir, "results", "clustering")

results_dir_vtest <- file.path(analysis_dir, "results", "vtest")
if(!dir.exists(results_dir_vtest)){
  dir.create(results_dir_vtest, recursive=TRUE)
}

results_dir_cluster <- file.path(analysis_dir, "results", "cluster_anno")
if(!dir.exists(results_dir_cluster)){
  dir.create(results_dir_cluster, recursive=TRUE)
}

results_dir_chisq <- file.path(analysis_dir, "results", "chisq_test")
if(!dir.exists(results_dir_chisq)){
  dir.create(results_dir_chisq, recursive=TRUE)
}

heatmap_dir <- file.path(analysis_dir, "plots", "cluster_heatmap")
if(!dir.exists(heatmap_dir)){
  dir.create(heatmap_dir, recursive=TRUE)
}

source(file.path(analysis_dir, "utils", "calc_functions.R"))

#### Read in files necessary for analyses --------------------------------------
# histology file
histology_df <- readr::read_tsv(opt$histology, guess_max=100000)

# CC parameter selection
cc_data_match <- readr::read_tsv(opt$cc_data_match)

# use pathway file to get gene list
pathways <- readr::read_tsv(opt$pathways) 
genes_in_pathways <-  pathways$genes %>%
  str_split(",") %>% unlist() %>% unique()

# file matching short and long histology
short_long_match <- readr::read_tsv(opt$short_long_match)

######### now perform the analysis on each disease type of interest
for (i in 1:nrow(cc_data_match)){
  
  # obtain parameter information from the match file
  cg_of_interest <- cc_data_match[i,]$cancer_group
  distance_param <- cc_data_match[i,]$distance_param
  cluster_param <- cc_data_match[i,]$cluster_param
  cluster_n <- cc_data_match[i,]$cluster_n
  
  # read in relevant cluster information 
  input_file <- file.path(input_dir, 
                          cg_of_interest, 
                          paste0(distance_param, "_", cluster_param, "_CC.Rdata"))
  
  load(input_file)
  
  # extract cluster assignment for cluster params
  cluster_info <- result_CC[[cluster_n]][["consensusClass"]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column()
  # clean up
  colnames(cluster_info) <- c("sample", "cluster_assigned")
  
  # read in filtered count cg coding pval file from 01 script 
  input_count <- file.path(input_dir, 
                            cg_of_interest, 
                            "transformed_diptest_coding_goi_filtered_counts.rds")
  filtered_count <- readRDS(input_count)
  
  #################### calculate v test for GOI + 1k dip test features
  # annotate cluster information
  filtered_count_cluster_all <- filtered_count %>%
    t() %>% 
    as.data.frame() %>%
    tibble::rownames_to_column('sample') %>%
    tidyr::gather(gene_symbol, gene_count, -sample) %>%
    dplyr::left_join(cluster_info) 
  
  # gather information about genes and their count number
  filtered_count_cluster_all <- filtered_count_cluster_all %>%
    dplyr::group_by(cluster_assigned, gene_symbol) %>%
    dplyr::mutate(cluster_gene_mean_score = mean(gene_count)) %>% # mean of gene count per cluster 
    ungroup() %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::mutate(gene_mean_score = mean(gene_count),
                  gene_variance = var(gene_count)) # global mean & variance per gene count
  
  # apply v.test function per gene
  vtest_output_all <- plyr::ddply(.data = filtered_count_cluster_all, 
                                  .variables = "gene_symbol", 
                                  .fun = function(x) compute.v.test(x, clustering_col = "cluster_assigned"))
  
  vtest_output_all <- vtest_output_all %>%
    dplyr::rename("geneSymbol" = "gene_symbol", "vtest_score" = "v_score") %>% 
    readr::write_tsv(file.path(results_dir_vtest, paste0(cg_of_interest, "_", 
                                                         distance_param, "_", 
                                                         cluster_param, "_k", 
                                                         cluster_n, 
                                                         "_cluster_all.tsv")))
  
  #################### calculate v test for GOI in metabolic pathways only 
  # filter to genes of interest
  filtered_count_interest <- filtered_count[rownames(filtered_count) %in% genes_in_pathways,]
  # annotate cluster information
  filtered_count_cluster <- filtered_count_interest %>%
    t() %>% 
    as.data.frame() %>%
    tibble::rownames_to_column('sample') %>%
    tidyr::gather(gene_symbol, gene_count, -sample) %>%
    dplyr::left_join(cluster_info) 
  
  # gather information about genes and their count number
  filtered_count_cluster <- filtered_count_cluster_all %>%
    dplyr::filter(gene_symbol %in%  genes_in_pathways) 
  
  # apply v.test function per gene
  vtest_output <- plyr::ddply(.data = filtered_count_cluster, 
                              .variables = "gene_symbol", 
                              .fun = function(x) compute.v.test(x, clustering_col = "cluster_assigned"))
  
  vtest_output <- vtest_output %>%
    dplyr::rename("geneSymbol" = "gene_symbol", "vtest_score" = "v_score") %>% 
    readr::write_tsv(file.path(results_dir_vtest, paste0(cg_of_interest, "_", 
                                                          distance_param, "_", 
                                                          cluster_param, "_k", 
                                                          cluster_n, 
                                                          "_cluster_goi_only.tsv")))
  
  ############################### plot heatmap with consensus clustering results
  # take log and then zscore of the matrix for better visualization
  filtered_count_interest_zscored <- zscore_transform(filtered_count_interest)
  
  # generate annotation file 
  cluster_anno <- histology_df %>% 
    dplyr::filter(Kids_First_Biospecimen_ID %in% cluster_info$sample) %>%
    dplyr::rename(sample = Kids_First_Biospecimen_ID) %>%
    dplyr::select(sample, harmonized_diagnosis) %>%
    dplyr::left_join(cluster_info) %>% 
    arrange(cluster_assigned, harmonized_diagnosis) %>%
    tibble::column_to_rownames("sample")
  cluster_anno$cluster_assigned <- as.factor(cluster_anno$cluster_assigned)
  
  # output the cluster information 
  cluster_anno %>% 
    tibble::rownames_to_column('Kids_First_Biospecimen_ID') %>% 
    readr::write_tsv(file.path(results_dir_cluster, paste0(cg_of_interest, "_", 
                                                           distance_param, "_", 
                                                           cluster_param, "_k", 
                                                           cluster_n, 
                                                           "_cluster_info.tsv")))
  
  # arrange the z-score transformed matrix in the same order 
  filtered_count_interest_zscored <- filtered_count_interest_zscored %>%
    dplyr::select(rownames(cluster_anno))
  
  # define breaks for heatmap
  breaks<- seq(-3,3, length.out=94)
  breaks <- c(-6,-5,-4,breaks,4, 5, 6)
  
  # output heatmap
  pheatmap::pheatmap(as.matrix(filtered_count_interest_zscored),
                     annotation_col = cluster_anno,
                     cluster_rows=TRUE,
                     cluster_cols=FALSE,
                     breaks=breaks,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     width = 10,
                     height = 8,
                     show_colnames = F,
                     filename = file.path(heatmap_dir, paste0(cg_of_interest, "_", 
                                                              distance_param, "_", 
                                                              cluster_param, "_k", 
                                                              cluster_n, 
                                                              "_cluster_heatmap.pdf")))
  
  ######## perform chisq test for harmonized diagnosis and cluster assigned 
  # chisq test
  chisq_results <- tidy(chisq.test(as.factor(cluster_anno$harmonized_diagnosis), 
                                   as.factor(cluster_anno$cluster_assigned))) %>%
    as.data.frame() %>% 
    dplyr::mutate(compare = "harmonized_diagnosis_vs_cluster") %>% 
    readr::write_tsv(file.path(results_dir_chisq, paste0(cg_of_interest, 
                                                         "_chisq_results.tsv")))
}

