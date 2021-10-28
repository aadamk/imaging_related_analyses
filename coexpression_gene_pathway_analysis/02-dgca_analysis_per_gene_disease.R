# Author: Run Jin
#
# GSNCA analysis comparing upper and lower quantile of gene expressions in each disease
# BiocManager::install("GSAR")
# BiocManager::install("GSVAdata")
# BiocManager::install("DGCA")
# BiocManager::install("GOstats")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GSAR"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("DGCA"))
suppressPackageStartupMessages(library(GOstats))
suppressPackageStartupMessages(library(HGNChelper))


#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-e","--expression"),type="character",
              help="gene expression rsem tpm file for OpenPedCan RNA samples (.rds) "),
  make_option(c("-l","--cg_gene_interest"),type="character",
              help="file containing gene of interest and matching cancer group (.tsv)"),
  make_option(c("-m","--short_long_match"),type="character",
              help="match between long and short names (.tsv)"),
  make_option(c("-g","--gtf_file"),type="character",
              help="gtf file for annotation (.gtf.gz)"),
  make_option(c("-o","--outfile"),type="character",
              help="path of the combined GO term analysis results (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "coexpression_gene_pathway_analysis")
score_results_dir <- file.path(analysis_dir, "results", "dgca_scores")
if(!dir.exists(score_results_dir)){
  dir.create(score_results_dir, recursive=TRUE)
}

go_term_results_dir <- file.path(analysis_dir, "results", "dgca_go_term")
if(!dir.exists(go_term_results_dir)){
  dir.create(go_term_results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "dgca_plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

#### Read in files necessary for analyses --------------------------------------
# histology file
histology_df <- readr::read_tsv(opt$histology, guess_max=100000)

# gene expression TPM file
expression_data <- readRDS(opt$expression)

# file matching short and long histology
short_long_match <- readr::read_tsv(opt$short_long_match)
# file containing genes of interest in each cancer group
cg_gene_interest <- readr::read_tsv(opt$cg_gene_interest)

# Gencode27 GTF file loading
gtf <- opt$gtf_file

#### Filter to the cancer group of interest ------------------------------------
cancer_group_list <- cg_gene_interest %>% 
  pull(short_name) %>% unique() %>%
  as.list()

cohort_df_list <- lapply(cancer_group_list, function(x){
  cancer_group_long <- short_long_match %>% 
    filter(short_name ==x) %>%
    pull(long_name)
  
  cohort_df_each <- histology_df %>% 
    # keep only initial CNV tumor
    dplyr::filter(sample_type == "Tumor") %>% 
    dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
    # do not include TCGA or GTEx
    dplyr::filter(cohort %in% c("PBTA","GMKF","TARGET")) %>% 
    dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
    dplyr::filter(cancer_group %in% cancer_group_long) %>%
    # exclude derived cell line
    dplyr::filter(composition != "Derived Cell Line") %>% 
    # keep only unique Kids First Participant 
    dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>% 
    dplyr::select(Kids_First_Biospecimen_ID, harmonized_diagnosis) %>%
    dplyr::mutate(short_name = x)
})

cohort_df <- do.call(rbind,cohort_df_list)

#### Subset gene count to sample of interest -----------------------------------
# get biospecimen ID's for samples 
cohort_bsid <- cohort_df %>% pull(Kids_First_Biospecimen_ID) %>% unique()

# subset to gene and sample of interest          
expression_of_interest <- expression_data %>% dplyr::select(cohort_bsid)

# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = gtf)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# filter expression count file to contain only protein coding gene
expression_of_interest_coding <- expression_of_interest[rownames(expression_of_interest) %in% gencode_gtf$gene_name,]

#### Run DGCA for the entire ------------------
# combined_results <- data.frame()

for(i in 1:nrow(cg_gene_interest)){
  # find the cancer group of interest
  cg_interest <- cg_gene_interest[i,1] %>% as.character()
  # find the gene of interest
  gene_interest <- cg_gene_interest[i,2] %>% as.character()
  # get the quantile of interest 
  quantile_interest <- cg_gene_interest[i,3] %>% as.character()
  
  # get BS ID for that particular cancer group
  cg_bsid_each <- cohort_df %>%
    filter(short_name == cg_interest) %>%
    pull(Kids_First_Biospecimen_ID) %>% unique()
  
  # filter expression to that particular cancer group
  expression_of_interest_coding_each <- expression_of_interest_coding %>%
    dplyr::select(cg_bsid_each)
  
  ######## assign 1 and 2 to all samples in the cancer group base on expression of goi
  expression_of_goi <- expression_of_interest_coding_each %>%
    tibble::rownames_to_column("Gene_symbol") %>% 
    filter(Gene_symbol == gene_interest) %>%
    tibble::column_to_rownames("Gene_symbol") %>%
    t() %>% 
    as.data.frame() 
  colnames(expression_of_goi) <-"gene_interest"
  
  # calculate the quantile and assign groups
  upper_quantile <- quantile(expression_of_goi$gene_interest, (1-(as.numeric(quantile_interest)/100)))
  lower_quantile <- quantile(expression_of_goi$gene_interest, (as.numeric(quantile_interest)/100))
  
  bs_id_quantile_df <- expression_of_goi %>%
    tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>% 
    mutate(group = case_when(
      gene_interest >= upper_quantile ~"upper", 
      gene_interest <= lower_quantile ~"lower",
      TRUE ~ "middle"
    )) %>% 
    filter(group != "middle") %>%
    dplyr::select(Kids_First_Biospecimen_ID, group) %>% 
    tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  
  # filter count matrix to those samples
  expression_of_interest_coding_each <- expression_of_interest_coding_each %>% 
    # by selecting, the matrix is sorted in column based on row order o bs_id_quantile_df
    dplyr::select(unlist(rownames(bs_id_quantile_df)))
  
  ######################### prepare to DGCA test for all pathways-filter out low expression genes 
  # filter lowly expressed genes by DGCA
  eoi_coding_each_filtered <- filterGenes(expression_of_interest_coding_each, 
                                          filterTypes = c("central", "dispersion"),
                                          filterDispersionType = "cv", 
                                          filterDispersionPercentile = 0.2,
                                          sequential= TRUE)
  
  ######## finally run DGCA test 
  # prepare design matrix
  design_matrix <- makeDesign(as.vector(bs_id_quantile_df$group))
  
  # define directory for the heatmap
  ddcor_res_100 <- ddcorAll(inputMat =as.matrix(eoi_coding_each_filtered),
                           design = design_matrix,
                           compare = c("upper", "lower"),
                           adjust = "none",
                           heatmapPlot = F,
                           nPerm = 0,
                           corrType = "spearman",
                           # specific filtering parameter for heatmap
                           filterCentralPercentile = 0.75, 
                           filterDispersionPercentile = 0.75,
                           nPairs=100)
  ddcor_res_100 %>%
    readr::write_tsv(file.path(score_results_dir, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_dgca_100_scores.tsv.gz" )))
  
  # write out all the results
  ddcor_res_all <- ddcorAll(inputMat =as.matrix(eoi_coding_each_filtered),
                            design = design_matrix,
                            compare = c("upper", "lower"),
                            adjust = "none",
                            heatmapPlot = FALSE,
                            nPerm = 0,
                            corrType = "spearman",
                            # specific filtering parameter for generating all results 
                            filterCentralPercentile = 0.3, 
                            filterDispersionPercentile = 0.3)

  ddcor_res_all %>%
    readr::write_tsv(file.path(score_results_dir, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_dgca_full_scores.tsv.gz" )))
  
  # ################# run the analysis for gene of interest 
  # if(!gene_interest %in% rownames(eoi_coding_each_filtered)){
  #   gene_rescue <- expression_of_interest_coding_each[gene_interest,]
  #   eoi_coding_each_filtered <- bind_rows(eoi_coding_each_filtered, gene_rescue)
  # }
  #   
  # ddcor_res_goi <- ddcorAll(inputMat =as.matrix(eoi_coding_each_filtered), 
  #                           design = design_matrix, 
  #                           compare = c("upper", "lower"), 
  #                           adjust = "none", 
  #                           heatmapPlot = FALSE, 
  #                           nPerm = 0, 
  #                           corrType = "spearman", 
  #                           splitSet = gene_interest)
  # 
  # # Generate GO results of enriched differential correlations by pathway - 
  # ddcorGO_res <-ddcorGO(ddcor_res_goi, 
  #                       universe = rownames(eoi_coding_each_filtered), 
  #                       gene_ontology = "all", 
  #                       HGNC_clean = TRUE, 
  #                       HGNC_switch = TRUE, 
  #                       annotation = "org.Hs.eg.db", 
  #                       calculateVariance = TRUE)
  # # write out results as RDS
  # ddcorGO_res %>% 
  #   saveRDS(file.path(go_term_results_dir, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_GO_by_dgca.rds" )))
  # 
  # # extract only GO term results 
  # gain_bp <- ddcorGO_res[[3]][[1]] %>% as.data.frame() %>% dplyr::mutate(change_dir = "gain_of_correlation_genes") %>%
  #   dplyr::rename(GOID = GOBPID) %>% filter(Pvalue < 0.05)
  # gain_mf <- ddcorGO_res[[3]][[2]] %>% as.data.frame() %>% dplyr::mutate(change_dir = "gain_of_correlation_genes") %>%
  #   dplyr::rename(GOID = GOMFID) %>%  filter(Pvalue < 0.05)
  # gain_cc <- ddcorGO_res[[3]][[3]] %>% as.data.frame() %>% dplyr::mutate(change_dir = "gain_of_correlation_genes") %>%
  #   dplyr::rename(GOID = GOCCID) %>% filter(Pvalue < 0.05)
  # 
  # loss_bp <- ddcorGO_res[[4]][[1]] %>% as.data.frame() %>% dplyr::mutate(change_dir = "loss_of_correlation_genes") %>%
  #   dplyr::rename(GOID = GOBPID) %>% filter(Pvalue < 0.05)
  # loss_mf <- ddcorGO_res[[4]][[2]] %>% as.data.frame() %>% dplyr::mutate(change_dir = "loss_of_correlation_genes") %>%
  #   dplyr::rename(GOID = GOMFID) %>% filter(Pvalue < 0.05)
  # loss_cc <- ddcorGO_res[[4]][[3]] %>% as.data.frame() %>% dplyr::mutate(change_dir = "loss_of_correlation_genes") %>%
  #   dplyr::rename(GOID = GOCCID) %>% filter(Pvalue < 0.05)
  # 
  #   
  # # combine only GO term results and write out as one file 
  # combined <- bind_rows(gain_bp, gain_mf, gain_cc, loss_bp, loss_mf, loss_cc) 
  # 
  # combined %>% 
  #   mutate(cancer_group = cg_interest) %>% 
  #   mutate(gene_parsed_by = gene_interest) %>% 
  #   mutate(percentile = quantile_interest) %>% 
  #   readr::write_tsv(file.path(go_term_results_dir, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_combined_GO_by_dgca.tsv" )))
  # 
  # combined_results <- bind_rows(combined_results, combined )
}

# # write out combined GO term results
# combined_results %>% 
#   readr::write_tsv(file.path(go_term_results_dir, opt$outfile))



