# Author: Run Jin
# GSVA analysis comparing upper and lower quantile of gene expressions in each disease
# BiocManager::install("BiocParallel")
# BiocManager::install("EGSEA")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GSVA"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("EGSEA"))

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
  make_option(c("-o","--out_file_score"),type="character",
              help="output file with SSGSEA scores of each pathway in each disease (.tsv)"),
  make_option(c("-p","--out_file_pval"),type="character",
              help="output file with pvals pathway significance based of SSGSEA scores per gene-disease (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
outfile_score <- opt$out_file_score
outfile_pval <- opt$out_file_pval

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "coexpression_gene_pathway_analysis")
scores_results_dir <- file.path(analysis_dir, "results", "ssgsea_scores")
if(!dir.exists(scores_results_dir)){
  dir.create(scores_results_dir, recursive = TRUE)
}

pval_results_dir <- file.path(analysis_dir, "results", "ssgsea_pval")
if(!dir.exists(pval_results_dir)){
  dir.create(pval_results_dir, recursive = TRUE)
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
expression_of_interest <- expression_data %>% dplyr::select(all_of(cohort_bsid)) 

# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = gtf)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# filter expression count file to contain only protein coding gene
expression_of_interest_coding <- expression_of_interest[rownames(expression_of_interest) %in% gencode_gtf$gene_name,]

#### Get the pathway information and select genes in the pathways---------------
# get entrezID for gene in expression pathway
expression_of_interest_coding_anno <- expression_of_interest_coding %>% 
  tibble::rownames_to_column("symbol") %>% 
  dplyr::select(symbol) %>%
  dplyr::mutate(eg = mapIds(org.Hs.eg.db, symbol, "ENTREZID", "SYMBOL")) 

# get GSCollectionSet object
kegg_build <- buildKEGGIdx(entrezIDs = expression_of_interest_coding_anno$eg, species = "human")

# get annotation to filter out disease related
kegg_build_anno <- kegg_build@anno %>% as.data.frame() %>% 
  dplyr::filter(Type != "Disease") %>% 
  dplyr::select(c("ID", "GeneSet", "Type")) %>%
  dplyr::rename(description = GeneSet, pathway=ID)

# get gene set IDs
kegg_build_genesets <- kegg_build_anno %>% pull(description)

# build df with pathways and entrezID of genes
kegg_build_genesets_list <- lapply(kegg_build_genesets, function(x){
  kegg_df <- kegg_build@idx[[x]] %>% as.data.frame() %>%
    mutate(description = x)
})
kegg_build_genesets_df <- do.call(rbind, kegg_build_genesets_list) %>%
  dplyr::left_join(kegg_build_anno) 
colnames(kegg_build_genesets_df) <- c("eg", "description", "pathway", "type")
kegg_build_genesets_df$eg <- as.character(kegg_build_genesets_df$eg)

# annotate gene symbol and ensemble IDs 
kegg_build_genesets_df <-kegg_build_genesets_df %>%
  dplyr::left_join(expression_of_interest_coding_anno) %>%
  # filter the pathway file to contain only symbols available in expression
  dplyr::filter(!is.na(symbol))

# generate a list of pathways and their matching genes for GSVA analysis 
pathway_names <- kegg_build_genesets_df %>% 
  pull(pathway) %>% unique() %>% as.list()

pathway_list <- lapply(pathway_names, function(x){
  genes <- kegg_build_genesets_df %>% 
    filter(pathway == x) %>% 
    pull(symbol) %>% unique()
  return(genes)
})
names(pathway_list) <- pathway_names 

################## Run GSVA scores for each gene, cancer group combination
ssgsea_scores_df_tidy <- data.frame()

for(i in 1:nrow(cg_gene_interest)){
  # find the cancer group of interest
  cg_interest <- cg_gene_interest[i,1] %>% 
    pull(short_name)
  # find the gene of interest
  gene_interest <- cg_gene_interest[i,2] %>% 
    pull(gene_of_interest)
  # get the quantile of interest 
  quantile_interest <- cg_gene_interest[i,3] %>% as.character()
  
  # get BS ID for that particular cancer group
  cg_bsid_each <- cohort_df %>%
    filter(short_name == cg_interest) %>%
    pull(Kids_First_Biospecimen_ID) %>% unique()
  
  # filter expression to that particular cancer group
  expression_of_interest_coding_each <- expression_of_interest_coding %>%
    dplyr::select(all_of(cg_bsid_each))
  
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
      gene_interest >= upper_quantile ~"1", 
      gene_interest <= lower_quantile ~"2",
      TRUE ~ "middle"
    )) %>% 
    filter(group != "middle") %>%
    dplyr::select(Kids_First_Biospecimen_ID, group) %>% 
    tibble::column_to_rownames("Kids_First_Biospecimen_ID")
  
  # filter count matrix to those samples
  expression_of_interest_coding_each <- expression_of_interest_coding_each %>% 
    # by selecting, the matrix is sorted in column based on row order o bs_id_quantile_df
    dplyr::select(unlist(rownames(bs_id_quantile_df)))
  # calculate z score for next steps
  expression_of_interest_coding_each_log2_matrix <- as.matrix( log2(expression_of_interest_coding_each + 1) )
  
  # We then calculate the Gaussian-distributed scores
  ssgsea_scores_each <- GSVA::gsva(expression_of_interest_coding_each_log2_matrix,
                                 pathway_list,
                                 method = "ssgsea",
                                 min.sz=5, 
                                 max.sz=500,## Arguments from K. Rathi
                                 parallel.sz = 8, # For the bigger dataset, this ensures this won't crash due to memory problems
                                 mx.diff = TRUE,
                                 BPPARAM=SerialParam(progressbar=T))        ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
  
  ### Clean scoring into tidy format
  ssgsea_scores_each_df <- as.data.frame(ssgsea_scores_each) %>%
    rownames_to_column(var = "pathway") 
  
  #first/last_bs needed for use in gather (we are not on tidyr1.0)
  first_bs <- head(colnames(ssgsea_scores_each), n=1)
  last_bs  <- tail(colnames(ssgsea_scores_each), n=1)

  ssgsea_scores_each_df_tidy <- ssgsea_scores_each_df %>%
    tidyr::gather(Kids_First_Biospecimen_ID, ssgsea_score, !!first_bs : !!last_bs) %>%
    dplyr::select(Kids_First_Biospecimen_ID, pathway, ssgsea_score)  %>% 
    dplyr::left_join(kegg_build_anno) %>%
    dplyr::mutate(cancer_group = cg_interest) %>% 
    dplyr::mutate(gene_parsed_by = gene_interest)
  
  # add group number fo the final table as well
  bs_id_quantile_df <- bs_id_quantile_df %>% tibble::rownames_to_column("Kids_First_Biospecimen_ID")
  ssgsea_scores_each_df_tidy <- ssgsea_scores_each_df_tidy %>%
    dplyr::left_join(bs_id_quantile_df)
  # write out individual scores
  ssgsea_scores_each_df_tidy %>% 
    arrange(ssgsea_score, descending = F) %>% 
    readr::write_tsv(file.path(scores_results_dir, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_ssgsea_scores.tsv")))
  
  # merge into a combined file
  ssgsea_scores_df_tidy <-  bind_rows(ssgsea_scores_df_tidy , ssgsea_scores_each_df_tidy)
}

# write out results
readr::write_tsv(ssgsea_scores_df_tidy, outfile_score)

########## Calculate t test statistics for scores 
combined_results <- data.frame()
for (i in 1:nrow(cg_gene_interest)){
  # find the cancer group of interest
  cg_interest <- cg_gene_interest[i,1] %>% 
    pull(short_name)
  # find the gene of interest
  gene_interest <- cg_gene_interest[i,2] %>% 
    pull(gene_of_interest)
  # get the quantile of interest 
  quantile_interest <- cg_gene_interest[i,3] %>% as.character()
  
  ssgsea_score_each <- ssgsea_scores_df_tidy %>% 
    filter(cancer_group == cg_interest) %>% 
    filter(gene_parsed_by == gene_interest) 
  
  # gather list of each pathway 
  pathway_list <- ssgsea_score_each %>% 
    pull(pathway) %>% unique
  
  ssgsea_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(ssgsea_results) <- c("pathway_id", "description","pval")
  for(j in 1:length(pathway_list)){
    pathway_interest <- pathway_list[j]
    
    # get description of pathway
    description <- ssgsea_score_each %>%  filter(pathway == pathway_interest) %>% 
      pull(description) %>% unique()
    
    # filter to contain only pathway of interest
    ssgsea_score_each_per_pathway <- ssgsea_score_each %>% 
      filter(pathway == pathway_interest)
    # group1 scores
    group1 <- ssgsea_score_each_per_pathway %>% filter(group == "1") %>% 
      pull(ssgsea_score)
    # group2 scores
    group2 <- ssgsea_score_each_per_pathway %>% filter(group == "2") %>% 
      pull(ssgsea_score)
    
    # compute t.test of the two
    res <- t.test(group1, group2)$p.value
    ssgsea_results[j,1] <- pathway_interest
    ssgsea_results[j,2] <- description
    ssgsea_results[j,3] <- res
    
    # write out individual file
    ssgsea_results %>% 
      arrange(pval, descending = F) %>%
      readr::write_tsv(file.path(pval_results_dir, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_ssgsea_pval.tsv" )))
  }
  # annotate gene and cancer group to the result
  ssgsea_results <- ssgsea_results %>% 
    mutate(cancer_group = cg_interest) %>% 
    mutate(gene_parsed_by = gene_interest)%>%
    mutate(percentile = quantile_interest)
  # combine them to a combined tsv file
  combined_results <- rbind(combined_results, ssgsea_results)
}

# write out combined results
combined_results %>%
  arrange(pval, descending = F) %>%
  readr::write_tsv(outfile_pval)

