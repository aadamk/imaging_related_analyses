# Author: Run Jin
# GSVA analysis comparing upper and lower quantile of gene expressions in each disease
# BiocManager::install("BiocParallel")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GSVA"))
suppressPackageStartupMessages(library("KEGGREST"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("BiocParallel"))

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
  make_option(c("-o","--out_file"),type="character",
              help="output file with GSVA scores of each pathway in each disease (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
outfile <- opt$out_file

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "coexpression_gene_pathway_analysis")
results_dir <- file.path(analysis_dir, "results", "gsva_scores")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
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
# get KEGG pathway and entrez gene ids
hsa_path_eg  <- keggLink("pathway", "hsa") %>% 
  tibble(pathway = ., eg = sub("hsa:", "", names(.)))

#annotated with the SYMBOL and ENSEMBL identifiers associated with each Entrez id
hsa_kegg_anno <- hsa_path_eg %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
    ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
  )
# get pathway description
hsa_pathways <- keggList("pathway", "hsa") %>% 
  tibble(pathway = names(.), description = .)

# combine to get a complete table
hsa_kegg_anno_pathways <- left_join(hsa_kegg_anno, hsa_pathways)

# generate a list of pathways and their matching genes for GSVA analysis 
pathway_names <- hsa_kegg_anno_pathways %>% 
  pull(pathway) %>% unique() %>% as.list()

pathway_list <- lapply(pathway_names, function(x){
  genes <- hsa_kegg_anno_pathways %>% 
    filter(pathway == x) %>% 
    pull(symbol) %>% unique()
  return(genes)
})
names(pathway_list) <- pathway_names 

################## Run GSVA scores for each gene, cancer group combination
gsva_scores_df_tidy <- data.frame()
for(i in 1:nrow(cg_gene_interest)){
  # find the cancer group of interest
  cg_interest <- cg_gene_interest[i,1] %>% 
    pull(short_name)
  # find the gene of interest
  gene_interest <- cg_gene_interest[i,2] %>% 
    pull(gene_of_interest)
  
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
  upper_quantile <- quantile(expression_of_goi$gene_interest, 0.75)
  lower_quantile <- quantile(expression_of_goi$gene_interest, 0.25)
  
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
  gsva_scores_each <- GSVA::gsva(expression_of_interest_coding_each_log2_matrix,
                                 pathway_list,
                                 method = "gsva",
                                 min.sz=1, max.sz=1500,## Arguments from K. Rathi
                                 parallel.sz = 8, # For the bigger dataset, this ensures this won't crash due to memory problems
                                 mx.diff = TRUE,
                                 BPPARAM=SerialParam(progressbar=T))        ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
  
  ### Clean scoring into tidy format
  gsva_scores_each_df <- as.data.frame(gsva_scores_each) %>%
    rownames_to_column(var = "pathway") 
  
  #first/last_bs needed for use in gather (we are not on tidyr1.0)
  first_bs <- head(colnames(gsva_scores_each), n=1)
  last_bs  <- tail(colnames(gsva_scores_each), n=1)
  
  gsva_scores_each_df_tidy <- gsva_scores_each_df %>%
    tidyr::gather(Kids_First_Biospecimen_ID, gsva_score, !!first_bs : !!last_bs) %>%
    dplyr::select(Kids_First_Biospecimen_ID, pathway, gsva_score)  %>% 
    dplyr::left_join(hsa_pathways) %>%
    dplyr::mutate(cancer_group = cg_interest) %>% 
    dplyr::mutate(gene_parsed_by = gene_interest)
  
  # add group number fo the final table as well
  bs_id_quantile_df <- bs_id_quantile_df %>% tibble::rownames_to_column("Kids_First_Biospecimen_ID")
  gsva_scores_each_df_tidy <- gsva_scores_each_df_tidy %>%
    dplyr::left_join(bs_id_quantile_df)
  # write out individual scores
  readr::write_tsv(gsva_scores_each_df_tidy, file.path(results_dir, paste0(cg_interest, "_parsed_by_", gene_interest, "_gsva_scores.tsv")))
  
  # merge into a combined file
  gsva_scores_df_tidy <-  bind_rows(gsva_scores_df_tidy , gsva_scores_each_df_tidy)
}

# write out results
readr::write_tsv(gsva_scores_df_tidy, outfile)

