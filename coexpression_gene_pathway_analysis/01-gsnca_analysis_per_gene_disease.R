# Author: Run Jin
#
# GSNCA analysis comparing upper and lower quantile of gene expressions in each disease
# BiocManager::install("GSAR")
# BiocManager::install("GSVAdata")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GSAR"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("KEGGREST"))

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
              help="gtf file for annotation (.gtf.gz)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "coexpression_gene_pathway_analysis")
results_dir <- file.path(analysis_dir, "results", "gsnca_pval")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "gsnca_plots")
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
# get pathway discription
hsa_pathways <- keggList("pathway", "hsa") %>% 
  tibble(pathway = names(.), description = .)

# combine to get a complete table
hsa_kegg_anno_pathways <- left_join(hsa_kegg_anno, hsa_pathways)

#### Run GSNCA for each combination in the CG Gene match file-------------------
combined_results <- data.frame()
for(i in 1:nrow(cg_gene_interest)){
  # find the cancer group of interest
  cg_interest <- cg_gene_interest[i,1] %>% 
    pull(short_name)
  # find the gene of interest
  gene_interest <- cg_gene_interest[i,2] %>% 
    pull(gene_of_interest)
  
  plots_dir_specific <- file.path(plots_dir, paste0(cg_interest, "_", gene_interest))
  if(!dir.exists(plots_dir_specific)){
    dir.create(plots_dir_specific)
  }
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
  
  ######################### prepare to GSNCA test for all pathways-filter out genes with sd < 0.01 and tpm <0.1 for each group 
  # separate expression to group 1 and group2
  bs_id_group1 <- bs_id_quantile_df %>% filter(group == "1") %>% rownames()
  bs_id_group2 <- bs_id_quantile_df %>% filter(group == "2") %>% rownames()
  expression_of_interest_coding_each_group1 <- expression_of_interest_coding_each %>% 
    dplyr::select(all_of(bs_id_group1))
  expression_of_interest_coding_each_group2 <- expression_of_interest_coding_each %>% 
    dplyr::select(all_of(bs_id_group2))
  
  # calculate sd for each gene and get gene names of sd < 0.0015
  gene_sd_low1 <- apply(expression_of_interest_coding_each_group1, 1, sd, na.rm = TRUE) %>%
    as.data.frame() %>% 
    # <0.01 will be filtered since GSNCA round to 3 digits
    filter(. <0.01) %>% 
    rownames()
  
  gene_sd_low2 <- apply(expression_of_interest_coding_each_group2, 1, sd, na.rm = TRUE) %>%
    as.data.frame() %>% filter(. <0.01) %>% rownames()
  
  gene_mean_low1 <- rowMeans(expression_of_interest_coding_each_group1, na.rm = TRUE) %>%
    as.data.frame() %>% filter(. <0.1) %>% rownames()
  
  gene_mean_low2 <- rowMeans(expression_of_interest_coding_each_group2, na.rm = TRUE) %>%
    as.data.frame() %>% filter(. <0.1) %>% rownames()
  
  # filter those genes out
  expression_of_interest_coding_each <- expression_of_interest_coding_each %>%
    tibble::rownames_to_column("Gene_symbol") %>% 
    filter(!Gene_symbol %in% gene_sd_low1) %>%
    filter(!Gene_symbol %in% gene_sd_low2) %>%
    filter(!Gene_symbol %in% gene_mean_low1) %>%
    filter(!Gene_symbol %in% gene_mean_low1) %>%
    tibble::column_to_rownames("Gene_symbol")

  ######## finally run GSNCA test for all pathways
  gsnca_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(gsnca_results) <- c("pathway_id", "pathway_description", "pvalue")
  # only look at pathways that has <=500 members
  pathway_list <- hsa_kegg_anno_pathways %>% 
    group_by(pathway) %>% 
    mutate(n=n()) %>% 
    filter(n <=500) %>% 
    pull(pathway) %>% 
    unique()
  
  for(j in 1:length(pathway_list)){
    pathway_of_interest <- pathway_list[j]
    description <- hsa_kegg_anno_pathways %>% 
      filter(pathway == pathway_of_interest) %>% 
      pull(description) %>% unique()
    
    # find genes in pathway of interest
    genes_in_pathway <- hsa_kegg_anno_pathways %>% 
      filter(pathway == pathway_of_interest) %>% 
      pull(symbol) %>% unique()
    
    # filter to genes in target pathway
    expr_target_pathway_each <- expression_of_interest_coding_each %>%
      tibble::rownames_to_column("gene_symbol") %>% 
      filter(gene_symbol %in% genes_in_pathway) %>% 
      tibble::column_to_rownames("gene_symbol") 
    
    # only run the test when there are >=10 genes in the pathway that is also in our expression table 
    if(nrow(expr_target_pathway_each)>=10){
      result_pval<-GSNCAtest(object=as.matrix(expr_target_pathway_each), 
                             # since the matrix is selected by order of row, the group will match
                             group=bs_id_quantile_df$group, 
                             nperm=1000, 
                             cor.method="pearson", 
                             check.sd=TRUE, 
                             min.sd=1e-3, 
                             max.skip=500
      )
      # store the pathway name and results in the results table
      gsnca_results[j,1] <- pathway_of_interest
      gsnca_results[j,2] <- result_pval
      gsnca_results[j,3] <- description
      
      #### For pathways that has <-0.05 significance, we can plot that----------------
      if(result_pval<0.05){
        pdf(file = file.path(plots_dir_specific, paste0(cg_interest, "_parsed_by_", gene_interest, "_", pathway_of_interest, "_plot.pdf" )))
        plotMST2.pathway(object=as.matrix(expr_target_pathway_each),
                         # since the matrix is selected by order of row, the group will match
                         group=bs_id_quantile_df$group,
                         cor.method="pearson",
                         legend.size=0.9,
                         label.size=1.2,
                         name=paste0(pathway_of_interest, " ", description))
        dev.off()
      }
    }
  }
  # write out results
  gsnca_results %>% 
    readr::write_tsv(file.path(results_dir, paste0(cg_interest, "_parsed_by_", gene_interest, "_pathway_analysis.tsv" )))
  # add gene and cancer group to the table
  gsnca_results <- gsnca_results %>% 
    mutate(cancer_group = cg_interest) %>% 
    mutate(gene_parsed_by = gene_interest) 
  # combine the results to the larger data set
  combined_results <- bind_rows(combined_results, gsnca_results)
  
}
# write out combined results
readr::write_tsv(combined_results, file.path(results_dir, "combined_pathway_analysis.tsv" ))



