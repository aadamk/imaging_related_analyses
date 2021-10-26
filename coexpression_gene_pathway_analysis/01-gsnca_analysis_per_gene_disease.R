# Author: Run Jin
#
# GSNCA analysis comparing upper and lower quantile of gene expressions in each disease
# BiocManager::install("GSAR")
# BiocManager::install("GSVAdata")
# BiocManager::install("DGCA")
# BiocManager::install("EGSEA")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GSAR"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("EGSEA"))
suppressPackageStartupMessages(library("DGCA"))


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

#### Run GSNCA for each combination in the CG Gene match file-------------------
combined_results <- data.frame()
for(i in 1:nrow(cg_gene_interest)){
  # find the cancer group of interest
  cg_interest <- cg_gene_interest[i,1] %>% as.character()
  # find the gene of interest
  gene_interest <- cg_gene_interest[i,2] %>% as.character()
  # get the quantile of interest 
  quantile_interest <- cg_gene_interest[i,3] %>% as.character()
  
  plots_dir_specific <- file.path(plots_dir, paste0(cg_interest, "_", gene_interest, "_", quantile_interest))
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
  
  ######################### prepare to GSNCA test for all pathways-filter out low expression genes 
  # filter lowly expressed genes by DGSA
  eoi_coding_each_filtered <- filterGenes(expression_of_interest_coding_each, 
                                           filterTypes = c("central", "dispersion"),
                                           filterDispersionType = "cv", 
                                           filterDispersionPercentile = 0.2,
                                           sequential= TRUE)
  
  # need to additionally filter on sd to allow next steps
  group1_bsid <- bs_id_quantile_df %>% filter(group=="1") %>% rownames()
  group2_bsid <- bs_id_quantile_df %>% filter(group=="2") %>% rownames()
  
  eoi_group1 <- eoi_coding_each_filtered %>% dplyr::select(all_of(group1_bsid)) 
  eoi_group1_sd <- apply(eoi_group1, 1, sd, na.rm = TRUE)
  eoi_group1_sd_filter <- which(eoi_group1_sd <= 0.015) %>% as.data.frame() %>% rownames()
  
  eoi_group2 <- eoi_coding_each_filtered %>% dplyr::select(all_of(group2_bsid))
  eoi_group2_sd <- apply(eoi_group2, 1, sd, na.rm = TRUE)
  eoi_group2_sd_filter <- which(eoi_group2_sd <= 0.015) %>% as.data.frame() %>% rownames()
  
  eoi_coding_each_filtered <- eoi_coding_each_filtered %>% 
    tibble::rownames_to_column("gene_symbol") %>% 
    dplyr::filter(!gene_symbol %in% eoi_group1_sd_filter) %>%
    dplyr::filter(!gene_symbol %in% eoi_group2_sd_filter) %>%
    tibble::column_to_rownames("gene_symbol")
  
  ######## finally run GSNCA test for all pathways
  gsnca_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(gsnca_results) <- c("pathway_id", "pathway_description", "pvalue")
  # only look at pathways that has <=500 members
  pathway_list <- kegg_build_genesets_df %>% 
    group_by(pathway) %>% 
    mutate(n=n()) %>% 
    filter(n <=500) %>% 
    pull(pathway) %>% 
    unique()
  
  for(j in 1:length(pathway_list)){
    pathway_of_interest <- pathway_list[j]
    description <- kegg_build_genesets_df %>% 
      filter(pathway == pathway_of_interest) %>% 
      pull(description) %>% unique()
    
    # find genes in pathway of interest
    genes_in_pathway <- kegg_build_genesets_df %>% 
      filter(pathway == pathway_of_interest) %>% 
      pull(symbol) %>% unique()
    
    # filter to genes in target pathway
    eoi_target_pathway_each <- eoi_coding_each_filtered %>%
      tibble::rownames_to_column("gene_symbol") %>% 
      filter(gene_symbol %in% genes_in_pathway) %>% 
      tibble::column_to_rownames("gene_symbol") 
    
    # only run the test when there are >=10 genes in the pathway that is also in our expression table 
    if(nrow(eoi_target_pathway_each)>=10){
      result_pval<-GSNCAtest(object=as.matrix(eoi_target_pathway_each), 
                             # since the matrix is selected by order of row, the group will match
                             group=bs_id_quantile_df$group, 
                             nperm=1000, 
                             cor.method="spearman", 
                             check.sd=TRUE, 
                             min.sd=1e-3, 
                             max.skip=500
      )
      # store the pathway name and results in the results table
      gsnca_results[j,1] <- pathway_of_interest
      gsnca_results[j,3] <- result_pval
      gsnca_results[j,2] <- description
      
      #### For pathways that has <-0.05 significance, we can plot that----------------
      if(result_pval<0.05){
        pdf(file = file.path(plots_dir_specific, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_", pathway_of_interest, "_plot.pdf" )))
        plotMST2.pathway(object=as.matrix(eoi_target_pathway_each),
                         # since the matrix is selected by order of row, the group will match
                         group=bs_id_quantile_df$group,
                         cor.method="spearman",
                         legend.size=0.9,
                         label.size=1.2,
                         name=paste0(pathway_of_interest, " ", description))
        dev.off()
      }
    }
  }
  # write out results
  gsnca_results <- gsnca_results %>% 
    dplyr::filter(!is.na(pathway_id)) %>%
    arrange(pvalue, descending = FALSE) 
  gsnca_results %>% 
    readr::write_tsv(file.path(results_dir, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_pathway_analysis.tsv" )))
  # add gene and cancer group to the table
  gsnca_results <- gsnca_results %>% 
    mutate(cancer_group = cg_interest) %>% 
    mutate(gene_parsed_by = gene_interest) %>%
    mutate(percentile = quantile_interest)
  # combine the results to the larger data set
  combined_results <- bind_rows(combined_results, gsnca_results)
  
}
# write out combined results
readr::write_tsv(combined_results, file.path(results_dir, "combined_gscna_pathway_analysis.tsv" ))



