# Author: Run Jin
# Calculate GSVA scores for all manually curated pathways
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GSVA"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("limma"))

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-e","--expression"),type="character",
              help="gene expression tpm data from OpenPedCan RNA samples (.rds) "),
  make_option(c("-l","--cg_gene_interest"),type="character",
              help="file containing gene of interest and matching cancer group (.tsv)"),
  make_option(c("-m","--short_long_match"),type="character",
              help="match between long and short names (.tsv)"),
  make_option(c("-p","--pathway"),type="character",
              help="manually curated pathways and genes of interest (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "intracellular_glutamine_pool")
results_dir <- file.path(analysis_dir, "results", "ssgsea_output")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots", "pathway_barplot")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

#### Read in files necessary for analyses --------------------------------------
# histology file
histology_df <- readr::read_tsv(opt$histology, guess_max=100000) %>% 
  # keep only initial CNV tumor
  dplyr::filter(sample_type == "Tumor") %>% 
  dplyr::filter(tumor_descriptor == "Initial CNS Tumor") %>%
  # do not include TCGA or GTEx
  dplyr::filter(cohort %in% c("PBTA","GMKF","TARGET")) %>% 
  dplyr::filter(experimental_strategy=="RNA-Seq") %>% 
  # exclude derived cell line
  dplyr::filter(composition != "Derived Cell Line") %>% 
  # keep only unique Kids First Participant 
  dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE)

# gene expression TPM file
expression_data <- readRDS(opt$expression)

# file matching short and long histology
short_long_match <- readr::read_tsv(opt$short_long_match)
# file containing genes of interest in each cancer group
cg_gene_interest <- readr::read_tsv(opt$cg_gene_interest)

# file containing pathways of interest and their genes
pathway_df <- readr::read_tsv(opt$pathway)

#### Filter to the cancer group of interest ------------------------------------

ssgsea_scores_df_tidy <- data.frame()
combined_results <- data.frame()

for(i in 1:nrow(cg_gene_interest)){
  
  # get the information of cohort of interest 
  short_name_interest <- cg_gene_interest[i,]$short_name
  gene_of_interest <- cg_gene_interest[i,]$gene_of_interest
  quantile_interest <- cg_gene_interest[i,]$percentile
  
  # find the matching cancer group name 
  cancer_group_long <- short_long_match %>% 
    filter(short_name == short_name_interest) %>%
    pull(long_name)
  
  # find the histology associated with the cancer group of interest 
  cohort_df_each <- histology_df  %>% 
    dplyr::filter(cancer_group %in% cancer_group_long) %>% 
    dplyr::select(Kids_First_Biospecimen_ID, harmonized_diagnosis) %>%
    dplyr::mutate(short_name = short_name_interest)
  
  #### Handle gene expression matrix -----------------------------------
  # first filter to our cohort of interest 
  expression_of_interest <- expression_data[,cohort_df_each$Kids_First_Biospecimen_ID]
  
  # annotate high vs. low expression group to the samples in the cohort 
  expression_of_gi <- expression_of_interest[rownames(expression_of_interest) == gene_of_interest, ] %>%
    t() %>% as.data.frame()
  colnames(expression_of_gi) <- "gene_of_interest"
  
  # calculate the quantile and assign groups
  upper_quantile <- quantile(expression_of_gi$gene_of_interest, (1-(as.numeric(quantile_interest)/100)))
  lower_quantile <- quantile(expression_of_gi$gene_of_interest, (as.numeric(quantile_interest)/100))
  
  bs_id_quantile_df <- expression_of_gi %>%
    tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>% 
    mutate(group = case_when(
      gene_of_interest >= upper_quantile ~"high", 
      gene_of_interest <= lower_quantile ~"low",
      TRUE ~ "middle"
    )) %>% 
    filter(group != "middle") %>%
    dplyr::select(Kids_First_Biospecimen_ID, group) %>% 
    dplyr::left_join(cohort_df_each) %>%
    dplyr::arrange(Kids_First_Biospecimen_ID)
  
  # subset to gene and sample of interest          
  expression_of_interest <- expression_data %>% 
    dplyr::select(bs_id_quantile_df$Kids_First_Biospecimen_ID) 
  
  #### Calculate geometric mean for all markers as well
  pathway_list <- pathway_df$genes %>%
    str_split(",") 
  names(pathway_list) <- pathway_df$pathway
  
  # get markers
  markers <- pathway_df$genes %>%
    str_split(",") %>% unlist() %>% unique()
  
  # filter to only markers
  expression_markers <- expression_of_interest[rownames(expression_of_interest) %in% markers,]
  # calculate z score for next steps
  expression_markers_log2_matrix <- as.matrix( log2(expression_markers + 1) )
  
  # We then calculate the Gaussian-distributed scores
  ssgsea_scores_each <- GSVA::gsva(expression_markers_log2_matrix,
                                   pathway_list,
                                   method = "ssgsea",
                                   min.sz=1, ## has to use 1 since this is a manual list 
                                   max.sz=500,## Arguments from OMPARE
                                   parallel.sz = 8, # For the bigger dataset, this ensures this won't crash due to memory problems
                                   mx.diff = TRUE,
                                   BPPARAM=SerialParam(progressbar=T))        ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)

  
  # arrange the expression matrix to match design matrix
  ssgsea_scores_each <- ssgsea_scores_each %>%
    as.data.frame() %>%
    dplyr::select(bs_id_quantile_df$Kids_First_Biospecimen_ID)
  
  bs_id_quantile_df$group <- as.factor(bs_id_quantile_df$group)
  bs_id_quantile_df$group <- relevel(bs_id_quantile_df$group, "low")
  
  # build model matrix
  mod <- model.matrix(~ bs_id_quantile_df$group)
  fit <- lmFit(as.matrix(ssgsea_scores_each), mod)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=2, n=Inf) 

  # generate results with directions
  ssgsea_results <- tt %>% 
    dplyr::mutate(direction = ifelse((logFC>0), "up", "down")) %>%
    tibble::rownames_to_column("pathway") %>%
    dplyr::select(pathway, logFC, P.Value, adj.P.Val, direction) %>%
    readr::write_tsv(file.path(results_dir, paste0(short_name_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_of_interest, "_ssgsea_pval.tsv" )))
  
  #### prepare the dataframe for plots
  ssgsea_results_plot <- ssgsea_results %>% 
    mutate(log_score = (-1)*log10(P.Value)) %>%
    arrange(log_score) %>%
    slice_head(n=25)
  ssgsea_results_plot$pathway <- factor(ssgsea_results_plot$pathway, levels = ssgsea_results_plot$pathway)
  
  # plots barplot
  pdf(file = file.path(plots_dir, paste0(short_name_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_of_interest, "_pathway_barplot.pdf" )))
  p <- ggplot(ssgsea_results_plot, aes(x = pathway, 
                                       y = log_score,
                                       fill = direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() +
    xlab("") + 
    ylab("-log10 P-Value") + 
    scale_fill_manual(name = "Direction", values = c("down" = "forest green", "up" = "red")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
    ggtitle(paste0("Diff. Expr. Pathways by ssGSEA in \n", short_name_interest, " parsed by ", quantile_interest, " quantile \n", gene_of_interest))
  
  print(p)
  dev.off()
  
  ###### prepare the scores for output
  ssgsea_results <- ssgsea_results %>% 
    dplyr::mutate(cancer_group = short_name_interest) %>%
    dplyr::mutate(gene_parsed_by = gene_of_interest) %>% 
    dplyr::mutate(percentile = quantile_interest) 
  
  combined_results <- bind_rows(combined_results, ssgsea_results)
  
  ### Clean scoring into tidy format
  ssgsea_scores_each_df <- as.data.frame(ssgsea_scores_each) %>%
    rownames_to_column(var = "pathway")
  
  #first/last_bs needed for use in gather (we are not on tidyr1.0)
  first_bs <- head(colnames(ssgsea_scores_each), n=1)
  last_bs  <- tail(colnames(ssgsea_scores_each), n=1)
  
  ssgsea_scores_each_df_tidy <- ssgsea_scores_each_df %>%
    tidyr::gather(Kids_First_Biospecimen_ID, ssgsea_score, !!first_bs : !!last_bs) %>%
    dplyr::select(Kids_First_Biospecimen_ID, pathway, ssgsea_score)  %>%
    dplyr::mutate(cancer_group = short_name_interest) %>%
    dplyr::mutate(gene_parsed_by = gene_of_interest)
  
  # add group number fo the final table as well
  ssgsea_scores_each_df_tidy <- ssgsea_scores_each_df_tidy %>%
    dplyr::left_join(bs_id_quantile_df)
  # write out individual scores
  ssgsea_scores_each_df_tidy %>%
    arrange(ssgsea_score, descending = F) %>%
    readr::write_tsv(file.path(results_dir, paste0(short_name_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_of_interest, "_ssgsea_scores.tsv")))
  
  # merge into a combined file
  ssgsea_scores_df_tidy <-  bind_rows(ssgsea_scores_df_tidy , ssgsea_scores_each_df_tidy)
}

# # write out results
# readr::write_tsv(ssgsea_scores_df_tidy, outfile_score)
# # write out combined results
# combined_results %>% readr::write_tsv(outfile_pval)


