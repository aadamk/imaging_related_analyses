# Author: Run Jin
# Barplot function Adapted from OMPARE
# Script to plot significant altered pathways from various pathway analysis methods

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("stringr"))

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-s","--ssgsea"),type="character",
              help="combined results file from ssGSEA analysis (.tsv)"),
  make_option(c("-g","--gsnca"),type="character",
              help="combined results file from GSNCA analysis (.tsv)"),
  make_option(c("-d","--dgca"),type="character",
              help="combined results file from DGCA analysis (.tsv)"),
  make_option(c("-l","--cg_gene_interest"),type="character",
              help="file containing gene of interest and matching cancer group (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))

#### function to create drug barplots -------------------------------------------
pathway_barplots <- function(dat, xlab, ylab, top = 20, title){
  dat <- dat %>%
    dplyr::select(xlab, ylab) %>%
    unique() %>%
    filter(get(ylab) != 0) %>%
    arrange(get(ylab), descending = FALSE) %>%
    slice_head(n = top) %>%
    as.data.frame()
  
  dat <- dat %>% 
    mutate(log_score = (-1)*log10(get(ylab))) %>%
    arrange(log_score, descending = TRUE)
  dat[,xlab] <- factor(dat[,xlab], levels = unique(dat[,xlab]))
  
  p <- ggplot(dat, aes(x = get(xlab), 
                       y = log_score,
                       fill = log_score)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() +
    xlab("") + 
    ylab("-log10 Adj. P-Value") + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
    ggtitle(title) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
    guides(fill = "none")

  return(p)
}

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "coexpression_gene_pathway_analysis")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

####### read in files
# pathway analysis results
ssgsea_results <- readr::read_tsv(opt$ssgsea)
gsnca_results <- readr::read_tsv(opt$gsnca)
dgca_results <- readr::read_tsv(opt$dgca)

# file containing genes of interest in each cancer group
cg_gene_interest <- readr::read_tsv(opt$cg_gene_interest)

########## generate plots
for (i in 1:nrow(cg_gene_interest)){
  # find the cancer group of interest
  cg_interest <- cg_gene_interest[i,1] %>% as.character()
  # find the gene of interest
  gene_interest <- cg_gene_interest[i,2] %>% as.character()
  # get the quantile of interest 
  quantile_interest <- cg_gene_interest[i,3] %>% as.character()
  
  ############ plot ssGSEA results
  # filter ssGSEA results
  ssgsea_results_each <- ssgsea_results %>%
    dplyr::filter(cancer_group == cg_interest) %>%
    dplyr::filter(gene_parsed_by == gene_interest) %>%
    dplyr::filter(percentile == quantile_interest)

  if(length(ssgsea_results_each)>1){
    # define directory
    plots_dir_ssgsea <- file.path(plots_dir, "ssgsea_pathway_barplot")
    if(!dir.exists(plots_dir_ssgsea )) {dir.create(plots_dir_ssgsea)}

    # plots barplot
    pdf(file = file.path(plots_dir_ssgsea, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_pathway_barplot.pdf" )))
    p<-pathway_barplots(dat = ssgsea_results_each,
                  xlab = "description", ylab = "pval",
                  top = 20,
                  title = paste0("Diff. Expr. Pathways by ssGSEA in ", cg_interest, "\n parsed by ", quantile_interest, " quantile ", gene_interest))
    print(p)
    dev.off()
  }

  ############ plot GSNCA results
  # filter GSNCA results
  gsnca_results_each <- gsnca_results %>%
    dplyr::filter(cancer_group == cg_interest) %>%
    dplyr::filter(gene_parsed_by == gene_interest) %>%
    dplyr::filter(percentile == quantile_interest)

  if(length(gsnca_results_each)>1){
    # define directory
    plots_dir_gsnca <- file.path(plots_dir, "gsnca_pathway_barplot")
    if(!dir.exists(plots_dir_gsnca )) {dir.create(plots_dir_gsnca)}

    # plots barplot
    pdf(file = file.path(plots_dir_gsnca, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_pathway_barplot.pdf" )))
    p<-pathway_barplots(dat = gsnca_results_each,
                  xlab = "pathway_description", ylab = "pvalue",
                  top = 20,
                  title = paste0("Diff. Expr. Pathways by GSNCA in ", cg_interest, "\n parsed by ", quantile_interest, " quantile ", gene_interest))
    print(p)
    dev.off()
  }
  
  ############ plot DGCA results
  # filter DGCA results 
  dgca_results_each <- dgca_results %>% 
    dplyr::filter(cancer_group == cg_interest) %>%
    dplyr::filter(gene_parsed_by == gene_interest) %>%
    dplyr::filter(percentile == quantile_interest) 
  
  ontology_list <- dgca_results %>% pull(Ontology) %>% unique()
  
  for (j in 1:length(ontology_list)) {
    ontology_interest <- ontology_list[j]
    dgca_bp_gain <- dgca_results_each %>%
      filter(change_dir == "gain_of_correlation_genes") %>%
      filter(Ontology == ontology_interest)
    
    dgca_bp_loss <- dgca_results_each %>%
      filter(change_dir == "loss_of_correlation_genes") %>%
      filter(Ontology == ontology_interest)
    
    # plot for gain of correlation
    if(nrow(dgca_bp_gain)>1){
      # define directory
      plots_dir_dgca <- file.path(plots_dir, "dgca_pathway_barplot", "gain_corr")
      if(!dir.exists(plots_dir_dgca )) {dir.create(plots_dir_dgca, recursive=TRUE)}
      
      # plots barplot
      pdf(file = file.path(plots_dir_dgca, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_gain_corr_", ontology_interest, "_pathway_barplot.pdf" )))
      p<-pathway_barplots(dat = dgca_bp_gain, 
                          xlab = "Term", ylab = "Pvalue",
                          top = 20, 
                          title = paste0("Gain of Diff. Expr. Pathways by DGCA in \n", cg_interest, " parsed by ", quantile_interest, " quantile ", gene_interest))
      print(p)
      dev.off()
    }
    
    # plot for loss of correlation
    if(nrow(dgca_bp_loss)>1){
      # define directory
      plots_dir_dgca <- file.path(plots_dir, "dgca_pathway_barplot", "loss_corr")
      if(!dir.exists(plots_dir_dgca )) {dir.create(plots_dir_dgca, recursive=TRUE)}
      
      # plots barplot
      pdf(file = file.path(plots_dir_dgca, paste0(cg_interest, "_parsed_by_", quantile_interest, "_quantile_", gene_interest, "_loss_corr_", ontology_interest, "_pathway_barplot.pdf" )))
      p<-pathway_barplots(dat = dgca_bp_loss, 
                          xlab = "Term", ylab = "Pvalue",
                          top = 20, 
                          title = paste0("Loss of Diff. Expr. Pathways by DGCA in \n", cg_interest, " parsed by ", quantile_interest, " quantile ", gene_interest))
      print(p)
      dev.off()
    }
  }
}


