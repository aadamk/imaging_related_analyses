# Author: Run Jin
# Output one representative tree using the party package
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("party")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"),
  make_option(c("-e","--path_to_meta"),type="character",
              help="path to the file directory that contains combined meta and gene expression data")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cg_list <-unlist(strsplit(opt$cg_interest,","))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "relative_importance")

plots_dir <- file.path(analysis_dir, "plots", "rf_party")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

meta_exp_input_dir <- opt$path_to_meta

# run the tree generation function for all cancer groups 
for (i in 1:length(cg_list)){
  x <- cg_list[i]
  
  # read in the combined data with metadata + gene expression 
  meta_exp_combined <- readr::read_tsv(file.path(meta_exp_input_dir, 
                                                 paste0("meta_exp_combined_data_in_", x, ".tsv")))
  
  # get the gene variables 
  gene_variables <- meta_exp_combined[, 7:37] %>% colnames()
  
  ########### prepare the formula and data based on which cancer group we are looking at
  if(x == "LGG"){
    # piece together a model
    survival.formula <- paste0("Surv(PFS_days, PFS_status_recode)", " ~ ",
                                paste0(gene_variables, collapse = " + "),
                                sep = "")
    # make the combined data compatible with the analysis 
    combined_data_sub <- meta_exp_combined %>% 
      dplyr::select(-c("Kids_First_Biospecimen_ID", "Kids_First_Participant_ID", 
                       "OS_status_recode", "OS_days", "PFS_status")) %>% 
      dplyr::filter(!is.na(PFS_days))
  
  } else{
    # piece together a model
    survival.formula <- paste0("Surv(OS_days, OS_status_recode)", " ~ ",
                                paste0(gene_variables, collapse = " + "),
                                sep = "")
    # make the combined data compatible with the analysis 
    combined_data_sub <- meta_exp_combined %>% 
      dplyr::select(-c("Kids_First_Biospecimen_ID", "Kids_First_Participant_ID", 
                       "PFS_status_recode", "PFS_days", "OS_status")) %>% 
      dplyr::filter(!is.na(OS_days))
  }
    
  ####### run the analysis with optimal parameter
  ctree_ctrl_obj <- ctree_control(teststat = "max",
                                  testtype = "Teststatistic",
                                  mincriterion = 0.95,
                                  minbucket = 5,
                                  stump = FALSE,
                                  maxsurrogate = 0,
                                  mtry = 0, 
                                  savesplitstats = TRUE, 
                                  maxdepth = 0, 
                                  remove_weights = FALSE)
  
  # get the results
  ctree_results <- ctree(as.formula(survival.formula), 
                         combined_data_sub, 
                         subset = NULL,
                         weights = NULL,
                         controls = ctree_ctrl_obj, 
                         xtrafo = ptrafo, 
                         ytrafo = ptrafo,
                         scores = NULL)
  
  # plot the tree and save the plot 
  pdf(file.path(plots_dir, paste0("rf_party_plot_in_", x, ".pdf")))
  plot(ctree_results)
  dev.off()

}
