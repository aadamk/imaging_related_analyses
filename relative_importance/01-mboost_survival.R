# Author: Run Jin
# Perform mboost survival analysis to measure the relative importance of SLC1A5 
# relative to known amino acid transporters with respect to disease prognosis
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("mboost")
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-c","--count"),type="character",
              help="gene count data from OpenPedCan RNA samples (.rds) "),
  make_option(c("-l","--cg_interest"),type="character",
              help="comma separated list of cancer groups of interest"),
  make_option(c("-m","--short_long_match"),type="character",
              help="match between long and short names (.tsv)")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
cg_list <-unlist(strsplit(opt$cg_interest,","))

#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "relative_importance")
results_dir <- file.path(analysis_dir, "results")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

#### Read in files necessary for analyses --------------------------------------
# histology file
histology_df <- readr::read_tsv(opt$histology, guess_max=100000)

# gene expected count file
count_matrix <- readRDS(opt$count)

# file matching short and long histology
short_long_match <- readr::read_tsv(opt$short_long_match)


