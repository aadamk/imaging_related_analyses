#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Rscript -e "rmarkdown::render('01-TPM_vs_harmonized_diagnosis.Rmd', clean = TRUE)"
# 
# Rscript -e "rmarkdown::render('02-lgg_cns_gtex.Rmd', clean = TRUE)"

# Rscript -e "rmarkdown::render('03a-km_w_logrank_survival_lgg.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('03b-km_w_logrank_survival_hgg.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('03c-km_w_logrank_survival_medullo.Rmd', clean = TRUE)"
# 
# Rscript -e "rmarkdown::render('04a-cox_reg_survival_lgg.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('04b-cox_reg_survival_hgg.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('04c-cox_reg_survival_medullo.Rmd', clean = TRUE)"
