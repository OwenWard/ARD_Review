library(here)
code_path <- here("Sahai_et_al_2018/Swupnil_Code/")

source(paste0(code_path,"src/load_data/surveys.R"))
source(paste0(code_path,"src/load_data/names.R"))

## this gives warnings but think works ok
source(paste0(code_path,"src/load_data/occs.R"))

source(paste0(code_path, "src/funcs/mix_kernel.R"))

## could change these to cmdstanr? 
## would need to change some of the code to process the results also
source(paste0(code_path, "src/load_data/compile_models.R"))

## made a small change in this for the last plot
## need to check this again
## somewhat slow to run
source(paste0(code_path, "src/fit/kernel_occ.R"))
source(paste0(code_path, "src/fit/kernel_name.R"))
## this one giving some warnings

source(paste0(code_path, "src/fit/kernel_spline_name.R"))

## error with this one, g_k_name not found
source(paste0(code_path, "src/fit/kernel_spline_comb.R"))

source(paste0(code_path, "src/fit/kernel_spline_occ.R"))

source(paste0(code_path, "src/fit/without_old/kernel_spline_name_woOld.R"))

## this has same issue as above
source("src/fit/without_old/kernel_spline_comb_woOld.R");


source(paste0(code_path, "src/fit/without_old/kernel_spline_occ_woOld.R"))
