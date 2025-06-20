library(here)
library(cmdstanr)
library(tidyverse)
library(loo)



null_01_log_lik <- readRDS(file = here("Summer_2025", "log_lik",
                                       "latent_null_01_log_lik.RDS"))
null_02_log_lik <- readRDS(file = here("Summer_2025", "log_lik",
                                       "latent_null_02_log_lik.RDS"))
zheng_06_log_lik <- readRDS(file = here("Summer_2025", "log_lik",
                                        "latent_zheng_06_log_lik.RDS"))


## process the 2015 fits
mccormick_15_log_lik <- matrix(NA, nrow = 4000, ncol = 15000)
num_data <- 15000
sample_size <- num_data/10

for(k in 1:10){
  curr_cv <- readRDS(file = here("Summer_2025", "log_lik",
                                 paste0("latent_2015_log_lik_", k, ".RDS")))
  
  mccormick_15_log_lik[, (sample_size*(k-1)+1):(k * sample_size)] <- curr_cv
}



(elpd_kfold_null_1 <- elpd(null_01_log_lik))
(elpd_kfold_null_2 <- elpd(null_02_log_lik))
(elpd_kfold_zheng <- elpd(zheng_06_log_lik))
(elpd_kfold_mccormick_15 <- elpd(mccormick_15_log_lik))

loo_compare(elpd_kfold_null_1, elpd_kfold_null_2, 
            elpd_kfold_zheng, #elpd_kfold_mccormick,
            elpd_kfold_mccormick_15)
