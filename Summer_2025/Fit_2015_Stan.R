library(rotasym)
library(tidyverse)
#library(tidymodels) # to use coord_obs_pred
library(cmdstanr)
library(here)
#library(bayesplot)
#library(posterior)
#library(grid)
#library(gridExtra)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())



stan_data <- readRDS(here("stan_models", "stan_data_2015.RDS"))

mod_2015 <- cmdstan_model(stan_file = stan_file_2015)

stan_fit_2015 <- mod_2015$sample(data = stan_data,
                                 seed = 123,
                                 chains = 4,
                                 iter_sampling = 1000,
                                 iter_warmup = 1000,
                                 parallel_chains = 4,
                                 refresh = 100)

stan_fit_2015$save_object(file = here("stan_models", "2015_cluster_fit.RDS"))
