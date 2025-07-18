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


data <- readRDS(here("Summer_2025", "ard_latent_mod.RDS"))
stan_data <- list(N = data$n_sample, 
                  K = data$n_subpop,
                  y = data$y_sim,
                  n_known = length(data$known_pops),
                  idx = data$known_pops,
                  p = 3,
                  known_prev = sum(data$true_subpops[data$G1_ind]/data$n_population))

stan_file_2015 <- here("stan_models", "mc_cormick_and_zheng_2015.stan")
mod_2015 <- cmdstan_model(stan_file = stan_file_2015)

stan_fit_2015 <- mod_2015$sample(data = stan_data,
                                 seed = 123,
                                 chains = 4,
                                 iter_sampling = 1000,
                                 iter_warmup = 1000,
                                 parallel_chains = 4,
                                 refresh = 100)

stan_fit_2015$save_object(file = here("stan_models",
                                      "ard_latent_2015_cluster_fit.RDS"))

# 
# fit <- readRDS(here("stan_models", "2015_cluster_fit.RDS"))
# # 
# fit$summary(variables = c("eta"))
# 
# mcmc_trace(fit$draws(variables = "eta"))
# 
# 
# ppc_2015 <- construct_ppc(fit, y_sim)
# 
# ppc_fit_2015 <- plot_ests(ppc_2015$ppc_draws,
#                           ppc_2015$y_tibble,
#                           prop_val = 5)
# ppc_fit_2015 + 
#   labs(title = "McCormick + Zheng, 2015") 
# 
# 
# 
# loo_est_2015 <- fit$loo(cores = 4)
# loo_est_2015
# 
# loo_compare(loo_est, loo_est_2015)

