#####  May 21 2025
##### K-fold CV for Null Models for ARD

library(rotasym)
library(tidyverse)
library(cmdstanr)
library(here)
library(bayesplot)
library(posterior)
library(loo)
options(mc.cores = parallel::detectCores())

set.seed(100)

data <- readRDS(here("Summer_2025", "ard_latent_mod.RDS"))
stan_data <- list(N = data$n_sample, 
                  K = data$n_subpop,
                  y = data$y_sim,
                  n_known = length(data$known_pops),
                  idx = data$known_pops,
                  p = 3,
                  known_prev = sum(data$true_subpops[data$G1_ind]/data$n_population))


# stan_data <- readRDS(here("stan_models", "stan_data_2015.RDS"))
y_sim <- stan_data$y
G1_ind <- stan_data$idx
Pg1 <- stan_data$known_prev
# y_folds <- kfold_split_random(K = 10, N = nrow(y_sim) * ncol(y_sim))
# saveRDS(y_folds, file = here("Summer_2025", "log_lik", "folds.RDS"))
y_folds <- readRDS(file = here("Summer_2025", "log_lik", "latent_sim_folds.RDS"))



# Then for 20010 Model ------------------------------------------------------

stan_file_mccormick <- here("stan_models", "mc_cormick_et_al_2010_cv.stan")
mod_mccormick_2010 <- cmdstan_model(stan_file_mccormick)
mccormick_2010_test <- cmdstan_model(here("stan_models/", "mc_cormick_et_al_2010_gq.stan"))

num_data <- 15000
sample_size <- num_data/10
mccormick_10_log_lik <- matrix(NA, nrow = 4000, ncol = 15000)

num_egos <- 6
num_alters <- 6
n_sample <- data$n_sample
sample_egos <- sample(1:num_egos, size = n_sample, replace = TRUE)
n_population <- data$n_population
n_subpop <- data$n_subpop
# true_subpop_size <- readRDS(here("Summer_2025", "log_lik", "true_subpop.RDS"))
true_subpop_size <- data$true_subpops

## could potentially use this beta, assuming each supop 
## equally spread across the ego groups, will give the right subpop proportions
beta_sim <- matrix(NA,
                   nrow = num_alters,
                   ncol = n_subpop)
for(i in 1:n_subpop){
  beta_sim[, i] <- rep(true_subpop_size[i]/(n_population), num_egos)
}




for(k in 1:10){
  test_idx  <- which(y_folds == k)
  train_idx <- which(y_folds != k)
  train_y <- y_sim
  train_y[test_idx] <- NA
  test_y  <- matrix(NA, nrow = nrow(y_sim), ncol = ncol(y_sim))
  test_y[test_idx] <- y_sim[test_idx]
  
  obs_mask_int <- (!is.na(train_y)) + 0L
  dummy <- 9999 
  train_y_stan <- train_y
  train_y_stan[is.na(train_y)] <- dummy
  test_y_stan <- test_y
  test_y_stan[is.na(test_y)] <- dummy
  
  stan_data_train <- list(E = num_egos,
                          A = num_alters,
                          N = nrow(train_y),
                          K = ncol(train_y),
                          y = train_y_stan,
                          ego = sample_egos,
                          Beta = beta_sim,
                          theta_d = c(6.2, .5),
                          theta_o = c(3, 2),
                          alpha = rep(1, num_alters)/num_alters,
                          p = 1,
                          obs_mask = obs_mask_int)
  stan_data_test <- list(E = num_egos,
                         A = num_alters,
                         N = nrow(test_y),
                         K = ncol(test_y),
                         y = test_y_stan,
                         ego = sample_egos,
                         Beta = beta_sim,
                         theta_d = c(6.2, .5),
                         theta_o = c(3, 2),
                         alpha = rep(1, num_alters)/num_alters,
                         p = 1,
                         obs_mask = 1 - obs_mask_int)
  
  fit_train <- mod_mccormick_2010$sample(data = stan_data_train,
                                     seed = 123,
                                     chains = 4,
                                     iter_sampling = 1000,
                                     iter_warmup = 1000,
                                     parallel_chains = 4,
                                     refresh = 0, sig_figs = 10)
  
  fit_gq <- mccormick_2010_test$generate_quantities(fit_train,
                                                data = stan_data_test,
                                                seed = 123)
  ll_test <- fit_gq$draws(variables = "log_lik") |> 
    as_draws_matrix() 
  
  good_col_idx <- which(!apply(ll_test, 2, function(x) any(is.nan(x))))
  log_lik_draws    <- ll_test[, good_col_idx, drop = FALSE] 
  mccormick_10_log_lik[, (sample_size*(k-1)+1):(k * sample_size)] <- log_lik_draws
  ## if I just combine these stacked column after column, they won't be in the
  ## same order as originally but will be consistent across all the folds
  ## (if I use the same folds)
}


saveRDS(mccormick_10_log_lik, file = here("Summer_2025", "log_lik",
                                      "latent_mccormick_10_log_lik.RDS"))
