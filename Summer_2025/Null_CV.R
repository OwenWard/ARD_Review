#####  May 21 2025
##### K-fold CV for Null Models for ARD

library(rotasym)
library(tidyverse)
library(cmdstanr)
library(here)
library(bayesplot)
library(posterior)
library(grid)
library(gridExtra)
options(mc.cores = parallel::detectCores())




stan_data <- readRDS(here("stan_models", "stan_data_2015.RDS"))

y_sim <- stan_data$y_sim


y_folds <- kfold_split_random(K = 10, N = nrow(y_sim))


stan_file_null_01 <- here("stan_models", "null_model_01_scaled.stan")

mod_null_01 <- cmdstan_model(stan_file_null_01)

null_01_test <- cmdstan_model(here("stan_models/", "null_model_01_gq.stan"))
num_data <- 15000
sample_size <- num_data/10
null_01_log_lik <- matrix(NA, nrow = 4000, ncol = 15000)

for(k in 1:10){
  train_y <- y_sim[y_folds != k, ]
  test_y <- y_sim[y_folds == k, ]
  stan_data_train <- list(N = nrow(train_y),
                          K = ncol(train_y),
                          y = train_y,
                          n_known = length(G1_ind),
                          idx = G1_ind,
                          known_prev = Pg1)
  stan_data_test <- list(N = nrow(test_y),
                         K = ncol(test_y),
                         y = test_y,
                         n_known = length(G1_ind),
                         idx = G1_ind,
                         known_prev = Pg1)
  fit_train <- mod_null_01$sample(data = stan_data_train,
                                  seed = 123,
                                  chains = 4,
                                  iter_sampling = 1000,
                                  iter_warmup = 1000,
                                  parallel_chains = 4,
                                  refresh = 0)
  
  fit_gq <- null_01_test$generate_quantities(fit_train,
                                             data = stan_data_test,
                                             seed = 123)
  log_lik_draws <- fit_gq$draws(variables = "log_lik") |> 
    as_draws_matrix()
  null_01_log_lik[, (sample_size*(k-1)+1):(k * sample_size)] <- log_lik_draws
  ## if I just combine these stacked column after column, they won't be in the
  ## same order as originally but will be consistent across all the folds
  ## (if I use the same folds)
}


saveRDS(null_01_log_lik, file = here("Summer_2025", "log_lik",
                                     "null_01_log_lik.RDS"))




# Then Repeat for Second Null Model ---------------------------------------



stan_file_null_02 <- here("stan_models", "null_model_02_scaled.stan")

mod_null_02 <- cmdstan_model(stan_file_null_02)

null_02_test <- cmdstan_model(here("stan_models/", "null_model_02_gq.stan"))

null_02_log_lik <- matrix(NA, nrow = 4000, ncol = 15000)

for(k in 1:10){
  train_y <- y_sim[y_folds != k, ]
  test_y <- y_sim[y_folds == k, ]
  stan_data_train <- list(N = nrow(train_y),
                          K = ncol(train_y),
                          y = train_y,
                          n_known = length(G1_ind),
                          idx = G1_ind,
                          known_prev = Pg1)
  stan_data_test <- list(N = nrow(test_y),
                         K = ncol(test_y),
                         y = test_y,
                         n_known = length(G1_ind),
                         idx = G1_ind,
                         known_prev = Pg1)
  fit_train <- mod_null_02$sample(data = stan_data_train,
                                  seed = 123,
                                  chains = 4,
                                  iter_sampling = 100,
                                  iter_warmup = 100,
                                  parallel_chains = 4,
                                  refresh = 0)
  
  fit_gq <- null_02_test$generate_quantities(fit_train,
                                             data = stan_data_test,
                                             seed = 123)
  log_lik_draws <- fit_gq$draws(variables = "log_lik") |> 
    as_draws_matrix()
  null_02_log_lik[, (sample_size*(k-1)+1):(k * sample_size)] <- log_lik_draws
  ## if I just combine these stacked column after column, they won't be in the
  ## same order as originally but will be consistent across all the folds
  ## (if I use the same folds)
}




# Figuring out matrix folds -----------------------------------------------

test_idx  <- which(y_folds == k)               # linear indices in fold k
train_idx <- which(y_folds != k)

train_y <- y_sim
train_y[test_idx] <- NA
test_y  <- matrix(NA, nrow = nrow(y_sim), ncol = ncol(y_sim))
test_y[test_idx] <- y_sim[test_idx]

obs_mask_int <- (!is.na(train_y)) + 0L

dummy <- 9999 
train_y_stan <- train_y
train_y_stan[is.na(train_y)] <- dummy

dummy <- 9999 
test_y_stan <- test_y
test_y_stan[is.na(test_y)] <- dummy

stan_file_null_01 <- here("stan_models", "null_model_01_scaled_cv.stan")
mod_null_01 <- cmdstan_model(stan_file_null_01)
null_01_test <- cmdstan_model(here("stan_models/", "null_model_01_gq.stan"))

stan_data_train <- list(N = nrow(train_y),
                        K = ncol(train_y),
                        y = y_train_stan,
                        n_known = length(G1_ind),
                        idx = G1_ind,
                        known_prev = Pg1,
                        obs_mask = obs_mask_int)
stan_data_test <- list(N = nrow(test_y),
                       K = ncol(test_y),
                       y = test_y_stan,
                       n_known = length(G1_ind),
                       idx = G1_ind,
                       known_prev = Pg1,
                       obs_mask = 1 - obs_mask_int)
fit_train <- mod_null_01$sample(data = stan_data_train,
                                seed = 123,
                                chains = 4,
                                iter_sampling = 1000,
                                iter_warmup = 1000,
                                parallel_chains = 4,
                                refresh = 0)

fit_gq <- null_01_test$generate_quantities(fit_train,
                                           data = stan_data_test,
                                           seed = 123)

ll_test <- fit_gq$draws(variables = "log_lik") |> 
  as_draws_matrix() 

good_col_idx <- which(!apply(ll_test, 2, function(x) any(is.nan(x))))
good_cols    <- ll_test[, good_col_idx, drop = FALSE]   # subset matrix
