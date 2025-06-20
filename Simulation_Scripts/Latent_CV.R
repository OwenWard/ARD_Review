### May 21 - Do K fold CV for Latent Space Model


library(tidyverse)
library(cmdstanr)
library(here)
library(bayesplot)
library(posterior)
##library(grid)
##library(gridExtra)
options(mc.cores = parallel::detectCores())

set.seed(100)

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
fold_id <- jobid


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
y_folds <- readRDS(file = here("Summer_2025", "log_lik", "latent_sim_folds.RDS"))


## maybe save these folds to ensure use the same every time


stan_data <- readRDS(here("stan_models", "stan_data_2015.RDS"))
stan_file_2015 <- here("stan_models", "mc_cormick_and_zheng_2015_cv.stan")
mod_2015 <- cmdstan_model(stan_file = stan_file_2015)
mod_2015_test <- cmdstan_model(stan_file = here("stan_models",
                                                "mc_cormick_and_zheng_2015_gq.stan"))


k <- fold_id
# for(k in 1:10){
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

stan_data_train <- list(N = nrow(train_y),
                        K = ncol(train_y),
                        p = 3,
                        y = train_y_stan,
                        n_known = length(G1_ind),
                        idx = G1_ind,
                        known_prev = Pg1,
                        obs_mask = obs_mask_int)
stan_data_test <- list(N = nrow(test_y),
                       K = ncol(test_y),
                       p = 3,
                       y = test_y_stan,
                       n_known = length(G1_ind),
                       idx = G1_ind,
                       known_prev = Pg1,
                       obs_mask = 1 - obs_mask_int)

fit_train <- mod_2015$sample(data = stan_data_train,
                                seed = 123,
                                chains = 4,
                                iter_sampling = 1000,
                                iter_warmup = 1000,
                                parallel_chains = 4,
                                sig_figs = 15,
                                refresh = 100)
fit_train$save_object(file = here("Summer_2025", "log_lik",
                                  paste0("latent_2015_fit_", k, ".RDS")))

fit_gq <- mod_2015_test$generate_quantities(fit_train,
                                           data = stan_data_test,
                                           sig_figs = 15,
                                           seed = 123)

ll_test <- fit_gq$draws(variables = "log_lik") |> 
  as_draws_matrix() 

good_col_idx <- which(!apply(ll_test, 2, function(x) any(is.nan(x))))
log_lik_draws    <- ll_test[, good_col_idx, drop = FALSE]   # subset matrix
#   null_01_log_lik[, (sample_size*(k-1)+1):(k * sample_size)] <- log_lik_draws
#   ## if I just combine these stacked column after column, they won't be in the
#   ## same order as originally but will be consistent across all the folds
#   ## (if I use the same folds)
# }


saveRDS(log_lik_draws, file = here("Summer_2025", "log_lik",
                                     paste0("latent_2015_log_lik_", k, ".RDS")))
