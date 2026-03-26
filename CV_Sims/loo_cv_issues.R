library(tidyverse)
library(cmdstanr)
library(here)
library(bayesplot)
library(posterior)
options(mc.cores = parallel::detectCores())
library(loo)
library(matrixStats)




# Simulate ARD Like Data ------------------------------------------
set.seed(100)
I <- 100   # number of rows
J <- 10    # number of columns
N <- I * J

sigma_a_true <- 5
sigma_b_true <- 0.3

a_true <- rnorm(I, 0, sigma_a_true)
b_true <- rnorm(J, 0, sigma_b_true)

dat <- expand.grid(row = 1:I, col = 1:J)
dat <- dat[order(dat$row, dat$col), ]
for(k in 1:nrow(dat)) {
  dat$y[k] <- rpois(n = 1, lambda = exp(a_true[dat$row[k]] + b_true[dat$col[k]]))
}

# dat$y <- rpois(n = 1, lambda = exp(a_true + b_true))
# quick look
head(dat)

# matrix form for convenience later
Y_mat <- matrix(dat$y, nrow = I, ncol = J, byrow = TRUE)

stan_data <- list(
  n_i = I,
  n_k = J,
  y = Y_mat
)

mod_od <- cmdstan_model("Revision_March_2026/overdispersed_nsum.stan")

fit_od <- mod_od$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 200
)

loo_od <- fit_od$loo(cores = 4)
loo_od

####


## Fit actual Poisson model also with a fixed value of
## the sigma_d parameter

stan_data_pois <- list(
  n_i = I,
  n_k = J,
  y = Y_mat,
  sigma_alpha = 10 ## and 1,5, 10, 20, 50
  # alphas = rep(0, I)
  # alphas = a_true
)

mod_poisson <- cmdstan_model("Revision_March_2026/poisson.stan")

fit_pois <- mod_poisson$sample(
  data = stan_data_pois,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 200
)

# loo_pois_sigma_5 <- loo_pois

loo_pois_sigma_10 <- fit_pois$loo(cores = 4)
loo_pois_sigma_10 ## 

loo_compare(loo_pois_sigma_1, loo_pois_sigma_5, loo_pois_sigma_10,
            loo_pois_sigma_20,
            loo_pois_sigma_50)

## this selects the wrong model, not reliable for model selection here

## similarly, if we look at the draws of the likelihood for any ARD
## entry as we vary sigma_alpha, they can dramatically vary

pois_draws <- as_draws_df(fit_pois$draws())

sig_al <- stan_data_pois$sigma_alpha

pois_draws |> select(starts_with("log_lik")) |> 
  select(`log_lik[1,1]`) |> 
  rename(log_lik  = `log_lik[1,1]`) |> 
  ggplot(aes(log_lik)) +
  geom_histogram() +
  labs(title = paste("Sigma alpha = ", sig_al))





# Fit the Park Models -----------------------------------------------------

## here we fit the park models to this simulated data 
## making small modifications to the code needed for newer versions
## of stan. We see that all of these give large values of the k
## diagnostic

### trying to check the park model for this
mod_park <- cmdstan_model("Revision_March_2026/park_poisson.stan")
mod_park <- cmdstan_model("Revision_March_2026/rand_mix_park.stan")
mod_park <- cmdstan_model("Revision_March_2026/park_social.stan")

park_data <- list(
  N = nrow(dat),
  I = I,
  J = J,
  ii = dat$row,
  jj = dat$col,
  y = dat$y,
  pops = rep(0.1, J),
  D = 1,
  sds = rep(1, J)
  # 
  # sigma_d = 0.50 ## and 10, 5, 1, 0.1, 0.001
  # alphas = rep(0, I)
  # alphas = a_true
)


fit_park <- mod_park$sample(
  data = park_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 200
)

loo_park <- fit_park$loo(cores = 4)
loo_park
