#### October 15th 2024

## In this example we will simulate data where
## there is a kernel function driving the mixing probabilities
## and we will fit a range of models to this

## Null Models
## 2006 Model
## 2010 Model
## 2019 Model

## we will compare each of these in terms of 
## estimated degree 
## posterior predictive checks


# Setup and Load Packages -------------------------------------------------

library(rotasym)
library(tidyverse)
library(tidymodels) # to use coord_obs_pred
library(cmdstanr)
library(here)
library(bayesplot)
library(posterior)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())

source(here("Summer_2024/", "helper_model_checking.R"))


# Simulate ARD ------------------------------------------------------------


sim_pars <- readRDS(here("Summer_2024", "complete_mixing_sim.RDS"))

n <- sim_pars$n
num_subpop <- sim_pars$k

true_degree <- sim_pars$degree

simulate_mixing <- function(d, omega, M, beta, ego) {
  y <- matrix(nrow = length(d), ncol = ncol(beta))
  
  for (i in 1:nrow(y)) {
    for (j in 1:ncol(y)) {
      mu_ij <- d[i] * M[ego[i], ] %*% beta[, j]
      y[i, j] <- MASS::rnegbin(1, mu_ij, omega[j] * mu_ij)
    }
  }
  
  return(y)
}


y_sim <- simulate_mixing(d = sim_pars$degree,
                         omega = sim_pars$omega,
                         M = sim_pars$M,
                         beta = sim_pars$beta,
                         ego = sim_pars$ego)

## here each subpopulation is a name

# Fit Null Models ---------------------------------------------------------

stan_data_null <- list(N = nrow(y_sim),
                       K = ncol(y_sim),
                       y = y_sim)

stan_file_null_01 <- here("Summer_2024/", "null_model_01.stan")

mod_null_01 <- cmdstan_model(stan_file_null_01)

stan_fit_null_01 <- mod_null_01$sample(data = stan_data_null,
                                       seed = 123,
                                       chains = 4,
                                       iter_sampling = 1000,
                                       iter_warmup = 1000,
                                       parallel_chains = 4,
                                       refresh = 100)

stan_fit_null_01$summary(variables = c("log_d", "beta"))

## look at the estimated degree

est_degrees_null_01 <- stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("log_d")) |> 
  mutate(draw = row_number())  |> 
  mutate(degree = exp(log_d))


true_degrees <- tibble(node = 1:nrow(y_sim),
                       true_degree = samp_degree)

true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")

est_degrees_null_01 |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  labs(x = "Expected Degree", title = "Posterior Estimates",
       subtitle = "Null Erdos Renyi Model", y = "") 


## look at simple ppc
ppc_null_1 <- construct_ppc(stan_fit_null_01, y_sim)

ppc_fit_null_1 <- plot_ests(ppc_null_1$ppc_draws, 
                            ppc_null_1$y_tibble,
                            prop_val = 5)
ppc_fit_null_1 + 
  labs(title = "Null Erdos Renyi Model") 

## Then fit the second null model
stan_file_null_02 <- here("Summer_2024/", "null_model_02.stan")

mod_null_02 <- cmdstan_model(stan_file_null_02)

stan_fit_null_02 <- mod_null_02$sample(data = stan_data_null,
                                       seed = 123,
                                       chains = 4,
                                       iter_sampling = 1000,
                                       iter_warmup = 1000,
                                       parallel_chains = 4,
                                       refresh = 100)

stan_fit_null_02$summary(variables = c("beta", "log_d[1]", "log_d[2]"))

est_degrees_null_02 <- stan_fit_null_02$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("log_d")) |> 
  mutate(draw = row_number())  |>
  pivot_longer(cols = starts_with("log_d"),
               names_to = "node", 
               values_to = "log_d") |> 
  mutate(degree = exp(log_d),
         node_id = parse_number(node)) 

est_degrees_null_02 |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  labs(title = "Estimated Degree Distribution",
       subtitle = "Null Model, varying Degrees")

true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")

ppc_null_2 <- construct_ppc(stan_fit_null_02, y_sim)

ppc_fit_null_2 <- plot_ests(ppc_null_2$ppc_draws, 
                            ppc_null_2$y_tibble,
                            prop_val = 1)
ppc_fit_null_2 + 
  labs(title = "Null Model, varying Degrees") 


# Fit 2006 Model ----------------------------------------------------------

## to fit this need to specify the beta priors here
## can do because we know the true degrees

var0 <- apply(y_sim, 1, var)
(var0 <- which(var0 == 0))
if (length(var0) > 0) {
  y_sim <- y_sim[-var0,]
  samp_degree_2006 <- samp_degree[-var0]
  I <- nrow(y_sim)
}



## need prior for the beta_k's here
## estimate it as
## this model doesn't account for the distance 
mu_beta <- log(apply(y_sim, 2, sum)/sum(true_degree))
## no way to specify sd_beta
sd_beta <- rep(1, length(mu_beta))

stan_data_2006 <- list(N = nrow(y_sim),
                       K = ncol(y_sim),
                       y = y_sim,
                       mu_beta = mu_beta,
                       sigma_beta = sd_beta)

stan_file_2006 <- here("Summer_2024", "zheng_et_al_2006.stan")

mod_2006 <- cmdstan_model(stan_file_2006)

stan_fit_2006 <- mod_2006$sample(data = stan_data_2006,
                                 seed = 123,
                                 chains = 4,
                                 iter_sampling = 1000,
                                 iter_warmup = 1000,
                                 parallel_chains = 4,
                                 refresh = 100)

stan_fit_2006$summary(variables = c("mu_alpha", 
                                    "sigma_alpha",
                                    "beta"))

mcmc_trace(stan_fit_2006$draws(variables = c("mu_alpha",
                                             "sigma_alpha",
                                             "beta")))

est_degrees_2006 <- stan_fit_2006$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("alpha")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("alpha"),
               names_to = "node",
               values_to = "est") |> 
  mutate(degree = exp(est),
         node = parse_number(node)) 


true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")

est_degrees_2006 |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  labs(x = "Expected Degree", title = "Posterior Estimates",
       subtitle = "2006 Model", y = "") 

ppc_2006 <- construct_ppc(stan_fit_2006, y_sim)

ppc_fit_2006 <- plot_ests(ppc_2006$ppc_draws,
                          ppc_2006$y_tibble,
                          prop_val = 1)
ppc_fit_2006 + 
  labs(title = "Zheng et al, 2006") 

# Fit 2010 Model ----------------------------------------------------------

## this will be the correct model here so should do reasonably well
## also know the true egos and the beta here

num_egos <- 6
num_alters <- 8

mix_data_stan <- list(E = num_egos,
                      A = num_alters,
                      K = ncol(y_sim),
                      N = nrow(y_sim),
                      y = y_sim,
                      ego = sim_pars$ego,
                      Beta = sim_pars$beta,
                      theta_d = c(6.2, .5),
                      theta_o = c(3, 2),
                      alpha = rep(1, num_alters)/num_alters,
                      p = 1)

stan_file_2010 <- here("Summer_2024", "mccormick_et_al_2010.stan")
mod_2010 <- cmdstan_model(stan_file = stan_file_2010)

stan_fit_2010 <- mod_2010$sample(data = mix_data_stan,
                                 seed = 123,
                                 chains = 4,
                                 iter_sampling = 1000,
                                 iter_warmup = 1000,
                                 parallel_chains = 4,
                                 refresh = 100)

## then do model checking and ppc for this also
stan_fit_2010$summary(variables = c("M", "omega"))

## plot the estimated degree distributions from this
est_degrees_2010 <- stan_fit_2010$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("log_d")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("log_d"),
               names_to = "node",
               values_to = "log_estimate") |> 
  mutate(node = parse_number(node)) |> 
  mutate(est_degree = exp(log_estimate)) 

est_degrees_2010 |> 
  ggplot(aes(est_degree)) +
  geom_histogram() +
  labs(title = "Posterior Estimates",
       subtitle = "2010 Model")

## compare to true distribution, seems to over estimate
true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")

## then do ppc for this model also
ppc_2010 <- construct_ppc(stan_fit_2010, y_sim)

ppc_fit_2010 <- plot_ests(ppc_2010$ppc_draws,
                          ppc_2010$y_tibble,
                          prop_val = 3)
ppc_fit_2010 + 
  labs(title = "McCormick et al, 2010") 



# Fit 2019 Model ----------------------------------------------------------


## just need to figure out how to adjust the kernel to use for this
## where we don't have to assume the subpops are names
