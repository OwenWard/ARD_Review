
# Setup and Load Packages -------------------------------------------------

library(tidyverse)
library(stringr)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(here)
library(scales)
library(gtools)
library(rstan)
library(tidymodels)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())
source(here("Summer_2024/", "helper_model_checking.R"))

# Simulate with Age information also --------------------------------------

data <- read_csv(here("Summer_2024/", "sim_pars", "omni.csv"))
names_data <- data[, 2:13]

gender <- as.factor(data$Gender)
weights <- data$wgt
genders <- levels(gender) 
age <- as.factor(data$Age)
ages <- levels(age)
age_numeric <- data$age

ego_sex_age <- rep(0, nrow(data))
ego_sex_age[(gender == genders[2]) & (age %in% ages[1:2])] <- 1
# Male 18-24
ego_sex_age[(gender == genders[2]) & (age %in% ages[3:10])] <- 2
# Male 25-64
ego_sex_age[(gender == genders[2]) & (age %in% ages[11:13])] <- 3
# Male 65+
ego_sex_age[(gender == genders[1]) & (age %in% ages[1:2])] <- 4
# Female 18-24
ego_sex_age[(gender == genders[1]) & (age %in% ages[3:10])] <- 5
# Female 25-64
ego_sex_age[(gender == genders[1]) & (age %in% ages[11:13])] <- 6
# Female 65+

pop_raw <- read.csv(here("Summer_2024/", "sim_pars", "pop_age_2013.csv"),
                    stringsAsFactors = TRUE)[-c(1, 102), ]

pop_age_sex <- c(sum(pop_raw$Male[2:18]), sum(pop_raw$Male[19:25]), 
                 sum(pop_raw$Male[26:65]), sum(pop_raw$Male[66:100]))
pop_age_sex <- c(pop_age_sex, 
                 sum(pop_raw$Female[2:18]),
                 sum(pop_raw$Female[19:25]), 
                 sum(pop_raw$Female[26:65]),
                 sum(pop_raw$Female[66:100]))
pop_age_sex <- t(as.matrix(pop_age_sex))
colnames(pop_age_sex) <- c("M_1-17", "M_18-24", "M_25-64", "M_65+",
                           "F_1-17", "F_18-24", "F_25-64", "F_65+")


source(here("Summer_2024", "sim_pars", "names.R"))
## this gives beta_names along with mu_k_name sigma_k_name and sum_prob_k_name

degree <- read.csv(here("Summer_2024", "sim_pars", "degree_mix.csv"))
degree_sim <- degree$CombSexAge
n_names <- ncol(names_data)
omega_sim <- runif(n_names, 0.85, 1)
mixing_sim <- readRDS(here("Summer_2024", "sim_pars", "mixing_sim.RDS"))

M_sim <- matrix(c(rdirichlet(1, mixing_sim[1, ]),
                  rdirichlet(1, mixing_sim[2, ]),
                  rdirichlet(1, mixing_sim[3, ]),
                  rdirichlet(1, mixing_sim[4, ]),
                  rdirichlet(1, mixing_sim[5, ]),
                  rdirichlet(1, mixing_sim[6, ])),
                nrow = 6, ncol = 8, byrow = TRUE)

sim_pars_age <- list(n = length(degree_sim),
                     k = ncol(names_data),
                     degree = degree_sim,
                     omega = omega_sim,
                     M = M_sim,
                     beta = beta_names,
                     ego = ego_sex_age)

data_pars <- list(age = age_numeric,
                  gender = gender,
                  weight = weights,
                  mu_k = mu_k_name,
                  sigma_k = sigma_k_name,
                  sum_prob = sum_prob_k_name)

## can just save these and use them again then

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

y_sim <- simulate_mixing(d = sim_pars_age$degree,
                         omega = sim_pars_age$omega,
                         M = sim_pars_age$M,
                         beta = sim_pars_age$beta,
                         ego = sim_pars_age$ego)


# Fit 2010 Mixing Model ---------------------------------------------------
num_egos <- 6
num_alters <- 8

mix_data_stan <- list(E = num_egos,
                      A = num_alters,
                      K = ncol(y_sim),
                      N = nrow(y_sim),
                      y = y_sim,
                      ego = sim_pars_age$ego,
                      Beta = sim_pars_age$beta,
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

## compare to true distribution
true_degrees <- tibble(node = 1:nrow(y_sim),
                       true_degree = sim_pars_age$degree)

true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")

## then do ppc for this model also
ppc_2010 <- construct_ppc(stan_fit_2010, y_sim)

ppc_fit_2010 <- plot_ests(ppc_2010$ppc_draws,
                          ppc_2010$y_tibble,
                          prop_val = 1)
ppc_fit_2010 + 
  labs(title = "McCormick et al, 2010") 


# Fit the 2019 Spline Model -----------------------------------------------

age_grid <- seq(min(data_pars$age), max(data_pars$age), 1);
knots <- quantile(age_grid, probs = seq(0, 1, .1));
N_K <- length(knots);
D <- 3;


## all the arguments to the stan fit are
stan_data_2019 <- list(K= ncol(y_sim),
                  N = nrow(y_sim),  
                  age_mean = age_mean_2014, 
                  Y = y_sim,  
                  w = data_pars$weight,
                  age = data_pars$age, 
                  g_n = as.numeric(data_pars$gender), 
                  g_k = c(rep(1, 6), rep(2, 6)), 
                  mu_k = mu_k_name,
                  sigma_k = sigma_k_name, 
                  sum_prob_k = sum_prob_k_name, 
                  mu_d = 6, sigma_d = 0.6, 
                  alpha_omega = 4.5,
                  beta_omega = 0.5, 
                  mu_lambda = log(100),
                  sigma_lambda = 0.5, 
                  alpha_rho = c(5, 5), 
                  mu_beta = c(6, 0.1, -3),
                  sigma_beta = rep(1, 3),
                  recall_power = 0, 
                  degree_regression = 1,
                  N_K = N_K,
                  N_S = length(age_grid),
                  D = D,
                  knots = knots,
                  X = age_grid)


stan_file_2019 <- here("Summer_2024", "2019_degree_kernel_spline.stan")

mod_2019 <- cmdstan_model(stan_file_2019)

stan_fit_2019 <- mod_2019$sample(data = stan_data_2019,
                                 seed = 123,
                                 chains = 4,
                                 iter_sampling = 1000,
                                 iter_warmup = 1000,
                                 parallel_chains = 4,
                                 refresh = 100)

stan_fit_2019$summary(variables = c("log_d[1]", "log_d[2]",
                                    "log_d[100]", "y_sim[1,2]"))

## then compare the estimated degrees and the ppc for this also

est_degrees_2019 <- stan_fit_2019$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("d[")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("d["),
               names_to = "node",
               values_to = "est_degree") |> 
  mutate(node = parse_number(node)) 

est_degrees_2019 |> 
  ggplot(aes(est_degree)) +
  geom_histogram() +
  labs(title = "Posterior Estimates",
       subtitle = "2019 Model")

## compare to true distribution
true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")


## then do the ppc down here

ppc_2019 <- construct_ppc(stan_fit_2019, y_sim)

ppc_fit_2019 <- plot_ests(ppc_2010$ppc_draws,
                          ppc_2010$y_tibble,
                          prop_val = 1)
ppc_fit_2019 + 
  labs(title = "Sahai et al, 2019") 

