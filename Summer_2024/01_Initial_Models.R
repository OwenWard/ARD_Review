#### July 5th 2024 ####

## Initial simulation studies where we simulate data from and
## fit the model of Zheng et al 2006 and McCormick et al 2010



# Setup -----

library(tidyverse)
library(stringr)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(here)
library(scales)
library(gtools)
library(rstan)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())


theme_no_y <- function(){ 
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

# Simulate Data -----

N <- 500
K <- 32

mu_alpha <- 5 ## possible choices are -2, 5, 7.5,
mu_beta <- -5
sigma_alpha <- 1 ## possible choices are 1, 5
sigma_beta <- 1
alpha <- rnorm(N, mu_alpha, sigma_alpha) ## log gregariousness of individuals
beta <- rnorm(K, mu_beta, sigma_beta) ## log prevalence of subgroups

omega_inv <- runif(K, 0.1, 0.95)
omega <- 1 / omega_inv
y <- array(dim = c(N, K))
for (i in 1:N) {
  for (k in 1:K) {
    xi_i_k <- exp(alpha[i] + beta[k]) / (omega[k] - 1)
    y[i, k] <- rnbinom(1,
                       size = xi_i_k,
                       prob = 1 / omega[k]
    )
  }
}



# Fit Stan Model -----

mu_beta <- -5 ## specify some parameters to make identifiable
sigma_beta <- 1
stan_data <- list(N = N,
                  K = K,
                  y = y)

stan_file <- here("Summer_2024", "zheng_et_al_06.stan")

mod <- cmdstan_model(stan_file)

stan_fit <- mod$sample(data = stan_data,
                       seed = 123,
                       chains = 4,
                       iter_sampling = 500,
                       iter_warmup = 500,
                       parallel_chains = 4,
                       refresh = 100)


# Initial Model Checking of this Model -----

stan_fit$summary(variables = c("mu_alpha", "sigma_alpha"))

mcmc_hist(stan_fit$draws(), pars = c("mu_alpha", "sigma_alpha"))

mcmc_trace(stan_fit$draws(), pars = c("mu_alpha", "sigma_alpha"))

## model is ok but clear identifiability issue



# Initial Posterior Predictive Checks -------------------------------------


## look at posterior predictive checking of the out degrees of the
## nodes in the network

est_out_degrees <- stan_fit$draws() |> 
  as_draws_df() |> 
  select(starts_with("out")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("out"), values_to = "degree") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+"))) 


true_out_degrees <- tibble(node_id = 1:N, degree = apply(y, 1, sum)) 

num_nodes <- c(1:20)

est_out_degrees |> 
  filter(node_id %in% num_nodes) |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  geom_vline(data = true_out_degrees |> filter(node_id %in% num_nodes), 
             mapping = aes(xintercept = degree), col = "red") +
  facet_wrap(~node_id, scales = "free", ncol = 5) +
  labs(title = "Full Model, which has poor convergence",
       x = "Sample Out Degree", y = element_blank()) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 2)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

## these look reasonable here for a model that fits quite well
## even with beta's being estimated these still look good
## could probably modify this to a scenario where it doesn't work

## see where the biggest differences are

worst <- est_out_degrees |> 
  left_join(true_out_degrees |> rename(true_degree = degree), 
            by = "node_id") |> 
  mutate(diff = true_degree - degree) |> 
  group_by(node_id) |> 
  summarise(mean_diff = mean(diff^2)) |>
  arrange(-abs(mean_diff)) |>
  slice_max(order_by = mean_diff, n = 10) |> 
  pull(node_id)

est_out_degrees |> 
  filter(node_id %in% worst) |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  geom_vline(data = true_out_degrees |> filter(node_id %in% worst), 
             mapping = aes(xintercept = degree), col = "red") +
  facet_wrap(~node_id, scales = "free", ncol = 5) +
  labs(title = "Full Model, which has poor convergence",
       subtitle = "Out degrees which are poorly predicted") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 2)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


## even for these poor examples the out degrees look quite good
## is this due to sparsity or something?

## look at the first row of y instead

est_y <- stan_fit$draws() |> 
  as_draws_df() |> 
  select(starts_with("y_sim")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("y_sim"), values_to = "count") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+")),
         sub_pop_id = str_extract(name, pattern = "\\d+]"),
         sub_pop_id = as.numeric(str_replace(sub_pop_id, "\\]", ""))) |> 
  filter(node_id == 1) |> 
  filter(sub_pop_id %in% num_nodes)


matrix_df <- as.data.frame(as.table(y))
colnames(matrix_df) <- c("node_id", "sub_pop_id", "count")
matrix_df$node_id <- as.numeric(matrix_df$node_id)
matrix_df$sub_pop_id <- as.numeric(matrix_df$sub_pop_id)
true_y <- as_tibble(matrix_df)

est_y |> 
  ggplot(aes(count)) +
  geom_histogram() +
  geom_vline(data = true_y |> filter(node_id == 1 & sub_pop_id %in% num_nodes),
             mapping = aes(xintercept = count), col = "red") +
  facet_wrap(~sub_pop_id, scales = "free", ncol = 5) +
  labs(title = "Posterior Predictive of Entries of y, Full Model") +
  scale_x_continuous(breaks = breaks_pretty(n = 2))


## distribution of average gregariousness also
est_greg <- stan_fit$draws() |> 
  as_draws_df() |> 
  select(starts_with("alpha")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("alpha"), values_to = "est_alpha") |> 
  mutate(est_greg = exp(est_alpha))

est_greg |> 
  group_by(draw) |> 
  summarise(avg = mean(est_greg)) |> 
  ggplot(aes(avg)) +
  geom_histogram() +
  labs(y = element_blank(),
       x = "Mean Gregariousness Parameter", 
       title = "Full Model") +
  # geom_vline(data = true_out_degrees |> summarise(avg = mean(degree)), 
  #            mapping = aes(xintercept = avg), col = "red") +
  geom_vline(aes(xintercept = mean(exp(alpha))), col = "red") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


## not sure the average is accurate, given the potential long tail
## will instead just plot the distribution of gregariousness, compare
## it to the truth

true_greg <- tibble(greg = exp(alpha))

est_greg |> 
  ggplot(aes(est_greg)) +
  geom_histogram(aes(y = after_stat(density))) +
  geom_histogram(data = true_greg, 
                 mapping = aes(x = greg,
                               y = after_stat(density)),
                 fill = "red", alpha = 0.75) +
  # scale_x_log10() +
  labs(x = "Log Gregariousness") +
  NULL

## these not comparable at all here now


# Repeat for the simpler Bayesian Model -----------------------------------


stan_data_simple <- list(N = N,
                         K = K,
                         y = y,
                         mu_beta = mu_beta,
                         sigma_beta = sigma_beta)


stan_file_simple <- here("Summer_2024", "zheng_et_al_06_simple.stan")

mod_simple <- cmdstan_model(stan_file_simple)

stan_fit_simple <- mod_simple$sample(data = stan_data_simple,
                       seed = 123,
                       chains = 4,
                       # iter_sampling = 500,
                       # iter_warmup = 500,
                       parallel_chains = 4,
                       refresh = 100)

stan_fit_simple$summary(variables = c("mu_alpha", "sigma_alpha"))
mcmc_hist(stan_fit_simple$draws(), pars = c("mu_alpha", "sigma_alpha"))
mcmc_trace(stan_fit_simple$draws(), pars = c("mu_alpha", "sigma_alpha"))

est_out_degrees_simple <- stan_fit_simple$draws() |> 
  as_draws_df() |> 
  select(starts_with("out")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("out"), values_to = "degree") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+"))) 

## try identify the most wrong ones
worst_simple <- est_out_degrees_simple |> 
  left_join(true_out_degrees |> rename(true_degree = degree), 
            by = "node_id") |> 
  mutate(diff = true_degree - degree) |> 
  group_by(node_id) |> 
  summarise(mean_diff = mean(diff^2)) |>
  arrange(-abs(mean_diff)) |> 
  slice_max(order_by = mean_diff, n = 10) |> 
  pull(node_id)

est_out_degrees_simple |> 
  filter(node_id %in% worst_simple) |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  geom_vline(data = true_out_degrees |> filter(node_id %in% worst_simple), 
             mapping = aes(xintercept = degree), col = "red") +
  facet_wrap(~node_id, scales = "free", ncol = 5) +
  labs(title = "Model with Known Beta distribution") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 2)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

est_y_simple <- stan_fit_simple$draws() |> 
  as_draws_df() |> 
  select(starts_with("y_sim")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("y_sim"), values_to = "count") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+")),
         sub_pop_id = str_extract(name, pattern = "\\d+]"),
         sub_pop_id = as.numeric(str_replace(sub_pop_id, "\\]", ""))) |> 
  filter(node_id == 1) |> 
  filter(sub_pop_id %in% num_nodes)

est_y_simple |> 
  ggplot(aes(count)) +
  geom_histogram() +
  geom_vline(data = true_y |> filter(node_id == 1 & sub_pop_id %in% num_nodes),
             mapping = aes(xintercept = count), col = "red") +
  facet_wrap(~sub_pop_id, scales = "free", ncol = 5) +
  labs(title = "Posterior Predictive of Entries of y, Simple Model") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


## plot estimated degrees after transforming draws of alpha

est_greg_simple <- stan_fit_simple$draws() |> 
  as_draws_df() |> 
  select(starts_with("alpha")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("alpha"), values_to = "est_alpha") |> 
  mutate(est_greg = exp(est_alpha))

est_greg_simple |> 
  group_by(draw) |> 
  summarise(avg = mean(est_greg)) |> 
  ggplot(aes(avg)) +
  geom_histogram() +
  labs(y = element_blank(),
       x = "Mean Gregariousness Parameter", 
       title = "Simple Model") +
  geom_vline(aes(xintercept = mean(exp(alpha))), col = "red") +
  # geom_vline(data = true_out_degrees |> summarise(avg = mean(degree)), 
  #            mapping = aes(xintercept = avg), col = "red") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


est_greg_simple |> 
  ggplot(aes(est_greg)) +
  geom_histogram(aes(y = after_stat(density))) +
  geom_histogram(data = true_greg, 
                 mapping = aes(x = greg,
                               y = after_stat(density)),
                 fill = "red", alpha = 0.5) +
  scale_x_log10() +
  labs(x = "Log Gregariousness") +
  NULL

## for this model the distribution looks pretty good, seems to make sense

# Things to think about ---------------------------------------------------

## Other comparisons, such as plotting estimated against true alpha, beta

## Other possible ppc.



# Trying to Recreate Figure 2 in Sahai et al ------------------------------

## to do this need to specify 14 names it seems

load(here("Sahai_et_al_2018/", "Swupnil_Code", "data",
          "mix_sim", "mix_mat.Rdata"))

## this only loads in mixing.sim, no other data/variables, 
## which has 6 rows and 8 columns


n_names <- 14

omega_sim <- runif(n_names, 0.85, 1)

M_sim <- matrix(c(rdirichlet(1, mixing.sim[1, ]),
                  rdirichlet(1, mixing.sim[2, ]),
                  rdirichlet(1, mixing.sim[3, ]),
                  rdirichlet(1, mixing.sim[4, ]),
                  rdirichlet(1, mixing.sim[5, ]),
                  rdirichlet(1, mixing.sim[6, ])),
                     nrow = 6, ncol = 8, byrow = TRUE)

## after running source(paste0(code_path,"src/load_data/occs.R"))
## there is then an object called beta_names,
## along with ego_sex_age

code_path <- here("Sahai_et_al_2018/Swupnil_Code/")
source(paste0(code_path,"src/load_data/surveys.R"))
source(paste0(code_path,"src/load_data/names.R"))
source(paste0(code_path,"src/load_data/occs.R"))

## think there are 6 ego groups and 8 alter groups

beta_sim <- matrix(rexp(14 * 8, 1/mean(beta_names[beta_names>0])),
                   nrow = 8, ncol = 14)
colnames(beta_sim) <- paste("Name", 1:ncol(beta_sim))
rownames(beta_sim) <- paste("Alter", 1:nrow(beta_sim))

ego_sim <- ego_sex_age


## need to find beta_names and degree which has
## CombSexAge in it
## I think this should be some of the real data, but not sure which part
## think it's from degree_mix.csv

degree <- read.csv(here("Sahai_et_al_2018/", "Swupnil_Code",
                        "estimates", "degree_mix.csv"))

degree_sim <- degree$CombSexAge

source(here("Sahai_et_al_2018", "Swupnil_code", "src",
            "funcs", "mix_matrix.R"))

names_data_sim <- simulate_mixing(degree_sim, omega_sim,
                                  M_sim, beta_sim, ego_sim)

names_data_sim


# for(k in c(4,6,8,10,12,14)) {
k <- 4 ## the number of names used in the fitting
mcmc_data_sim_k <- list(E = 6, A = 8, K = k,
                        N = nrow(names_data_sim),
                        y = names_data_sim[, c(1:k)],
                        ego = ego_sim,
                        Beta = beta_sim[, c(1:k)],
                        theta_d = c(6.2, .5),
                        theta_o = c(3, 2),
                        alpha = rep(1, 8)/8,
                        p = 1)

## not sure about the value of p here, can be in [0.5, 1]
## then load in the stan code for this model

## should update this to cmdstanr

source(here("Sahai_et_al_2018/", "Swupnil_code", "archive", 
            "Degree_Mixing_Code.Stan.R"))

degree_mixing_fit <- stan_model(model_code = degree_mixing_code,
                                model_name = "Degree")

fit_comb_k <- sampling(degree_mixing_fit, 
                       data = mcmc_data_sim_k,
                       iter = 1000, chains = 4)

plot_mix_comp(fit_comb_k, M_sim, paste(k, "Names"))
# }



fit_comb_k


## updating this stan code 

degree_mix_mod <- cmdstan_model(stan_file = here("Summer_2024",
                                                 "Degree_Mixing.stan"))

fit <- degree_mix_mod$sample(data = mcmc_data_sim_k)


## then can look at histogram of the degrees
## here I think the degree are their true degrees, because 
## mixing accounts for proportion across population

fit$summary()


degree_draws <- fit$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("log_d")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("log_d"), values_to = "log_degree") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+"))) 


true_deg <- tibble(degree = degree_sim) |> 
  mutate(log_degree = log(degree), node_id = row_number())

num_nodes <- 20

degree_draws |> 
  filter(node_id <= num_nodes) |> 
  ggplot(aes(log_degree)) +
  geom_histogram() + 
  geom_vline(data = true_deg |> filter(node_id <= num_nodes), 
             mapping = aes(xintercept = log_degree),
             col = "red") +
  facet_wrap(~node_id, scales = "free") +
  theme_no_y()


## could check the proportion

degree_draws |> 
  group_by(node_id) |> 
  summarise(lower = quantile(log_degree, prob = 0.025),
            upper = quantile(log_degree, prob = 0.975)) |> 
  left_join(true_deg, by = "node_id") |> 
  rename(true_log_deg = log_degree) |> 
  rowwise() |> 
  mutate(in_int = ifelse(true_log_deg > lower & true_log_deg < upper,
                         1, 0)) |> 
  ungroup() |> 
  summarise(sum(in_int)/n())

## so basically 95% coverage here
