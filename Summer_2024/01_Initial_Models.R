#### July 5th 2024 ####

## Initial simulation studies where we simulate data from and
## fit the model of Zheng et al 2006



# Setup -----

library(tidyverse)
library(stringr)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(here)
library(scales)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())

# Simulate Data -----

N <- 200
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
                  # mu_beta = mu_beta,
                  # sigma_beta = sigma_beta)


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

# stan_fit$summary()

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

num_nodes <- 20

est_out_degrees |> 
  filter(node_id <= num_nodes) |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  geom_vline(data = true_out_degrees |> filter(node_id <= num_nodes), 
             mapping = aes(xintercept = degree), col = "red") +
  facet_wrap(~node_id, scales = "free", ncol = 5) +
  labs(title = "Full Model, which has poor convergence") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 2))

## these look reasonable here for a model that fits quite well
## even with beta's being estimated these still look good
## could probably modify this to a scenario where it doesn't work


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
  filter(sub_pop_id <= num_nodes)


matrix_df <- as.data.frame(as.table(y))
colnames(matrix_df) <- c("node_id", "sub_pop_id", "count")
matrix_df$node_id <- as.numeric(matrix_df$node_id)
matrix_df$sub_pop_id <- as.numeric(matrix_df$sub_pop_id)


true_y <- as_tibble(matrix_df)

est_y |> 
  ggplot(aes(count)) +
  geom_histogram() +
  geom_vline(data = true_y |> filter(node_id == 1 & sub_pop_id <= num_nodes),
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
  geom_vline(data = true_out_degrees |> summarise(avg = mean(degree)), 
             mapping = aes(xintercept = avg), col = "red") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


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
                       iter_sampling = 500,
                       iter_warmup = 500,
                       parallel_chains = 4,
                       refresh = 100)


mcmc_hist(stan_fit_simple$draws(), pars = c("mu_alpha", "sigma_alpha"))
mcmc_trace(stan_fit_simple$draws(), pars = c("mu_alpha", "sigma_alpha"))


est_out_degrees_simple <- stan_fit_simple$draws() |> 
  as_draws_df() |> 
  select(starts_with("out")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("out"), values_to = "degree") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+"))) 


est_out_degrees_simple |> 
  filter(node_id <= num_nodes) |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  geom_vline(data = true_out_degrees |> filter(node_id <= num_nodes), 
             mapping = aes(xintercept = degree), col = "red") +
  facet_wrap(~node_id, scales = "free", ncol = 5) +
  labs(title = "Model with Known Beta distribution") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 2))

est_y_simple <- stan_fit_simple$draws() |> 
  as_draws_df() |> 
  select(starts_with("y_sim")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("y_sim"), values_to = "count") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+")),
         sub_pop_id = str_extract(name, pattern = "\\d+]"),
         sub_pop_id = as.numeric(str_replace(sub_pop_id, "\\]", ""))) |> 
  filter(node_id == 1) |> 
  filter(sub_pop_id < num_nodes)


est_y_simple |> 
  ggplot(aes(count)) +
  geom_histogram() +
  geom_vline(data = true_y |> filter(node_id == 1 & sub_pop_id <= num_nodes),
             mapping = aes(xintercept = count), col = "red") +
  facet_wrap(~sub_pop_id, scales = "free", ncol = 5) +
  labs(title = "Posterior Predictive of Entries of y, Simple Model")


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
  geom_vline(data = true_out_degrees |> summarise(avg = mean(degree)), 
             mapping = aes(xintercept = avg), col = "red") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


# Things to think about ---------------------------------------------------

## Other comparisons, such as plotting estimated against true alpha, beta

## Other possible ppc.