library(rotasym)
library(tidyverse)
library(tidymodels) # to use coord_obs_pred
library(cmdstanr)
library(here)
library(bayesplot)
library(posterior)
library(grid)
library(gridExtra)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())

source(here("helper", "helper_model_checking.R"))
source(here("helper", "helper_latent_surface_model.R"))
source(here("helper", "helper_functions_latent_ard.R"))
source(here("helper", "helper_plots.R"))
set.seed(100)
# Simulate the Positions and Overall Network -----------------------------------

n_population <- 1e5
n_subpop <- 15   # number of subpops observed
n_sample <- 1000 # the people where ARD is recorded


###
true_subpops <- round(runif(n_subpop, min = n_sample/5, max = n_population / 20))
perc_subpop <- true_subpops/n_population
###

subpop_centers <- r_vMF(n = n_subpop, mu = c(1, 0, 0), kappa = 0.5)

xi <- 2.5
eta_vec <- rgamma(n = n_subpop, shape = 1, rate = 1)

mu <- c(0, 0, 1)
kappa <- 2
sample_pos <- r_vMF(n = n_sample, mu = mu, kappa = kappa)
alpha <- rnorm(n = n_sample, mean = 4, sd = 1)

dist <- sample_pos %*% t(subpop_centers)
theta <- acos(dist)

y <- matrix(NA, nrow = n_sample, ncol = n_subpop)
for (i in 1:n_sample) {
  for (k in 1:n_subpop) {
    num <- c_vMF(p = 3, kappa = xi) * c_vMF(p = 3, kappa = eta_vec[k])
    denom_term <- sqrt(xi ^ 2 + eta_vec[k] ^ 2 +
                         2 * xi * eta_vec[k] * cos(theta[i,k])) 
    denom <- c_vMF(p = 3, 0) * c_vMF(p = 3, kappa = denom_term)
    rate <- exp(alpha[i]) * perc_subpop[k] * num/denom
    y[i, k] <- rpois(n = 1, lambda = rate)
  }
}
y_sim <- y

known_pops <- 1:14
Pg1 <- sum(true_subpops[known_pops]/n_population)
G1_ind <- known_pops


sim_data <- list(y_sim = y_sim,
                 n_sample = n_sample,
                 n_population = n_population,
                 n_subpop = n_subpop,
                 true_alpha = alpha,
                 known_pops = known_pops,
                 G1_ind = G1_ind,
                 true_subpops = true_subpops)

saveRDS(sim_data, file = here("Summer_2025", "ard_latent_mod.RDS"))

stan_data_null <- list(N = nrow(y),
                       K = ncol(y),
                       y = y,
                       n_known = length(G1_ind),
                       idx = G1_ind,
                       known_prev = Pg1)

stan_file_null_01 <- here("stan_models", "null_model_01_scaled.stan")
mod_null_01 <- cmdstan_model(stan_file_null_01)

stan_fit_null_01 <- mod_null_01$sample(data = stan_data_null,
                                       seed = 123,
                                       chains = 4,
                                       iter_sampling = 1000,
                                       iter_warmup = 1000,
                                       parallel_chains = 4,
                                       refresh = 100)

stan_fit_null_01$summary(variables = c("scaled_beta"))
mcmc_trace(stan_fit_null_01$draws(), pars = "scaled_log_d")

samp_degree <- exp(alpha)

stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_log_d")) |> 
  mutate(draw = row_number())  |> 
  mutate(degree = exp(scaled_log_d)) |> 
  ggplot(aes(degree)) +
  geom_histogram()

hist(samp_degree)

subpop_info <- tibble(subpop = 1:15,
                      size = true_subpops)

stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("b[")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("b"),
               names_to = "par", 
               values_to = "sample") |> 
  mutate(subpop = parse_number(par),
         subpop_size = n_population * sample) |> 
  # filter(!(subpop %in% G1_ind)) |> 
  ggplot(aes(subpop_size)) +
  geom_histogram() +
  facet_wrap(~subpop, scales = "free", nrow = 3) +
  geom_vline(data = subpop_info, aes(xintercept = size), col = "red")

y_sim <- y
## compute ppc
ppc_null_1 <- construct_ppc(stan_fit_null_01, y_sim)

(ppc_fit_null_1 <- plot_ests(ppc_null_1$ppc_draws, 
                             ppc_null_1$y_tibble,
                             prop_val = 1))

rm(ppc_null_1)
ppc_1_null_1 <- ppc_fit_null_1 + 
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()
ppc_1_null_1



# Fit Second Null Model ---------------------------------------------------



stan_file_null_02 <- here("stan_models", "null_model_02_scaled.stan")
mod_null_02 <- cmdstan_model(stan_file_null_02)
stan_fit_null_02 <- mod_null_02$sample(data = stan_data_null,
                                       seed = 123,
                                       chains = 4,
                                       iter_sampling = 1000,
                                       iter_warmup = 1000,
                                       parallel_chains = 4,
                                       refresh = 100)

stan_fit_null_02$summary(variables = c("scaled_beta",
                                       "scaled_log_d[1]",
                                       "scaled_log_d[2]"))

## check the population sizes
b_draws <- stan_fit_null_02$draws() |> 
  as_draws_df() |> 
  select(starts_with("b["))

size_ests <- b_draws * n_population
head(rowSums(size_ests[, G1_ind]))
sum(true_subpops[1:14])

stan_fit_null_02$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_log_d")) |> 
  mutate(draw = row_number())  |>
  pivot_longer(cols = starts_with("scaled_log_d"),
               names_to = "node", 
               values_to = "log_d") |> 
  mutate(degree = exp(log_d),
         node_id = parse_number(node)) |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  labs(title = "Estimated Degree Distribution",
       subtitle = "Null Model, varying Degrees")


stan_fit_null_02$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("b[")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("b"),
               names_to = "par", 
               values_to = "sample") |> 
  mutate(subpop = parse_number(par),
         subpop_size = n_population * sample) |> 
  # filter(!(subpop %in% G1_ind)) |> 
  ggplot(aes(subpop_size)) +
  geom_histogram() +
  facet_wrap(~subpop, scales = "free", nrow = 3) +
  geom_vline(data = subpop_info, aes(xintercept = size), col = "red")

ppc_null_2 <- construct_ppc(stan_fit_null_02, y_sim)

(ppc_fit_null_2 <- plot_ests(ppc_null_2$ppc_draws, 
                             ppc_null_2$y_tibble,
                             prop_val = 1))
rm(ppc_null_2)
ppc_1_null_2 <- ppc_fit_null_2 +
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()

ppc_1_null_2

## show the estimated degree from null models and compare to true distribution
null_1_deg <- stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_log_d")) |> 
  mutate(draw = row_number())  |> 
  mutate(degree = exp(scaled_log_d)) |> 
  select(draw, degree) |> 
  mutate(model = "Null Model 1")

null_2_deg <- stan_fit_null_02$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_log_d")) |> 
  mutate(draw = row_number())  |> 
  pivot_longer(cols = starts_with("scaled_log_d"),
               names_to = "node", 
               values_to = "log_d") |> 
  mutate(degree = exp(log_d),
         node_id = parse_number(node)) |> 
  select(draw, degree) |> 
  mutate(model = "Null Model 2")

true_deg <- tibble(draw = 1,
                   degree = samp_degree) |> 
  mutate(model = "True Degree")


null_deg_plot <- bind_rows(null_1_deg,
                           null_2_deg,
                           true_deg) |> 
  ggplot(aes(x = degree, color = model)) +
  # geom_density(size = 1, key_glyph = "path") +
  stat_density(aes(y = ..scaled..), geom = "line",
               key_glyph = "path",
               position = "identity") +
  guides(color = guide_legend(override.aes = list(fill = NA)),
         override.aes = list(linetype = 1, size = 1)) +
  labs(color = "", x = "Degree", y = "") +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(0, 750)) +
  theme_single_legend() 

null_deg_plot


# Fit Zheng 2006 Model ----------------------------------------------------

stan_file_zheng <- here("stan_models", "zheng_et_al_2006_scaled.stan")
mod_zheng <- cmdstan_model(stan_file_zheng)
stan_fit_zheng <- mod_zheng$sample(data = stan_data_null,
                                   seed = 123,
                                   chains = 4,
                                   iter_sampling = 1000,
                                   iter_warmup = 1000,
                                   parallel_chains = 4,
                                   refresh = 100)

stan_fit_zheng$summary(variables = c("sigma_alpha", "scaled_beta"))

## check it recovers the groups correctly
b_draws <- stan_fit_zheng$draws() |> 
  as_draws_df() |> 
  select(starts_with("scaled_beta"))

size_ests <- exp(b_draws) * n_population
head(rowSums(size_ests[, G1_ind]))
sum(true_subpops[1:14])


est_degrees_2006 <- stan_fit_zheng$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_alpha")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("scaled_alpha"),
               names_to = "node",
               values_to = "log_estimate") |> 
  mutate(node = parse_number(node)) |> 
  mutate(est_degree = exp(log_estimate)) 

## look at posterior estimates of degree and subpop size
stan_fit_zheng$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_alpha")) |> 
  mutate(draw = row_number())  |> 
  pivot_longer(cols = starts_with("scaled_alpha"), 
               names_to = "par", values_to = "sample") |> 
  mutate(node = parse_number(par)) |> 
  mutate(degree = exp(sample)) |> 
  ggplot(aes(degree)) +
  geom_histogram()



stan_fit_zheng$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_beta[")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("scaled_beta"),
               names_to = "par", 
               values_to = "sample") |> 
  mutate(subpop = parse_number(par),
         subpop_size = n_population * exp(sample)) |> 
  # filter(!(subpop %in% G1_ind)) |> 
  ggplot(aes(subpop_size)) +
  geom_histogram() +
  facet_wrap(~subpop, scales = "free", nrow = 3) +
  geom_vline(data = subpop_info, aes(xintercept = size), col = "red")



ppc_zheng <- construct_ppc(stan_fit_zheng, y_sim)

(ppc_fit_zheng <- plot_ests(ppc_zheng$ppc_draws, 
                            ppc_zheng$y_tibble,
                            prop_val = 1))

ppc_plot_zheng <- ppc_fit_zheng + 
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()


# Fit McCormick 2010 ------------------------------------------------------

## to fit this model need to specify ego and alter groups 
## in both the nodes in ard sample and also in the ard subpopulations
## also need to specify beta, which corresponds to proportions
## of the alter groups in a subpop
## will assume all alter groups equally weighted across all populations
## USING TRUE subpop sizes

num_egos <- 6
num_alters <- 6

## as we don't have ego alter mixing structure seems we should just
## simulate this somehow here based on assigned subpopulations maybe


## assign egos randomly

sample_egos <- sample(1:num_egos, size = n_sample, replace = TRUE)

beta_sim <- matrix(NA,
                   nrow = num_alters,
                   ncol = n_subpop)
## use this beta, assuming each supop 
## equally spread across the ego groups, will give the right subpop proportions
for(i in 1:n_subpop){
  beta_sim[, i] <- rep(true_subpops[i]/(n_population), num_egos)
}


mix_data_stan <- list(E = num_egos,
                      A = num_alters,
                      K = ncol(y_sim),
                      N = nrow(y_sim),
                      y = y_sim,
                      ego = sample_egos,
                      Beta = beta_sim,
                      theta_d = c(6.2, .5),
                      theta_o = c(3, 2),
                      alpha = rep(1, num_alters)/num_alters,
                      p = 1)

stan_file_2010 <- here("stan_models", "mc_cormick_et_al_2010.stan")
mod_2010 <- cmdstan_model(stan_file = stan_file_2010)

stan_fit_2010 <- mod_2010$sample(data = mix_data_stan,
                                 seed = 123,
                                 chains = 4,
                                 iter_sampling = 1000,
                                 iter_warmup = 1000,
                                 parallel_chains = 4,
                                 refresh = 100)

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
true_deg |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")

## then do ppc for this model also
ppc_2010 <- construct_ppc(stan_fit_2010, y_sim)

ppc_fit_2010 <- plot_ests(ppc_2010$ppc_draws,
                          ppc_2010$y_tibble,
                          prop_val = 1)
rm(ppc_2010)
ppc_fit_2010 + 
  labs(title = "McCormick et al, 2010") 
