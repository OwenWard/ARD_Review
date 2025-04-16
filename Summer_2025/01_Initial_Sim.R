##### April 15th 2025 ####


## Modifying and extending the simulation setting to consider 
## more realistic data, along with

## 1. Correctly rescale the estimates so that they make sense
## 2. Consider PPC for the fitted models
## 3. Consider LOO-CV and similar metrics for model comparison




# Setup and Load Packages -------------------------------------------------

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

perc_subpop <- 0.01  # prevalence of subpop in population, same for all subpops
perc_subpop <- round(rgamma(n = n_subpop, shape = 1, rate = 10), 3)
num_subpop <- round(n_population * perc_subpop) # approx size of subpopos in pop

mu <- c(0, 0, 1)
kappa <- 2
latent_positions <- r_vMF(n = n_population, mu = mu, kappa = kappa)

all_greg <- rep(NA, n_population)

gen_greg <- -12 ## the population mean greg, if not in a subpop
sd_greg <- 1

## sample the gregariousness params of ARD sample
## these don't need to be integer necessarily
samp_greg <- rnorm(n_sample, mean = gen_greg, sd = sd_greg)

## the latent positions of the ARD sample
sample_ids <- sample(1:n_population, n_sample)

samp_pos <- latent_positions[sample_ids, ]

all_greg[sample_ids] <- samp_greg


# Simulate the Subpopulation Members and Structure------------------------------

## first simulate the subpopulation centers uniformly on the sphere

subpop_centers <- r_vMF(n = n_subpop, mu = c(1, 0, 0), kappa = 0.5)
sub_pop_id <- rep(NA, n_population)

## assume that different subpops have common dist for greg
## differences in mean gregariousness and smaller sd than
## the population
subpop_greg_means <- seq(from = 5, to = 7.5, length.out = n_subpop)
subpop_greg_sd <- 0.25

true_subpop_size <- rep(0, n_subpop)

lambda <- 10 ## scaling parameter in the distance
for(k in 1:n_subpop){
  
  ## randomly assign nodes to subpopulation based on distance
  ## to subpop center
  curr_center <- subpop_centers[k, ]
  
  dot_prod <- latent_positions %*% curr_center
  dist <- acos(dot_prod)
  ## these distances between 0 and 2pi
  prob_sub <- exp(- lambda * dist)
  prob_sub <- prob_sub /sum(prob_sub) * n_population * perc_subpop[k]
  prob_sub <- ifelse(prob_sub > 1, 1, prob_sub)
  sub_member <- rbinom(n_population, 1, p = prob_sub)
  cat(paste0(sum(sub_member), "\n"))  
  true_subpop_size[k] <- sum(sub_member)
  ## this number should be the approximate size of the subpopulation
  ## in the population, because want to randomly assign all nodes based on
  ## distance to center
  curr_subpop <- latent_positions[sub_member == 1, ]
  
  subpop_greg <- rnorm(n = sum(sub_member),
                       mean = subpop_greg_means[k],
                       sd = subpop_greg_sd)
  all_greg[sub_member == 1] <- subpop_greg
  sub_pop_id[sub_member == 1] <- k
  ## then need to store these subpopgreg parameters for those members of
  ## the pop
  
}


remaining_greg <- rnorm(n = sum(is.na(all_greg)), mean = gen_greg, sd = sd_greg)
all_greg[is.na(all_greg)] <- remaining_greg



# Simulate the True Degrees and Corresponding ARD-------------------------------

y_sim <- matrix(NA, nrow = n_sample, ncol = n_subpop)
samp_degree <- rep(NA, n_sample)
gamma <- -1 ## smaller values lead to larger degree

for(i in 1:n_sample){
  pop_id <- sample_ids[i]
  curr_greg <- samp_greg[i]
  curr_loc <- samp_pos[i, ]
  node_dist <- latent_positions %*% curr_loc
  exp_term <- curr_greg + all_greg + gamma * as.vector(node_dist)
  prob_edge <- 1/(1 + exp(-exp_term)) 
  true_edge <- rbinom(n_population, size = 1, prob = as.vector(prob_edge))
  samp_degree[i] <- sum(true_edge)
  ## then get the ARD data out of this
  ard_info <- tibble(id = 1:n_population, 
                     sub_pop_member = sub_pop_id,
                     edge = true_edge) |>
    drop_na() |> 
    filter(id != pop_id) |>  # to avoid self loops
    group_by(sub_pop_member) |> 
    summarise(n = sum(edge)) 
  y_sim[i, ard_info$sub_pop_member] <- ard_info$n
}

y_sim

## the true degree of the ard sample
hist(samp_degree, breaks = 50)

## how many large values are there here, in the ard sample
length(samp_degree[samp_degree > 1500])
summary(samp_degree)

## check the sparsity of the ard matrix
sum(y_sim > 0)/(n_sample * n_subpop)
summary(as.vector(y_sim))
hist(y_sim, breaks = 50)


# Plot the latent positions of Subpopulations -----------------------------
colnames(subpop_centers) <- c("x", "y", "z")
subpop_info <- as_tibble(subpop_centers) |> 
  mutate(id = row_number()) |> 
  mutate(theta = atan2(y, x),
         theta = ifelse(theta < 0, theta + 2 * pi, theta),
         phi = acos(z)) 

colnames(latent_positions) <- c("x", "y", "z")

plot_subpops <- as_tibble(latent_positions) |> 
  mutate(subpop = sub_pop_id) |> 
  # mutate(node = row_number()) |>  ## these lines to only show 
  # filter(node %in% sample_ids) |> ## nodes in ard sample
  drop_na() |> 
  mutate(theta = atan2(y, x),
         theta = ifelse(theta < 0, theta + 2 * pi, theta),
         phi = acos(z)) |>
  # slice_sample(n = 1000) |> 
  ggplot(aes(theta, phi, colour = as.factor(subpop))) +
  geom_point(alpha = 0.1) +
  geom_point(data = subpop_info, 
             mapping = aes(theta, phi, colour = as.factor(id)),
             alpha = 1.5, size = 5, pch = 4, stroke = 2) +
  theme(legend.position = "none") +
  labs(x = expression(theta), y = expression(phi)) +
  scale_x_continuous(breaks = c(pi, 2*pi), 
                     limits = c(0, 2*pi),
                     expand = c(0, 0.1),
                     labels = c(expression(pi), expression(2 * pi))) +
  scale_y_continuous(breaks = c(pi/2, pi),
                     expand = c(0, 0.1),
                     labels = c(expression(pi/2), expression(pi))) +
  theme_single_y()

plot_subpops



# Specify the Known Populations -------------------------------------------

known_pops <- c(1:5)


# Fit Simple Null Models --------------------------------------------------

stan_data_null <- list(N = nrow(y_sim),
                       K = ncol(y_sim),
                       y = y_sim)

stan_file_null_01 <- here("stan_models", "null_model_01.stan")

mod_null_01 <- cmdstan_model(stan_file_null_01)

stan_fit_null_01 <- mod_null_01$sample(data = stan_data_null,
                                       seed = 123,
                                       chains = 4,
                                       iter_sampling = 1000,
                                       iter_warmup = 1000,
                                       parallel_chains = 4,
                                       refresh = 100)

stan_fit_null_01$summary(variables = c("log_d", "beta"))

#### can we rescale these using the known subpopulations here?




## look at the estimated degree

stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("log_d")) |> 
  mutate(draw = row_number())  |> 
  mutate(degree = exp(log_d)) |> 
  ggplot(aes(degree)) +
  geom_histogram()

## look at simple ppc
ppc_null_1 <- construct_ppc(stan_fit_null_01, y_sim)

ppc_fit_null_1 <- plot_ests(ppc_null_1$ppc_draws, 
                            ppc_null_1$y_tibble,
                            prop_val = 1)

ppc_1_null_1 <- ppc_fit_null_1 + 
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()

ggsave(filename = here("Summer_2024", "figures",
                       "latent_ppc_null_1.png"),
       plot = ppc_1_null_1,
       dpi = 600,
       height = 5, width = 7)


## Fit the second null model, allowing varying degrees

stan_file_null_02 <- here("helper", "null_model_02.stan")
mod_null_02 <- cmdstan_model(stan_file_null_02)
stan_fit_null_02 <- mod_null_02$sample(data = stan_data_null,
                                       seed = 123,
                                       chains = 4,
                                       iter_sampling = 1000,
                                       iter_warmup = 1000,
                                       parallel_chains = 4,
                                       refresh = 100)

stan_fit_null_02$summary(variables = c("beta", "log_d[1]", "log_d[2]"))

stan_fit_null_02$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("log_d")) |> 
  mutate(draw = row_number())  |>
  pivot_longer(cols = starts_with("log_d"),
               names_to = "node", 
               values_to = "log_d") |> 
  mutate(degree = exp(log_d),
         node_id = parse_number(node)) |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  labs(title = "Estimated Degree Distribution",
       subtitle = "Null Model, varying Degrees")

ppc_null_2 <- construct_ppc(stan_fit_null_02, y_sim)

ppc_fit_null_2 <- plot_ests(ppc_null_2$ppc_draws, 
                            ppc_null_2$y_tibble,
                            prop_val = 1)

ppc_1_null_2 <- ppc_fit_null_2 +
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()

ggsave(filename = here("Summer_2024", "figures",
                       "latent_ppc_null_2.png"),
       plot = ppc_1_null_2,
       dpi = 600,
       height = 5, width = 7)


## show the estimated degree from null models and compare to true distribution
null_1_deg <- stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("log_d")) |> 
  mutate(draw = row_number())  |> 
  mutate(degree = exp(log_d)) |> 
  select(draw, degree) |> 
  mutate(model = "Null Model 1")

null_2_deg <- stan_fit_null_02$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("log_d")) |> 
  mutate(draw = row_number())  |> 
  pivot_longer(cols = starts_with("log_d"),
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

ggsave(filename = here("Summer_2024", "figures",
                       "latent_null_deg.png"),
       plot = null_deg_plot,
       dpi = 600,
       height = 5, width = 7)

## compute loo-cv for this between these models

library(loo)
loo_1 <- stan_fit_null_01$loo(cores = 4)
loo_2 <- stan_fit_null_02$loo(cores = 4)

loo_compare(loo_1, loo_2)
