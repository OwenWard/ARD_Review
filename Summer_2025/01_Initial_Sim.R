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

# Fit Simple Null Models, with scaling----------------------------------------


# Model 1 --------------------------------------
G1_ind <- known_pops
stan_data_null <- list(N = nrow(y_sim),
                       K = ncol(y_sim),
                       y = y_sim,
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

stan_fit_null_01$summary()
mcmc_trace(stan_fit_null_01$draws(), pars = "scaled_log_d")


## check the population sizes here for this
b_draws <- stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  select(starts_with("b["))

size_ests <- b_draws * n_population
head(rowSums(size_ests[, G1_ind]))
sum(true_subpop_size[1:5])

## look at posterior estimates of degree and subpop size
stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_log_d")) |> 
  mutate(draw = row_number())  |> 
  mutate(degree = exp(scaled_log_d)) |> 
  ggplot(aes(degree)) +
  geom_histogram()

hist(samp_degree)

## show the true subpopulation sizes
subpop_info <- tibble(subpop = 1:n_subpop,
       size = true_subpop_size)

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


## compute ppc
ppc_null_1 <- construct_ppc(stan_fit_null_01, y_sim)

(ppc_fit_null_1 <- plot_ests(ppc_null_1$ppc_draws, 
                            ppc_null_1$y_tibble,
                            prop_val = 1))

ppc_1_null_1 <- ppc_fit_null_1 + 
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()
ppc_1_null_1



# Fit Model 2 with Scaling  -----------------------------------------------


## Fit the second null model, allowing varying degrees

stan_file_null_02 <- here("stan_models", "null_model_02_scaled.stan")
mod_null_02 <- cmdstan_model(stan_file_null_02)
stan_fit_null_02 <- mod_null_02$sample(data = stan_data_null,
                                       seed = 123,
                                       chains = 4,
                                       iter_sampling = 4000,
                                       iter_warmup = 1000,
                                       parallel_chains = 4,
                                       refresh = 100)

stan_fit_null_02$summary(variables = c("scaled_beta", "scaled_log_d[1]",
                                       "scaled_log_d[2]"))

## check the population sizes
b_draws <- stan_fit_null_02$draws() |> 
  as_draws_df() |> 
  select(starts_with("b["))

size_ests <- b_draws * n_population
head(rowSums(size_ests[, G1_ind]))
sum(true_subpop_size[1:5])

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


## compute loo-cv for this between these models

library(loo)
loo_1 <- stan_fit_null_01$loo(cores = 4)
loo_2 <- stan_fit_null_02$loo(cores = 4)

loo_compare(loo_1, loo_2)


## maybe the fact that loo-cv not reasonable for these models is enough...
## need to dig into this a bit more




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
sum(true_subpop_size[1:5])

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

loo_est <- stan_fit_zheng$loo(cores = 4)
## these look good, which means we can use it for this data!!! (hopefully)




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

## create random alter centers and construct beta based on distance between
## these and the subpop centers (such that they are reasonably small)
## then that's all i need
alter_centers <- r_vMF(n = num_alters, mu = c(1, 1, 1), kappa = 0.05)

beta_sim <- matrix(NA,
                   nrow = num_alters,
                   ncol = n_subpop)

lambda_alter <- 5
## make this large enough so that prob sub smallish

for(i in 1:num_alters){
  ###
  curr_alter <- alter_centers[i, ]
  dot_prod <- subpop_centers %*% curr_alter
  dist <- acos(dot_prod)
  ## these distances between 0 and 2pi
  prob_sub <- exp(- lambda_alter * dist)
  beta_sim[i, ] <- prob_sub
}

## could potentially use this beta, assuming each supop 
## equally spread across the ego groups, will give the right subpop proportions
# for(i in 1:n_subpop){
#   beta_sim[, i] <- rep(true_subpop_size[i]/(n_population), num_egos)
# }


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
ppc_fit_2010 + 
  labs(title = "McCormick et al, 2010") 


## look at the subpopulation estimates
stan_fit_2010$summary(variables = c("M")) |> print(n = 36)

# dot_product(M[ego[n]], sub_col(Beta, 1, k, A))
## this is it, from stan
mix_data_stan$Beta
    

loo_est_2010 <- stan_fit_2010$loo(cores = 4)
loo_est_2010


loo_compare(loo_est, loo_est_2010)



# Fit McCormick and Zheng 2015 --------------------------------------------


dim(y_sim)

ls.dim <- 3
n <- dim(y_sim)[1]
n.iter <- 5000#3000
n.thin <- 10
m.iter <- 4
total.prop <- 0.25

## taking this from the github
muk.fix.ind <- sample(1:8, size = 4, replace = F)
muk.fix <- matrix(runif(12), nrow = 4, ncol = 3)
muk.fix <- sweep(muk.fix, MARGIN = 1, 1 / sqrt(rowSums(muk.fix^2)), `*`)

z.pos.init <- generateRandomInitial(n, ls.dim)
out <- f.metro(y_sim,
               total.prop = total.prop,
               n.iter = n.iter,
               m.iter = m.iter,
               n.thin = n.thin,
               z.pos.init = z.pos.init,
               muk.fix = muk.fix,
               ls.dim = ls.dim)

## Use these estimated model fits

posterior <- getPosterior(out, n.iter, m.iter, n.thin, n)
est.degrees <- posterior$est.degrees

colnames(est.degrees) <- paste0("node", 1:n)
est_degrees_2015 <- as_tibble(est.degrees) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("node"), names_to = "node",
               values_to = "est_log_degree") |> 
  mutate(degree = exp(est_log_degree))


est_degrees_2015 |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  labs(title = "Degree Estimates",
       subtitle = "2015 Model")

true_deg |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")


all_post <- get_all_post(out, n.iter, m.iter, n.thin, n, k)
num_sims <- 1000
## this can likely be parallelized, but works for now

sim_y <- function(all_post, n, k, num_sims, ls.dim = 3){
  y_sim <- array(NA, dim = c(n, k, num_sims))
  for(sim in 1:num_sims){
    ## get the current estimates of everything
    curr_latent_pos <- matrix(all_post$est.latent.pos[sim, ],
                              byrow = F,
                              nrow = n,
                              ncol = ls.dim)
    curr_mu <- matrix(all_post$est.mu.k[sim, ],
                      byrow = F, nrow = k, ncol = ls.dim)
    curr_eta <- all_post$est.eta[sim] 
    curr_degrees <- all_post$est.degrees[sim, ]
    curr_beta <- all_post$est.beta[sim, ]
    curr_eta_k <- all_post$eta.eta.k[sim, ]
    ## simulate the entries of y and populate y_sim
    for(i in 1:n){
      curr_deg <- curr_degrees[i]
      curr_pos <- curr_latent_pos[i, ]
      curr_prod <- acos(curr_pos %*% t(curr_mu) )
      sqrt_term <- sqrt(curr_eta ^ 2 + curr_eta_k^2 + 
                          2 * curr_eta * curr_eta_k * cos(curr_prod))
      cp_term <- cp_fcn(curr_eta) * cp_fcn(curr_eta_k) / (cp_fcn(0) * 
                                                            cp_fcn(sqrt_term))
      rate_vec <- exp(curr_deg) * exp(curr_beta) * cp_term
      y_sim[i, ,sim] <- mapply(rpois, n = 1, lambda = rate_vec)
    }
  }
  y_sim
}

y_sim_ppc <- sim_y(all_post, n, k, num_sims, ls.dim = 3)

dim(y_sim_ppc)

index_grid <- expand_grid(
  Index1 = seq_len(dim(y_sim_ppc)[1]),
  Index2 = seq_len(dim(y_sim_ppc)[2]),
  Index3 = seq_len(dim(y_sim_ppc)[3])
)

# Add the values from the array to the grid
tibble_df <- index_grid |> 
  mutate(Value = y_sim_ppc[cbind(Index1, Index2, Index3)])

ard_ppc_draws <- tibble_df |> 
  rename(node_id = Index1,
         sub_pop_id = Index2,
         draw = Index3,
         count = Value)

## then do the ppc with this data
ppc_fit_2015 <- plot_ests(ard_ppc_draws,
                          ppc_2010$y_tibble,
                          prop_val = 1)

ppc_fit_2015 + 
  labs(title = "McCormick et al, 2015")


## compute log-likelihood instead
## this a quick solution, can make it precise after
## needs to actually be a matrix maybe here, need to figure it 
## out

## need to fix this so it can be passed into loo
## does stan flatten the array of llh row wise or column wise?

llh_y <- function(all_post, n, k, num_sims, ls.dim = 3, y_sim){
  llh <- array(NA, dim = c(num_sims, n, k))
  for(sim in 1:num_sims){
    ## get the current estimates of everything
    curr_latent_pos <- matrix(all_post$est.latent.pos[sim, ],
                              byrow = F,
                              nrow = n,
                              ncol = ls.dim)
    curr_mu <- matrix(all_post$est.mu.k[sim, ],
                      byrow = F, nrow = k, ncol = ls.dim)
    curr_eta <- all_post$est.eta[sim] 
    curr_degrees <- all_post$est.degrees[sim, ]
    curr_beta <- all_post$est.beta[sim, ]
    curr_eta_k <- all_post$eta.eta.k[sim, ]
    ## simulate the entries of y and populate y_sim
    for(i in 1:n){
      curr_deg <- curr_degrees[i]
      curr_pos <- curr_latent_pos[i, ]
      curr_prod <- acos(curr_pos %*% t(curr_mu) )
      sqrt_term <- sqrt(curr_eta ^ 2 + curr_eta_k^2 + 
                          2 * curr_eta * curr_eta_k * cos(curr_prod))
      cp_term <- cp_fcn(curr_eta) * cp_fcn(curr_eta_k) / (cp_fcn(0) * 
                                                            cp_fcn(sqrt_term))
      rate_vec <- exp(curr_deg) * exp(curr_beta) * cp_term
      ## just need to change this
      llh[sim, i,] <- mapply(dpois, x = y_sim[i,],
                             lambda = rate_vec, log = TRUE)
    }
  }
  llh
}



llh_sim <- llh_y(all_post, n, k, num_sims, ls.dim = 3, y_sim)
## can then just modify this if I know the right way it does this...

dim(llh_sim)
mat <- apply(llh_sim, 1, function(x) as.vector(x)) |> t() 
## this flattens row by row
dim(mat)


loo_2015 <- loo(mat)
loo_2015


loo_compare(loo_2015, loo_est, loo_est_2010)

## should probably try to improve the mcmc for the 2015 model here for a 
## full comparison, but initial starting point

num_pars <- 1035
par_rhat <- rep(NA, num_pars)

for(i in 1:num_pars){
  draws <- out$sims[ , , i]
  par_rhat[i] <- rhat(draws)
}

summary(par_rhat)