#### October 2nd 2024 ####

## the goal of this script is to simulate data from the 
## latent position model of McCormick and Zheng, 2015
## with positions on the p+1 dimensional hypersphere


## then fit the 2006, 2010, 2015 and 2019 models to this simulated data 
## and examine the performance in terms of the estimated degree and 
## ppc for this data
## the latent position vectors drawn uniformly from the surface
## and that p = 2



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
source(here("Summer_2024/", "helper_latent_surface_model.R"))
source(here("Summer_2024/", "helper_functions_latent_ard.R"))

# Simulate the Positions and Overall Network -----------------------------------

n_population <- 1e5
n_subpop <- 15   # number of subpops observed
n_sample <- 1000 # the people where ARD is recorded

perc_subpop <- 0.01  # prevalence of subpop in population, same for all subpops
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

lambda <- 8 ## scaling parameter in the distance
for(k in 1:n_subpop){
  
  ## randomly assign nodes to subpopulation based on distance
  ## to subpop center
  curr_center <- subpop_centers[k, ]
  
  dot_prod <- latent_positions %*% curr_center
  dist <- acos(dot_prod)
  ## these distances between 0 and 2pi
  prob_sub <- exp(- lambda * dist)
  prob_sub <- prob_sub /sum(prob_sub) * n_population * perc_subpop
  prob_sub <- ifelse(prob_sub > 1, 1, prob_sub)
  sub_member <- rbinom(n_population, 1, p = prob_sub)
  cat(paste0(sum(sub_member), "\n"))  
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
gamma <- -3 ## smaller values lead to larger degree

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
  labs(title = "Positions of Subpopulations in Projected Space",
       x = expression(theta), y = expression(phi)) +
  scale_x_continuous(breaks = c(pi, 2*pi), 
                     labels = c(expression(pi), expression(2 * pi))) +
  scale_y_continuous(breaks = c(pi/2, pi), 
                     labels = c(expression(pi/2), expression(pi)))

plot_subpops



# Fit Simple Null Models --------------------------------------------------

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
ppc_fit_null_1 + 
  labs(title = "Null Erdos Renyi Model") 



## Fit the second null model, allowing varying degrees

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
ppc_fit_null_2 + 
  labs(title = "Null Model, varying Degrees") 

# Fit the 2006 ARD Model --------------------------------------------------

## to do this need to specify a strong prior for the beta terms
## in this model
## need to figure out how to do this, incorporating the known 
## truth about their propensity in the network

### eliminate the zero variance responses ###

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
mu_beta <- log(apply(y_sim, 2, sum)/sum(samp_degree))
# mu_beta <- subpop_greg_means
sd_beta <- rep(subpop_greg_sd, length(mu_beta))  
## take this as true sd of greg
## getting poor convergence here so may need to modify this

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



## Plot histogram of estimated degree a_i = exp(alpha_i)
## and compare it to the true known degrees for this sample

est_degrees_2006 <- stan_fit_simple$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("alpha")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("alpha"),
               names_to = "node",
               values_to = "est") |> 
  mutate(degree = exp(est),
         node = parse_number(node)) 


true_degrees <- tibble(node = 1:nrow(y_sim),
                       true_degree = samp_degree)

true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")

est_degrees |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  labs(x = "Expected Degree", title = "Posterior Estimates",
       subtitle = "2006 Model", y = "") 
## plot the ppc for this in terms of the proportion of y_{ik}=l 
## for different values of l


## construct data to be used for the ppc plot below

ppc_2006 <- construct_ppc(stan_fit_2006, y_sim)

ppc_fit_2006 <- plot_ests(ppc_2006$ppc_draws,
                          ppc_2006$y_tibble,
                          prop_val = 1)
ppc_fit_2006 + 
  labs(title = "Zheng et al, 2006") 

## note that when we do this with 0, no simulated y entries = 0 which messes
## up the plot a bit



# Fit the 2010 Mixing Model -----------------------------------------------


## to fit this model need to specify ego and alter groups 
## in both the nodes in ard sample and also in the ard subpopulations
## also need to specify beta, which corresponds to proportions
## of the alter groups in a subpop

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
                          prop_val = 1)
ppc_fit_2010 + 
  labs(title = "McCormick et al, 2010") 


# Fit the 2019 Model ------------------------------------------------------

## here we do not fit this model exactly but rather modify the 2010
## model to parameterise the rows of the beta matrix with a simple 
## kernel model



# Fit the Latent ARD Model ------------------------------------------------

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


# Do PPC using Latent ARD Model -------------------------------------------

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
  
true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution")
  

summary(est_degrees_2015$degree)
summary(true_degrees$true_degree)


## then simulate the y_ik entries using the draws

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
