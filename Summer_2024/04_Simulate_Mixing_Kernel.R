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
library(gtools)
library(gridExtra)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())

source(here("Summer_2024/", "helper_model_checking.R"))
source(here("Summer_2024/", "helper_plots.R"))
set.seed(100)

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


## plot the degree distribution

true_degrees <- tibble(node = 1:nrow(y_sim),
                       true_degree = sim_pars$degree)

deg_plot <- true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  # labs(title = "True Degree Distribution") +
  labs(x = "Node Degree", y = "") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  theme_single()

ggsave(filename = here("Summer_2024", "figures",
                       "mix_node_degree.png"),
       plot = deg_plot,
       dpi = 600,
       height = 5, width = 7)  


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
                       true_degree = true_degree)

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
                            prop_val = 1)

ppc_1_null_1 <- ppc_fit_null_1 + 
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()

ggsave(filename = here("Summer_2024", "figures",
                       "mix_1_ppc_null_1.png"),
       plot = ppc_1_null_1,
       dpi = 600,
       height = 5, width = 7)

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


## show the null model degree distributions
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
                   degree = true_degree) |> 
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
                     limits = c(0, 1000)) +
  theme_single_legend() 

ggsave(filename = here("Summer_2024", "figures",
                       "mix_1_null_deg.png"),
       plot = null_deg_plot,
       dpi = 600,
       height = 5, width = 7)

ppc_null_2 <- construct_ppc(stan_fit_null_02, y_sim)

ppc_fit_null_2 <- plot_ests(ppc_null_2$ppc_draws, 
                            ppc_null_2$y_tibble,
                            prop_val = 1)

ppc_1_null_2 <- ppc_fit_null_2 +
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()

ggsave(filename = here("Summer_2024", "figures",
                       "mix_1_ppc_null_2.png"),
       plot = ppc_1_null_2,
       dpi = 600,
       height = 5, width = 7)


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

ppc_1_2006 <- ppc_fit_2006 +
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()

ggsave(filename = here("Summer_2024", "figures",
                       "mix_1_ppc_2006.png"),
       plot = ppc_1_2006,
       dpi = 600,
       height = 5, width = 7)

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

## plot 2006 and 2010 degree estimates and the truth in one plot
comb_deg_plot <- est_degrees_2006 |> 
  select(node, degree) |> 
  mutate(model = "Zheng et al, 2006") |> 
  bind_rows(est_degrees_2010 |> select(node, degree = est_degree) |> 
              mutate(model = "McCormick et al, 2010")) |> 
    bind_rows(true_degrees |> 
                rename(degree = true_degree) |> 
                mutate(model = "True Degree")) |> 
    mutate(model = factor(model, 
                          levels = c("True Degree", "Zheng et al, 2006",
                                     "McCormick et al, 2010"))) |> 
    ggplot(aes(x = degree, color = model)) +
    geom_density(size = 1, key_glyph = "path") +
    guides(color = guide_legend(override.aes = list(fill = NA)),
           override.aes = list(linetype = 1, size = 1)) +
    labs(color = "", x = "Degree", y = "") +
    scale_x_continuous(expand = c(0, 5)) +
    theme_single_legend()

comb_deg_plot_
## then save this
ggsave(filename = here("Summer_2024", "figures",
                       "mix_1_node_06_10_degree_true.png"),
       plot = comb_deg_plot,
       dpi = 600,
       height = 5, width = 7) 

## then do ppc for this model also
ppc_2010 <- construct_ppc(stan_fit_2010, y_sim)

ppc_fit_2010 <- plot_ests(ppc_2010$ppc_draws,
                          ppc_2010$y_tibble,
                          prop_val = 1)

ppc_1_2010 <- ppc_fit_2010 +
  labs(title = expression("Proportion of " * y[ik] == 1),
       subtitle = "") +
  theme_single()

ggsave(filename = here("Summer_2024", "figures",
                       "mix_1_ppc_2010.png"),
       plot = ppc_1_2010,
       dpi = 600,
       height = 5, width = 7)



# Creating a grid of PPC plots --------------------------------------------

## first the 2006 model

ppc_fit_2006_0 <- plot_ests(ppc_2006$ppc_draws,
                            ppc_2006$y_tibble,
                            prop_val = 0)

ppc_0_2006 <- ppc_fit_2006_0 +
  labs(title = expression(Pr(y[ik] == 0)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))

ppc_fit_2006_1 <- plot_ests(ppc_2006$ppc_draws,
                            ppc_2006$y_tibble,
                            prop_val = 1)

ppc_1_2006 <- ppc_fit_2006_1 +
  labs(title = expression(Pr(y[ik] == 1)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))

ppc_fit_2006_3 <- plot_ests(ppc_2006$ppc_draws,
                            ppc_2006$y_tibble,
                            prop_val = 3)

ppc_3_2006 <- ppc_fit_2006_3 +
  labs(title = expression(Pr(y[ik] == 3)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))

ppc_fit_2006_5 <- plot_ests(ppc_2006$ppc_draws,
                            ppc_2006$y_tibble,
                            prop_val = 5)

ppc_5_2006 <- ppc_fit_2006_5 +
  labs(title = expression(Pr(y[ik] == 5)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))


x_label <- textGrob("Simulated",
                    gp = gpar(fontsize = 16))
y_label <- textGrob("Data",
                    rot = 90, gp = gpar(fontsize = 16))

grid_plot <- grid.arrange(
  arrangeGrob(
    y_label,                # Left y-axis label
    arrangeGrob(ppc_0_2006, ppc_1_2006,
                ppc_3_2006, ppc_5_2006, nrow = 1, ncol = 4), # Plots in a 2x2 grid
    ncol = 2,
    widths = unit.c(grobWidth(y_label) + unit(0.5, "line"),
                    unit(0.95, "npc") - grobWidth(y_label))
  ),
  x_label,                  # Bottom x-axis label
  nrow = 2,
  heights = unit.c(unit(0.65, "npc") - 
                     grobHeight(x_label),
                   grobHeight(x_label) - unit(0.5, "line"))
)

ggsave(filename = here("Summer_2024", "figures",
                       "mix_1_ppc_all_2006.png"),
       plot = grid_plot,
       dpi = 600,
       height = 5, width = 7)


## then the 2010 model

ppc_fit_2010_0 <- plot_ests(ppc_2010$ppc_draws,
                            ppc_2010$y_tibble,
                            prop_val = 0)

ppc_0_2010 <- ppc_fit_2010_0 +
  labs(title = expression(Pr(y[ik] == 0)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))

ppc_fit_2010_1 <- plot_ests(ppc_2010$ppc_draws,
                            ppc_2010$y_tibble,
                            prop_val = 1)

ppc_1_2010 <- ppc_fit_2010_1 +
  labs(title = expression(Pr(y[ik] == 1)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))

ppc_fit_2010_3 <- plot_ests(ppc_2010$ppc_draws,
                            ppc_2010$y_tibble,
                            prop_val = 3)

ppc_3_2010 <- ppc_fit_2010_3 +
  labs(title = expression(Pr(y[ik] == 3)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))

ppc_fit_2010_5 <- plot_ests(ppc_2010$ppc_draws,
                            ppc_2010$y_tibble,
                            prop_val = 5)

ppc_5_2010 <- ppc_fit_2010_5 +
  labs(title = expression(Pr(y[ik] == 5)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))



x_label <- textGrob("Simulated",
                    gp = gpar(fontsize = 16))
y_label <- textGrob("Data",
                    rot = 90, gp = gpar(fontsize = 16))

grid_plot <- grid.arrange(
  arrangeGrob(
    y_label,                # Left y-axis label
    arrangeGrob(ppc_0_2010, ppc_1_2010,
                ppc_3_2010, ppc_5_2010, nrow = 1, ncol = 4), # Plots in a 2x2 grid
    ncol = 2,
    widths = unit.c(grobWidth(y_label) + unit(0.5, "line"),
                    unit(0.95, "npc") - grobWidth(y_label))
  ),
  x_label,                  # Bottom x-axis label
  nrow = 2,
  heights = unit.c(unit(0.65, "npc") - 
                     grobHeight(x_label),
                   grobHeight(x_label) - unit(0.5, "line"))
)

ggsave(filename = here("Summer_2024", "figures",
                       "mix_1_ppc_all_2010.png"),
       plot = grid_plot,
       dpi = 600,
       height = 5, width = 7)
