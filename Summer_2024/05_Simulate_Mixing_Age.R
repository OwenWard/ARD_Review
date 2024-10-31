
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
library(gridExtra)
library(grid)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())
source(here("Summer_2024/", "helper_model_checking.R"))
source(here("Summer_2024/", "helper_plots.R"))
set.seed(100)
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
                                    "log_d[100]", "inv_omega", 
                                    "rho", "beta"))

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


## plot combined degree distribution from each of these models

comb_deg_plot <- est_degrees_2010 |> 
  select(node, degree = est_degree) |> 
  mutate(model = "McCormick et al. 2010") |> 
  bind_rows(est_degrees_2019 |> select(node, degree = est_degree) |> 
              mutate(model = "Sahai et al. 2019")) |> 
  bind_rows(true_degrees |> 
              select(node, degree = true_degree) |> 
              mutate(model = "True Degree")) |> 
    mutate(model = factor(model, 
                          levels = c("True Degree",
                                     "McCormick et al. 2010",
                                     "Sahai et al. 2019"))) |> 
    ggplot(aes(x = degree, color = model)) +
    geom_density(size = 1, key_glyph = "path") +
    guides(color = guide_legend(override.aes = list(fill = NA)),
           override.aes = list(linetype = 1, size = 1)) +
    labs(color = "", x = "Degree", y = "") +
    scale_x_continuous(expand = c(0, 100), limits = c(0, 5000)) +
    theme_single_legend()

comb_deg_plot
ggsave(filename = here("Summer_2024", "figures",
                       "mix_2_node_10_19_degree_true.png"),
       plot = comb_deg_plot,
       dpi = 600,
       height = 5, width = 7) 

## then do the ppc down here

ppc_2019 <- construct_ppc(stan_fit_2019, y_sim)

ppc_fit_2019 <- plot_ests(ppc_2010$ppc_draws,
                          ppc_2010$y_tibble,
                          prop_val = 1)
ppc_fit_2019 + 
  labs(title = "Sahai et al, 2019") 



# Construct the block of PPC Plots ----------------------------------------

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
                       "mix_2_ppc_all_2010.png"),
       plot = grid_plot,
       dpi = 600,
       height = 5, width = 7)


## 2019 model

ppc_fit_2019_0 <- plot_ests(ppc_2019$ppc_draws,
                            ppc_2019$y_tibble,
                            prop_val = 0)

ppc_0_2019 <- ppc_fit_2019_0 +
  labs(title = expression(Pr(y[ik] == 0)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))

ppc_fit_2019_1 <- plot_ests(ppc_2019$ppc_draws,
                            ppc_2019$y_tibble,
                            prop_val = 1)

ppc_1_2019 <- ppc_fit_2019_1 +
  labs(title = expression(Pr(y[ik] == 1)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))

ppc_fit_2019_3 <- plot_ests(ppc_2019$ppc_draws,
                            ppc_2019$y_tibble,
                            prop_val = 3)

ppc_3_2019 <- ppc_fit_2019_3 +
  labs(title = expression(Pr(y[ik] == 3)),
       subtitle = "") +
  theme_single_grid() +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2))

ppc_fit_2019_5 <- plot_ests(ppc_2019$ppc_draws,
                            ppc_2019$y_tibble,
                            prop_val = 5)

ppc_5_2019 <- ppc_fit_2019_5 +
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
    arrangeGrob(ppc_0_2019, ppc_1_2019,
                ppc_3_2019, ppc_5_2019, nrow = 1, ncol = 4), # Plots in a 2x2 grid
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
                       "mix_2_ppc_all_2019.png"),
       plot = grid_plot,
       dpi = 600,
       height = 5, width = 7)
