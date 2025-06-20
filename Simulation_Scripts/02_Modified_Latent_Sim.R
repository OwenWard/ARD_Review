library(rotasym)
library(tidyverse)
library(tidymodels) # to use coord_obs_pred
library(cmdstanr)
library(here)
library(bayesplot)
library(posterior)
library(grid)
library(gridExtra)
library(ggdist)
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


# Read in Data, Plot Model information ------------------------------------



sim_data <- readRDS(file = here("Summer_2025", "ard_latent_mod.RDS"))
y <- sim_data$y_sim
G1_ind <- sim_data$G1_ind
known_prev <- sum(sim_data$true_subpops[G1_ind])/sim_data$n_population
true_subpops <- sim_data$true_subpops
sim_data$n_subpop
Pg1 <- sum(true_subpops[G1_ind]/n_population)


true_alpha <- sim_data$true_alpha

true_degree <- tibble(alpha = true_alpha) |> 
  mutate(degree = exp(alpha)) |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  theme_single() +
  labs(y = element_blank(), x = "Sample Degree")

ggsave(filename = here("Summer_2025", "figures",
                       "latent_true_degree.png"),
       plot = true_degree,
       dpi = 600,
       height = 5, width = 5)

# Fit first null model ----------------------------------------------------



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

color_scheme_set("viridis")
trace_plot <- mcmc_trace(stan_fit_null_01$draws(), pars = "scaled_log_d") +
  labs(y = "Scaled Log Degree", x = "MCMC Iterations",
       title = "Erdos Renyi Model, Scaled")


ggsave(filename = here("Summer_2025", "figures",
                       "latent_trace.png"),
       plot = trace_plot,
       dpi = 600,
       height = 5, width = 5)

## Quickly fit unscaled model

stan_file_null_01_unscaled <- here("stan_models", "null_model_01.stan")
mod_null_01_unscaled <- cmdstan_model(stan_file_null_01_unscaled)

stan_fit_null_01_unscaled <- mod_null_01_unscaled$sample(data = stan_data_null,
                                       seed = 123,
                                       chains = 4,
                                       iter_sampling = 1000,
                                       iter_warmup = 1000,
                                       parallel_chains = 4,
                                       refresh = 100)

trace_plot2 <- mcmc_trace(stan_fit_null_01_unscaled$draws(), 
                          pars = "log_d") +
  labs(y = "Log Degree", x = "MCMC Iterations",
       title = "Erdos Renyi Model, No Scaling")

ggsave(filename = here("Summer_2025", "figures",
                       "latent_trace_unscaled.png"),
       plot = trace_plot2,
       dpi = 600,
       height = 5, width = 5)

##

samp_degree <- exp(alpha)

stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_log_d")) |> 
  mutate(draw = row_number())  |> 
  mutate(degree = exp(scaled_log_d)) |> 
  ggplot(aes(degree)) +
  geom_histogram()

hist(samp_degree)

subpop_info <- tibble(subpop = 1:sim_data$n_subpop,
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

ppc_1_plot <- with(ppc_null_1, plot_ests_all(ppc_null_1$ppc_draws,
                                             ppc_null_1$y_tibble, 
                                             prop_vals =  c(0, 1, 3, 5, 10)))

# (ppc_fit_null_1 <- plot_ests(ppc_null_1$ppc_draws, 
#                              ppc_null_1$y_tibble,
#                              prop_val = 1))

rm(ppc_null_1)
# ppc_1_null_1 <- ppc_fit_null_1 + 
#   labs(title = expression("Proportion of " * y[ik] == 1),
#        subtitle = "") +
#   theme_single()
# ppc_1_null_1
# 
# ppc_plot_null_1 <- plot_ests_all(ppc_y = ppc_null_1$ppc_draws,
#               true_y = ppc_null_1$y_tibble,
#               prop_vals = c(0, 1, 3, 5, 10))


# y_vals <- c(0, 1, 3, 5, 10)
# ppc_fit_null_1_0 <- plot_ests(ppc_null_1$ppc_draws, 
#                             ppc_null_1$y_tibble,
#                             prop_val = y_vals[1]) + 
#   labs(title = expression(y[ik] == 0),
#        x = element_blank(), y = element_blank(),
#        subtitle = "") + theme_single()
# ppc_fit_null_1_1 <- plot_ests(ppc_null_1$ppc_draws, 
#                               ppc_null_1$y_tibble,
#                               prop_val = y_vals[2]) + 
#   labs(title = expression(y[ik] == 1),
#        x = element_blank(), y = element_blank(),
#        subtitle = "") + theme_single()
# ppc_fit_null_1_3 <- plot_ests(ppc_null_1$ppc_draws, 
#                               ppc_null_1$y_tibble,
#                               prop_val = y_vals[3]) + 
#   labs(title = expression(y[ik] == 3),
#        x = element_blank(), y = element_blank(),
#        subtitle = "") + theme_single()
# ppc_fit_null_1_5 <- plot_ests(ppc_null_1$ppc_draws, 
#                               ppc_null_1$y_tibble,
#                               prop_val = y_vals[4]) + 
#   labs(title = expression(y[ik] == 5),
#        x = element_blank(), y = element_blank(),
#        subtitle = "") + theme_single()
# ppc_fit_null_1_10 <- plot_ests(ppc_null_1$ppc_draws, 
#                               ppc_null_1$y_tibble,
#                               prop_val = y_vals[5]) + 
#   labs(title = expression(y[ik] == 10),
#        x = element_blank(), y = element_blank(),
#        subtitle = "") + theme_single()
# 
# ppc_fit_null_1_10
# 
# grid.arrange(ppc_fit_null_1_0, 
#              ppc_fit_null_1_1,
#              ppc_fit_null_1_3,
#              ppc_fit_null_1_5,
#              # ppc_fit_null_1_10,
#              nrow = 1)


ggsave(filename = here("Summer_2025", "figures",
                       "latent_ppc_1.png"),
       plot = ppc_plot_null_1,
       dpi = 600,
       height = 5, width = 10)

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

# (ppc_fit_null_2 <- plot_ests(ppc_null_2$ppc_draws, 
#                              ppc_null_2$y_tibble,
#                              prop_val = 1))
# rm(ppc_null_2)
# ppc_1_null_2 <- ppc_fit_null_2 +
#   labs(title = expression("Proportion of " * y[ik] == 1),
#        subtitle = "") +
#   theme_single()
# 
# ppc_1_null_2

ppc_plot_null_2 <- plot_ests_all(ppc_y = ppc_null_2$ppc_draws,
                                 true_y = ppc_null_2$y_tibble,
                                 prop_vals = c(0, 1, 3, 5, 10))
rm(ppc_null_2)

ggsave(filename = here("Summer_2025", "figures",
                       "latent_ppc_2.png"),
       plot = ppc_plot_null_2,
       dpi = 600,
       height = 5, width = 10)


## show the estimated degree from null models and compare to true distribution
null_1_deg <- stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_log_d")) |> 
  mutate(draw = row_number())  |> 
  mutate(degree = exp(scaled_log_d)) |> 
  select(draw, degree) |> 
  mutate(model = "Null 1")

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
  mutate(model = "Null 2")

true_deg <- tibble(draw = 1,
                   degree = exp(true_alpha)) |> 
  mutate(model = "True")


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
                     limits = c(0, 800)) +
  scale_color_manual(values = color_vec) +
  theme_single_legend() 

null_deg_plot

ggsave(filename = here("Summer_2025", "figures",
                       "latent_null_degree.png"),
       plot = null_deg_plot,
       dpi = 600,
       height = 5, width = 5)


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

# stan_fit_zheng$summary(variables = c("sigma_alpha", "scaled_beta"))

## check it recovers the groups correctly
b_draws <- stan_fit_zheng$draws() |> 
  as_draws_df() |> 
  select(starts_with("scaled_beta"))

size_ests <- exp(b_draws) * n_population
head(rowSums(size_ests[, G1_ind]))
sum(true_subpops[1:14])


# est_degrees_2006 <- stan_fit_zheng$draws() |> 
#   as_draws_df() |> 
#   dplyr::select(starts_with("scaled_alpha")) |> 
#   mutate(draw = row_number()) |> 
#   pivot_longer(cols = starts_with("scaled_alpha"),
#                names_to = "node",
#                values_to = "log_estimate") |> 
#   mutate(node = parse_number(node)) |> 
#   mutate(est_degree = exp(log_estimate)) 

## look at posterior estimates of degree and subpop size
# stan_fit_zheng$draws() |> 
#   as_draws_df() |> 
#   dplyr::select(starts_with("scaled_alpha")) |> 
#   mutate(draw = row_number())  |> 
#   pivot_longer(cols = starts_with("scaled_alpha"), 
#                names_to = "par", values_to = "sample") |> 
#   mutate(node = parse_number(par)) |> 
#   mutate(degree = exp(sample)) |> 
#   ggplot(aes(degree)) +
#   geom_histogram()



# stan_fit_zheng$draws() |> 
#   as_draws_df() |> 
#   dplyr::select(starts_with("scaled_beta[")) |>
#   mutate(draw = row_number()) |> 
#   pivot_longer(cols = starts_with("scaled_beta"),
#                names_to = "par", 
#                values_to = "sample") |> 
#   mutate(subpop = parse_number(par),
#          subpop_size = n_population * exp(sample)) |> 
#   # filter(!(subpop %in% G1_ind)) |> 
#   ggplot(aes(subpop_size)) +
#   geom_histogram() +
#   facet_wrap(~subpop, scales = "free", nrow = 3) +
#   geom_vline(data = subpop_info, aes(xintercept = size), col = "red")



ppc_zheng <- construct_ppc(stan_fit_zheng, y_sim)

# (ppc_fit_zheng <- plot_ests(ppc_zheng$ppc_draws, 
#                             ppc_zheng$y_tibble,
#                             prop_val = 1))
# 
# ppc_plot_zheng <- ppc_fit_zheng + 
#   labs(title = expression("Proportion of " * y[ik] == 1),
#        subtitle = "") +
#   theme_single()

ppc_plot_zheng <- plot_ests_all(ppc_y = ppc_zheng$ppc_draws,
                                 true_y = ppc_zheng$y_tibble,
                                 prop_vals = c(0, 1, 3, 5, 10))

rm(ppc_zheng)
# ggsave(filename = here("Summer_2025", "figures",
#                        "latent_ppc_zheng.png"),
#        plot = ppc_plot_zheng,
#        dpi = 600,
#        height = 5, width = 10)



# Fit McCormick 2015 ------------------------------------------------------

stan_fit_2015 <- readRDS(file = here("stan_models",
                                      "ard_latent_2015_cluster_fit.RDS"))

# est_degrees_2015 <- stan_fit_2015$draws() |> 
#   as_draws_df() |> 
#   dplyr::select(starts_with("scaled_alpha")) |> 
#   mutate(draw = row_number()) |> 
#   pivot_longer(cols = starts_with("scaled_alpha"),
#                names_to = "node",
#                values_to = "log_estimate") |> 
#   mutate(node = parse_number(node)) |> 
#   mutate(est_degree = exp(log_estimate)) 

stan_fit_2015$draws() |> 
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


ppc_2015 <- construct_ppc(stan_fit_2015, y_sim)

ppc_plot_2015 <- plot_ests_all(ppc_y = ppc_2015$ppc_draws,
                                true_y = ppc_2015$y_tibble,
                                prop_vals = c(0, 1, 3, 5, 10))
rm(ppc_2015)

ggsave(filename = here("Summer_2025", "figures",
                       "latent_ppc_2015.png"),
       plot = ppc_plot_2015,
       dpi = 600,
       height = 5, width = 10)


# Additional Plots --------------------------------------------------------

RColorBrewer::brewer.pal(n = 4, name = "Set1")


gg_color_hue <- function(n) { 
  
  
  hues = seq(15, 375, length = n + 1) 
  
  
  hcl(h = hues, l = 65, c = 100)[1:n]
  
  
}

consistent_colours <- tibble(
  model = c("Null 1",
            "Null 2",
            "2006", 
            "2015", 
            "True"),
  colours = RColorBrewer::brewer.pal(n = 5, name = "Set1")
)

color_vec <- setNames(consistent_colours$colours, consistent_colours$model)


## Degree estimates from 2006 and 2015, alongside true degrees

deg_plot <- bind_rows(est_degrees_2006 |> 
                        rename(degree = est_degree) |> 
                        mutate(model = "2006"),
                      est_degrees_2015 |> 
                        rename(degree = est_degree) |> 
                        mutate(model = "2015"),
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
                     limits = c(0, 800),
                     breaks = scales::pretty_breaks(n = 3)) +
  scale_color_manual(values = color_vec) +
  theme_single_legend() 

deg_plot

ggsave(filename = here("Summer_2025", "figures",
                       "latent_est_degree.png"),
       plot = deg_plot,
       dpi = 600,
       height = 5, width = 5)


## Subpopulation estimates for unknown population for all 4 models

subpop_null_1 <- stan_fit_null_01$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("b[")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("b"),
               names_to = "par", 
               values_to = "sample") |> 
  mutate(subpop = parse_number(par),
         subpop_size = n_population * sample)  |> 
  filter(subpop == 15) |> 
  mutate(model = "Null 1")

subpop_null_2 <- stan_fit_null_02$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("b[")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("b"),
               names_to = "par", 
               values_to = "sample") |> 
  mutate(subpop = parse_number(par),
         subpop_size = n_population * sample)  |> 
  filter(subpop == 15) |> 
  mutate(model = "Null 2")

subpop_2006 <- stan_fit_zheng$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_beta[")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("scaled_beta"),
               names_to = "par", 
               values_to = "sample") |> 
  mutate(subpop = parse_number(par),
         subpop_size = n_population * exp(sample)) |> 
  filter(subpop == 15) |> 
  mutate(model = "2006")

subpop_2015 <- stan_fit_2015$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_beta[")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("scaled_beta"),
               names_to = "par", 
               values_to = "sample") |> 
  mutate(subpop = parse_number(par),
         subpop_size = n_population * exp(sample)) |> 
  filter(subpop == 15) |> 
  mutate(model = "2015")


subpop_data <- bind_rows(subpop_null_1, 
                         subpop_null_2,
                         subpop_2006,
                         subpop_2015) |> 
  mutate(model = factor(model, 
                        levels = c("Null 1", "Null 2",
                                   "2006", "2015")))

subpop_plot <- subpop_data |> 
  ggplot(aes(x = subpop_size, y = as.factor(subpop), color = model)) +
  stat_pointinterval(position = position_dodge(width = 0.2)) +
  labs(y = element_blank(), x = "Unknown Subpopulation Size", 
       color = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_color_manual(values = color_vec) +
  geom_vline(xintercept = true_subpops[15], col = "red", linetype = 2)
subpop_plot

ggsave(filename = here("Summer_2025", "figures",
                       "latent_subpop_est.png"),
       plot = subpop_plot,
       dpi = 600,
       height = 5, width = 5)

## Subpopulation estimates from 2015 model for all populations

true_vals <- tibble(
  subpop = 1:15,  # or matching your subpop indices
  true_subpop_size = true_subpops  # true values on the same scale as subpop_size
)

stan_fit_2015$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_beta[")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("scaled_beta"),
               names_to = "par", 
               values_to = "sample") |> 
  mutate(subpop = parse_number(par),
         subpop_size = n_population * exp(sample)) |> 
  group_by(subpop) |> 
  summarise(avg = mean(subpop_size))

subpop_y <- true_vals |> 
  arrange(true_subpop_size) |> 
  mutate(y_val = row_number(),
         shape = as.factor(ifelse(subpop == 15, 4, 1)))

subpops_plot <- stan_fit_2015$draws() |> 
  as_draws_df() |> 
  dplyr::select(starts_with("scaled_beta[")) |>
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("scaled_beta"),
               names_to = "par", 
               values_to = "sample") |> 
  mutate(subpop = parse_number(par),
         subpop_size = n_population * exp(sample)) |>
  left_join(subpop_y, by = "subpop") |> 
  ggplot(aes(x = subpop_size, y = y_val, group = subpop)) +
  stat_pointinterval(position = position_dodge(width = 0.2)) +
  # geom_boxplot() +
  geom_point(
    data = subpop_y,
    aes(x = true_subpop_size, #xend = true_subpop_size,
        y = y_val, shape = shape),# yend = y_val + 0.5),
    color = "red"
  ) + 
  labs(x = "Subpopulation Size", y = "Subpopulation") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

subpops_plot
ggsave(filename = here("Summer_2025", "figures",
                       "latent_subpop_all.png"),
       plot = subpops_plot,
       dpi = 600,
       height = 5, width = 5)
  
## then arrange these by largest, add little red lines for true values, etc



# Degree Distribution Plot ------------------------------------------------

degree_null_1 <- tibble(degree = exp(stan_fit_null_01$summary("scaled_log_d") |> 
                                       pull(mean)),
                        model_abbr = "er")
degree_null_2 <- tibble(degree = exp(stan_fit_null_02$summary("scaled_log_d") |> 
                                       pull(mean)),
                        model_abbr = "deg")
degree_od <- tibble(degree = exp(stan_fit_zheng$summary("scaled_alpha") |> 
                                   pull(mean)),
                    model_abbr = "overd")

degree_latent <- tibble(degree = exp(stan_fit_2015$summary("scaled_alpha") |> 
                                       pull(mean)),
                        model_abbr = "latent")

degs_all <- bind_rows(degree_null_1,
                      degree_null_2,
                      degree_od,
                      degree_latent)

hist_input <- data.frame(degree = exp(sim_data$true_alpha)) 
bin_data <- ggplot(hist_input, aes(x = degree)) +
  geom_histogram(bins = 50)

hist_data <- ggplot_build(bin_data)$data[[1]] %>%
  as_tibble() %>%
  mutate(scaled_count = count / max(count)) 

deg_hist_all <- degs_all %>% 
  mutate( model = factor( model_abbr,
                          levels = c("latent","overd","deg","er"),
                          labels = modelnames )) %>%
  filter( model_abbr != "er" ) %>%
  ggplot(aes(degree, col = model)) + 
  # geom_histogram( data = data.frame(degree = exp(sim_data$true_alpha)),
  #                 aes(degree, y = after_stat(density)),
  #                 col="lightgrey",alpha = .1, bins = 50) +
  # geom_histogram(data = hist_data, 
  #                # inherit.aes = FALSE,
  #                aes(x = xmin, y = scaled_count, width = xmax - xmin),
  #                # fill = "lightgrey",
  #          alpha = 0.1, colour = "lightgrey") +
  geom_histogram( data = data.frame(degree = exp(sim_data$true_alpha)),
                  aes(degree,y=after_stat(density)),
                  col="lightgrey",alpha=.1,bins=50)  + 
  # geom_density( data = data.frame(degree=exp(sim_data$true_alpha),
  #                                 model = "Observed"),
  #               aes(degree, y = after_stat(scaled)),
  #               lty = 2, key_glyph = "path" ) +
  geom_density( data = data.frame(degree = exp(sim_data$true_alpha),
                                  model="Observed"),
                aes(degree),
                lty = 2, key_glyph = "path" ) +
  geom_segment(data = degree_null_1,
               aes(x = degree,
                   xend = degree,
                   y=0, yend=Inf,
                   col="Erdos-Renyi")) +
  # geom_vline(xintercept = degree_null_1$degree, col = "Erdos-Renyi") +
  geom_density(key_glyph = "path") + 
  labs(x="Sample Degree", y="") + 
  theme_bw() + 
  scale_color_manual( values=c(gg_color_hue(4),"black")) + 
  guides(color = guide_legend(title = "Model",
                              override.aes = 
                                list(linetype=c(1,1,1,1,2),
                                     col=c(gg_color_hue(4),"black")) ) ) +
  theme(legend.position = "bottom")

final_deg_plot <- deg_hist_all + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(), legend.title = element_blank()) 

ggsave(filename = here("Summer_2025", "figures",
                       "lat_all_deg_plot.png"),
       plot = deg_hist_all,
       dpi = 600,
       height = 7.5, width = 7.5)

# Subpop Size Estimate Plot---------------------------------------------------

size_ests_plot_data_er <- mcmc_intervals_data(stan_fit_null_01$draws("b") |> 
                                                as_draws_df()) |> 
  mutate(parameter = c(rep("Known", 14), "Unknown"))
size_ests_plot_data_vd <- mcmc_intervals_data(stan_fit_null_02$draws("b") |> 
                                                as_draws_df()) |> 
  mutate(parameter = c(rep("Known", 14), "Unknown"))
size_ests_plot_data_od <- mcmc_intervals_data(stan_fit_zheng$draws("scaled_beta") |> 
                                                as_draws_df()) |> 
  mutate(parameter = c(rep("Known", 14), "Unknown")) |> 
  mutate(ll = exp(ll), l = exp(l), m = exp(m), h = exp(h), hh = exp(hh))
size_ests_plot_data_lat <- mcmc_intervals_data(stan_fit_2015$draws("scaled_beta") |> 
                                                as_draws_df()) |> 
  mutate(parameter = c(rep("Known", 14), "Unknown")) |> 
  mutate(ll = exp(ll), l = exp(l), m = exp(m), h = exp(h), hh = exp(hh))

size_ests_plotdata_all <- size_ests_plot_data_lat %>% 
  mutate( model = "Latent Space") %>%
  bind_rows( size_ests_plot_data_er %>% mutate( model = "Erdos Renyi" ),
             size_ests_plot_data_vd %>% mutate( model = "Varying Degree"),
             size_ests_plot_data_od %>% mutate( model = "Overdispersed") ) %>%
  mutate( model = factor( model, levels = modelnames ) )

size_true <- tibble(parameter = c(rep("Known", 14), "Unknown"),
                    true = true_subpops)

size_ests_all_plot <- size_ests_plotdata_all %>% 
  filter(parameter=="Unknown") %>%
  ggplot() +
  geom_segment( aes(x=ll * n_population, xend=hh * n_population,y=model,yend=model),
                linewidth = .5 ) +
  geom_segment( aes(x=l * n_population, xend=h * n_population ,y=model,yend=model),
                linewidth = 2 ) +
  geom_point( aes( x = m * n_population, y = model,
                   col=model), size = 3.5 ) +
  geom_vline(xintercept = size_true %>% 
               filter(parameter=="Unknown") %>% 
               pull(true), lty=2) +
  #scale_y_discrete()
  bayesplot_theme_get() + theme(legend.position = "bottom") +
  guides(col=guide_legend(title="Model")) +
  labs(y = "Model", x = "Size" )

size_ests_all_plot <- size_ests_all_plot + theme_classic() 

ggsave(filename = here("Summer_2025", "figures",
                       "size_ests_all_plot.png"),
       plot = size_ests_all_plot,
       dpi = 600,
       height = 5, width = 5)

# Final PPC Plots ---------------------------------------------------------

modelnames <- c("Latent Space", "Overdispersed", "Varying Degree", "Erdos Renyi")

ppc_all_plotdata <- ppc_plot_2015$final_plot_data |>
  mutate( model = "Latent Space") |>
  bind_rows( ppc_1_plot$final_plot_data |> mutate( model = "Erdos Renyi" ),
             ppc_plot_null_2$final_plot_data |> mutate( model = "Varying Degree" ),
             ppc_plot_zheng$final_plot_data |> mutate( model = "Overdispersed" )) |>
  mutate( model = factor( model, levels = rev(modelnames) ) )

ppc_all_plotdata_lims <- ppc_all_plotdata |>
  group_by(model,count) |>
  summarize( square_lim = max(true_prop,avg),
             ebar_ht = .075*square_lim )

ppc_all_plot <- ppc_all_plotdata |>
  filter( count %in% paste0("P(y_ik =",c(0,1,3,5,10),")") ) |>
  mutate( covers_true = (true_prop >= lower) & (true_prop <= upper)) |>
  left_join(ppc_all_plotdata_lims) |>
  ggplot(aes(y = true_prop, x = avg,
             col = covers_true)) +
  geom_point() +
  geom_point(aes(square_lim,square_lim),alpha=0) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, height = ebar_ht), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.25, lty = 2) + #col = "red") +
  labs(x = "Posterior predictive draws", y = "Observed data",
       #subtitle = paste0("Prop = ", prop_val)
  ) +
  facet_wrap(vars(model,count),scales="free") +
  theme_bw() + theme(legend.position="bottom",strip.background =element_rect(fill="gray95")) + 
  guides(col=guide_legend(title="95% Credible Interval Covers the Truth")) +
  NULL

ppc_all_plot

ggsave(filename = here("Summer_2025", "figures",
                       "latent_all_ppc.png"),
       plot = ppc_all_plot,
       dpi = 600,
       height = 10, width = 10)
