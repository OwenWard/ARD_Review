#### July 5th 2024 ####

## Use simulated data based on the Sahai et al paper and compare 
## multiple existing models 



# Setup -------------------------------------------------------------------

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



# Simulate Data -----------------------------------------------------------

## here we simulate data based on the real data used in Sahai et al


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

names_data_sim <- simulate_mixing(degree_sim, omega_sim,
                                  M_sim, beta_sim, ego_sim)

names_data_sim


y_sim <- names_data_sim




# Fit Zheng et al 2006 ----------------------------------------------------


stan_file <- here("Summer_2024", "zheng_et_al_06.stan")

model1 <- cmdstan_model(stan_file)


stan_data <- list(N = nrow(y_sim),
                  K = ncol(y_sim),
                  y = y_sim)

fit1 <- model1$sample(data = stan_data,
                       seed = 123,
                       chains = 4,
                       iter_sampling = 500,
                       iter_warmup = 500,
                       parallel_chains = 4,
                       refresh = 100)


## do some initial checking for this model

fit1$summary(variables = c("mu_beta", "sigma_beta", "sigma_alpha", "alpha[1]"))

mcmc_trace(fit1$draws(), 
           pars = c("mu_beta", "sigma_beta", "sigma_alpha", "alpha[1]"))


## do some initial ppc of the values of y

est_out_degrees <- fit1$draws() |> 
  as_draws_df() |> 
  select(starts_with("out")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("out"), values_to = "degree") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+"))) 


true_out_degrees <- tibble(node_id = 1:nrow(y_sim),
                           degree = apply(y_sim, 1, sum))

est_out_degrees |> 
  filter(node_id < 11) |> 
  ggplot(aes(degree)) +
  geom_histogram() +
  facet_wrap(~node_id, nrow = 2) +
  geom_vline(data = true_out_degrees |> filter(node_id < 11),
             aes(xintercept = degree),
             col = "red")


## can we plot estimated gregariousness against truth for this problem?
## can we estimate it under the true data model? don't think so


## recreate figure 10 from zheng et al

y_draws <- fit1$draws() |> 
  as_draws_df()

## first get the generated quantities of interest here

ppc_y <- y_draws |> 
  dplyr::select(starts_with("y_sim")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("y_sim"), values_to = "count") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+")),
         sub_pop_id = str_extract(name, pattern = "\\d+]"),
         sub_pop_id = as.numeric(str_replace(sub_pop_id, "\\]", ""))) 

ppc_y

true_y_mat <- y_sim

matrix_df <- as.data.frame(as.table(true_y_mat))
colnames(matrix_df) <- c("node_id", "sub_pop_id", "count")
matrix_df$node_id <- as.numeric(matrix_df$node_id)
matrix_df$sub_pop_id <- as.numeric(matrix_df$sub_pop_id)
true_y <- as_tibble(matrix_df)

ppc_data_fit1 <- ppc_y

plot_ests <- function(ppc_y, true_y, prop_val = 0) {
  ## takes in two tibbles, ppc_y and true_y, of
  ## the draws of each entry of y and of the true y 
  ## respectively
  
  ppc_prop <- ppc_y |> 
    group_by(sub_pop_id, draw) |> 
    summarise(prop = sum(count == prop_val)/n()) |> 
    group_by(sub_pop_id) |> 
    summarise(avg = mean(prop), 
              lower = quantile(prop, 0.025),
              upper = quantile(prop, 0.975))
  
  true_prop <- true_y |> 
    group_by(sub_pop_id) |> 
    summarise(true_prop = sum(count == prop_val)/n()) 
  
  
  final_plot <- ppc_prop |> 
    left_join(true_prop) |> 
    arrange(true_prop) |> 
    mutate(index = row_number()) |> 
    ggplot(aes(x = true_prop, y = avg)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.001, alpha = 0.5) +
    coord_flip() + 
    geom_abline(slope = 1, intercept = 0, alpha = 0.25, col = "red") +
    labs(x = "Simulated", y = "Data",
         subtitle = paste0("Prop = ", prop_val))
  
  final_plot
}


ppc_fit1 <- plot_ests(ppc_data_fit1, true_y, prop_val = 1)
ppc_fit1 + 
  labs(title = "Zheng et al, 2006")


## can we do any other checks here, given that this is such a simple model?






# Fit McCormick et al 2010 ------------------------------------------------


## the code in the swupnil paper to fit this model assumes the beta 
## matrix is known, which is easier
## it also assumes the egos groups are known, which I think is ok

k <- 8 ## the number of names used in the fitting
mcmc_data_sim_k <- list(E = 6, A = 8, K = k,
                        N = nrow(names_data_sim),
                        y = y_sim[, 1:k],
                        ego = ego_sim,
                        Beta = beta_sim[, c(1:k)],
                        theta_d = c(6.2, .5),
                        theta_o = c(3, 2),
                        alpha = rep(1, 8)/8,
                        p = 1)


degree_mix_mod <- cmdstan_model(stan_file = here("Summer_2024",
                                                 "Degree_Mixing.stan"))

fit2 <- degree_mix_mod$sample(data = mcmc_data_sim_k,
                             iter_warmup = 500,
                             iter_sampling = 500)


y_draws_2 <- fit2$draws() |> 
  as_draws_df()

## first get the generated quantities of interest here

ppc_data_fit2 <- y_draws_2 |> 
  dplyr::select(starts_with("y_sim")) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("y_sim"), values_to = "count") |> 
  mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+")),
         sub_pop_id = str_extract(name, pattern = "\\d+]"),
         sub_pop_id = as.numeric(str_replace(sub_pop_id, "\\]", ""))) 



ppc_fit2 <- plot_ests(ppc_data_fit2, true_y, prop_val = 1)
ppc_fit2 + 
  labs(title = "McCormick et al, 2010")



# Fit Sahai et al 2019 ----------------------------------------------------


