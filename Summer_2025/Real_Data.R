##### October 2025

## Fitting with real McCarthy data as a comparison

dat <- readRDS("Summer_2025/mccartyCount_dta.rds")
head(dat)
dim(dat)



###
### using the McCarty.RData file 
### which I got somewhere...

data("McCarty")

pop_size <- 2.8e+08
n_population <- pop_size
mccarty_props <- McCarty$known / pop_size

dim(dat)


######

mccartyCount_dta <- as.matrix(dat)
# set maximum ARD item responses as lower bounds for network size:
maximum_xiks <-
  apply(mccartyCount_dta, 1, max)

mccartyCount_MaltielRDM_StanData <-
  list(
    N = nrow(mccartyCount_dta),
    K = ncol(mccartyCount_dta),
    y = mccartyCount_dta,
    m = mccarty_props,
    L = maximum_xiks
  )




#######

## then feed

library(cmdstanr)
library(here)
library(bayesplot)
library(posterior)


G1_ind <- 1:29
Pg1 <- as.numeric(mccarty_props)

stan_data_null <- list(N = nrow(dat),
                       K = ncol(dat),
                       y = as.matrix(dat),
                       n_known = length(G1_ind),
                       idx = G1_ind,
                       known_prev = sum(Pg1))

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


subpop_info <- tibble(subpop = 1:29,
                      size = McCarty$known)

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
  filter(subpop < 10) |> 
  ggplot(aes(subpop_size)) +
  geom_histogram() +
  facet_wrap(~subpop, scales = "free", nrow = 3) +
  geom_vline(data = subpop_info |> filter(subpop < 10), aes(xintercept = size), col = "red")
