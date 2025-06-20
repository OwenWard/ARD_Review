#### April 29th 2025

## quick code to check if loo-cv works when simulate data from the 
## model that is fit
## still getting some issues, even though model fits well and
## recovers parameters


K <- 15
n <- 1000
set.seed(101)
n_pop <- n_population
b_s <- rbeta(K, shape1 = 1, shape2 = 5)
b_s
true_pop <- n_pop * b_s

log_d <- rnorm(1, sd = 5)

Y <- matrix(NA, nrow = n, ncol = K)

for(i in 1:n){
  for(k in 1:K){
    Y[i,k] <- rpois(1, lambda = exp(log_d)* b_s[k])
  }
}

Y

known_pops <- 1:5
G1_ind <- known_pops
stan_data_simple <- list(N = nrow(Y),
                       K = ncol(Y),
                       y = Y,
                       n_known = length(G1_ind),
                       idx = G1_ind,
                       known_prev = sum(b_s[known_pops]))

stan_file_null_01 <- here("stan_models", "null_model_01_scaled.stan")

mod_null_01 <- cmdstan_model(stan_file_null_01)

stan_fit_01 <- mod_null_01$sample(data = stan_data_simple,
                                  seed = 123,
                                  chains = 4,
                                  iter_sampling = 1000,
                                  iter_warmup = 1000,
                                  parallel_chains = 4,
                                  refresh = 100)

stan_fit_01$summary(variables = c("scaled_log_d", "scaled_beta"))
loo_1_check <- stan_fit_01$loo(cores = 4)
loo_1_check
