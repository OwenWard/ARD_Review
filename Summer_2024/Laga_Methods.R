### Feb 3rd 2025


# Playing around with package code ----------------------------------------


install.packages("networkscaleup")


library(networkscaleup)

# Analyze an example ard data set using the killworth function
data(example_data)

ard <- example_data$ard
subpop_sizes <- example_data$subpop_sizes
N <- example_data$N

mle.est <- killworth(ard,
  known_sizes = subpop_sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  N = N, model = "MLE"
)

pimle.est <- killworth(ard,
  known_sizes = subpop_sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  N = N, model = "PIMLE"
)

## Compare estimates with the truth
plot(mle.est$degrees, example_data$degrees)

data.frame(
  true = subpop_sizes[c(3, 5)],
  mle = mle.est$sizes,
  pimle = pimle.est$sizes
)

ard <- example_data$ard
subpop_sizes <- example_data$subpop_sizes
known_ind <- c(1, 2, 4)
N <- example_data$N
overdisp.est <- overdispersedStan(ard,
  known_sizes = subpop_sizes[known_ind],
  known_ind = known_ind,
  G1_ind = 1,
  G2_ind = 2,
  B2_ind = 4,
  N = N,
  chains = 1,
  cores = 1,
  warmup = 1000,
  iter = 2000
)
# Compare size estimates
round(data.frame(
  true = subpop_sizes,
  basic = colMeans(overdisp.est$sizes)
))
# Compare degree estimates
plot(example_data$degrees, colMeans(overdisp.est$degrees))
# Look at overdispersion parameter
colMeans(overdisp.est$omegas)



## fitting the raw stan model




library(cmdstanr)
library(tidyverse)
library(bayesplot)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())

# temp_stan_file <- write_stan_file(model_code)

model_laga <- cmdstan_model(stan_file = "Summer_2024/laga_stan/overdispersed.stan")

fit <- model_laga$sample(data = list(n_i = nrow(ard),
                              n_k = ncol(ard),
                              y = ard))

fit$summary() |> 
  print(n = 40)


mcmc_hist(fit$draws(variables = "betas"))
mcmc_trace(fit$draws(variables = "betas"))




# Fitting Laga Version to our Simulated Data ------------------------------

stan_data_null <- list(n_i = nrow(y_sim),
                       n_k = ncol(y_sim),
                       y = y_sim)


model1 <- cmdstan_model(stan_file = "Summer_2024/laga_stan/overdispersed.stan")

stan_fit_model1 <- model1$sample(data = stan_data_null,
                                 seed = 123,
                                 chains = 4,
                                 iter_sampling = 1000,
                                 iter_warmup = 1000,
                                 parallel_chains = 4,
                                 refresh = 100)

stan_fit_model1$summary(variables = "betas")
mcmc_trace(stan_fit_model1$draws(variables = "betas"))
## then rescale this after the fact

## need to convert these properly.. 
fit_draws <- stan_fit_model1$draws() |> 
  as_draws_df()

# draws <- stan_fit_model1$draws() |> as_draws()
betas <- fit_draws |> select(starts_with("betas"))
mu_beta <- fit_draws |> select(starts_with("mu_beta")) |> 
  pull(mu_beta)
alphas <- fit_draws |> select(starts_with("alphas"))
mu_alpha <- rep(NA, nrow(betas))

## how does this work if G1_ind is NULL?
G1_ind <- 1 ## setting this just to get the code working
N <- n_population
known_sizes <- true_subpop_size[G1_ind]
known_prevalences <- known_sizes / N
prevalences_vec <- rep(NA, n_subpop)
prevalences_vec[known_ind] <- known_prevalences
if (!is.null(G1_ind)) {
  Pg1 <- sum(prevalences_vec[G1_ind])
}

for (ind in 1:nrow(betas)) {
    C1 = log(sum(exp(betas[ind, G1_ind])/Pg1))
    C = C1
    alphas[ind, ] = alphas[ind, ] + C
    mu_alpha[ind] = C
    betas[ind, ] = betas[ind, ] - C
    mu_beta[ind] = mu_beta[ind] - C
  }

## then compare these estimates directly, use them to predict stuff, do
## ppcheck also

rescale_betas <- betas
rescale_alphas <- alphas
rescale_mu_beta <- mu_beta
rescale_mu_alpha <- mu_alpha
rescale_degrees <- exp(rescale_alphas)
rescale_sizes <- exp(rescale_betas) * n_population


true_degrees <- tibble(node = 1:nrow(y_sim),
                       true_degree = samp_degree)


est_degrees_2006 <- as_tibble(rescale_degrees) |> 
  mutate(draw = row_number()) |> 
  pivot_longer(cols = starts_with("alpha"),
               names_to = "node",
               values_to = "est") |> 
  mutate(node = parse_number(node))

est_degrees_2006 |> 
  ggplot(aes(est)) +
  geom_histogram() +
  labs(title = "2006 Model") +
  theme_single()

true_degrees |> 
  ggplot(aes(true_degree)) +
  geom_histogram() +
  labs(title = "True Degree Distribution") +
  labs(x = "Node Degree", y = "") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  theme_single()



est_degrees_2006 |> 
  rename(degree = est) |> 
  select(node, degree) |> 
  mutate(model = "Zheng et al, 2006") |> 
  bind_rows(true_degrees |> 
              rename(degree = true_degree) |> 
              mutate(model = "True Degree")) |> 
  ggplot(aes(x = degree, color = model)) +
  geom_density(size = 1) +
  guides(color = guide_legend(override.aes = list(fill = NA)),
         override.aes = list(linetype = 1, size = 1)) +
  labs(color = "", x = "Degree", y = "") +
  theme_single_legend()


## does this look better than the previous approach 
## seems to look quite similar here