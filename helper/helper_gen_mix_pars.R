mixing_sim <- readRDS(here("Summer_2024", "sim_pars", "mixing_sim.RDS"))


n_names <- 14

omega_sim <- runif(n_names, 0.85, 1)

M_sim <- matrix(c(rdirichlet(1, mixing_sim[1, ]),
                  rdirichlet(1, mixing_sim[2, ]),
                  rdirichlet(1, mixing_sim[3, ]),
                  rdirichlet(1, mixing_sim[4, ]),
                  rdirichlet(1, mixing_sim[5, ]),
                  rdirichlet(1, mixing_sim[6, ])),
                nrow = 6, ncol = 8, byrow = TRUE)

## after running source(paste0(code_path,"src/load_data/occs.R"))
## there is then an object called beta_names,
## along with ego_sex_age


## this loads a lot of data
source(here("Summer_2024/", "sim_pars", "surveys.R"))
source(here("Summer_2024/", "sim_pars", "names.R"))

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

degree <- read.csv(here("Summer_2024", "sim_pars", "degree_mix.csv"))

degree_sim <- degree$CombSexAge

degree_sim

sim_pars <- list(n = length(degree_sim),
                 k = n_names,
                 degree = degree_sim,
                 omega = omega_sim,
                 M = M_sim,
                 beta = beta_sim,
                 ego = ego_sim)

saveRDS(sim_pars, file = here("Summer_2024", 
                              "complete_mixing_sim.RDS"))