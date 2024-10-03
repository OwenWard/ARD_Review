#### October 2nd 2024 ####

## the goal of this script is to simulate data from the 
## latent position model of McCormick and Zheng, 2015
## with positions on the p+1 dimensional hypersphere


## the latent position vectors drawn uniformly from the surface
## and that p = 2



# Setup and Load Packages -------------------------------------------------

library(rotasym)




# Simulate the Positions and Overall Network -----------------------------------


n_population <- 1e5
n_subpop <- 15   # number of subpops observed
n_sample <- 1000 # the people where ARD is recorded

perc_subpop <- 0.01  # prevalence of subpop in population, same for all subpops
num_subpop <- round(n_population * perc_subpop) # approx size of subpopos in pop

mu <- c(0, 0, 1)
kappa <- 2
latent_positions <- r_vMF(n = n_population, mu = mu, kappa = kappa)


## sample the gregariousness params of ARD sample
## these don't need to be integer necessarily
## this the alpha of the 2006 paper?
samp_greg <- rnorm(n_sample, mean = -2, sd = 1)

## the latent positions of the ARD sample
samp_pos <- latent_positions[sample(1:n_population, n_sample), ]



# Simulate the Subpopulation Members --------------------------------------

## first simulate the subpopulation centers uniformly on the sphere


subpop_centers <- r_vMF(n = n_subpop, mu = c(1, 0, 0), kappa = 1)


## assume that different subpops have common dist for greg
## differences in mean gregariousness and smaller sd than
## the population
subpop_greg_means <- seq(from = -2.5, to = 3.5, length.out = n_subpop)
subpop_greg_sd <- 0.25

y_sim <- matrix(NA, nrow = n_sample, ncol = n_subpop)


gamma <- 1
lambda <- 1

for(k in 1:n_subpop){
  
  ## need to identify which nodes are in this subpopulation
  ## based on distance from the center
  
  curr_center <- subpop_centers[k, ]
  
  dot_prod <- latent_positions %*% curr_center
  dist <- acos(dot_prod)
  ## these distances between 0 and 2pi
  prob_sub <- exp(- lambda * dist)
  prob_sub <- prob_sub /sum(prob_sub) * n_population * perc_subpop
  prob_sub <- ifelse(prob_sub > 1, 1, prob_sub)
  sub_member <- rbinom(n_population, 1, p = prob_sub)
  cat(paste0(sum(sub_member), "\n"))  
  ## believe this number should be the approximate size of the subpopulation
  ## in the population, because want to randomly assign all nodes based on
  ## distance to center
  curr_subpop <- latent_positions[sub_member == 1, ]
  
  subpop_greg <- rnorm(n = sum(sub_member),
                       mean = subpop_greg_means[k],
                       sd = subpop_greg_sd)
  
  ## then need to store these subpopgreg parameters for those members of
  ## the pop
  
  for(i in 1:n_sample){
    curr_greg <- samp_greg[i]
    curr_loc <- samp_pos[i,]
    node_dist <- curr_subpop %*% curr_loc
    exp_term <- curr_greg + subpop_greg + gamma * node_dist
    ## these aren't normalized so need to figure this part out better...
    prob_edge <- 1/(1 + exp(-exp_term)) 
    true_edge <- rbinom(sum(sub_member), size = 1, prob = prob_edge)
    y_sim[i, k] <- sum(true_edge)
  }
  
}


y_sim

apply(y_sim, 2, sum)
## can tweak these parameters as needed to make them look more reasonable

## then can figure out the true alpha and betas from the way the data is
## simulated also
## because the n_sample nodes are a random sample from the population


# Fit the Latent ARD Model ------------------------------------------------


source("Summer_2024/latent_surface_model.R")
source("Summer_2024/helper_functions_latent_ard.R")


dim(y_sim)


ls.dim <- 3
n <- dim(y_sim)[1]
n.iter <- 300#3000
n.thin <- 10
m.iter <- 3
total.prop <- 0.25

## taking this from the github
muk.fix.ind <- sample(1:8, size = 4, replace = F)
muk.fix <- matrix(runif(12), nrow = 4, ncol = 3)
muk.fix <- sweep(muk.fix, MARGIN = 1, 1 / sqrt(rowSums(muk.fix^2)), `*`)


z.pos.init <- generateRandomInitial(n, ls.dim)
out <- f.metro(y_sim,
               total.prop = total.prop,
               n.iter = n.iter,
               # m.iter = m.iter,
               # n.thin = n.thin,
               z.pos.init = z.pos.init,
               muk.fix = muk.fix,
               ls.dim = ls.dim)

posterior <- getPosterior(out, n.iter, m.iter, n.thin, n)
est.degrees <- posterior$est.degrees
est.eta <- posterior$est.eta
est.latent.pos <- posterior$est.latent.pos
est.gi <- getGi(est.degrees, est.eta)


hist(est.gi) ## these seems similar to hist(est_greg)
## need to check it some more to be sure...
