# see DESCRIPTION file 
library(MASS)
library(MCMCpack)

nsum_pkg_path <- '/Users/annasmith/Library/CloudStorage/OneDrive-UniversityofKentucky/Research Projects/ARD Social Networks/NSUM/'
nsum_files <- list.files(paste0(nsum_pkg_path,"R/"))

for (file_i in nsum_files) source(paste0(nsum_pkg_path,"R/",file_i))

## load data
#data(McCarty)
load(paste0(nsum_pkg_path,"/data/McCarty.Rdata"))

McCarty
# $known: counts for each of 29 known subpopulations
# $unknown: (500,000) (estimated?) size of unknown population(s)
## K = total # of known/unknown populations = 30
# $N: (250,000,000) the (known) total population size
# $mu: (5.36) location parameter for the log-normal distribution of network degrees
# $sigma (0.8) scale parameter for the log-normal distribution of network degrees
# $rho dispersion parameters for the barrier effects (0<=rho[k]<=1) for each of the 30 total known/unknown subpop.s
# $tauK (1) multipliers for transmission biases (length = # of unknown subpop.s )

#McCarty_matchOwen <- McCarty
#McCarty_matchOwen$N <- 1e5

# McCarty model:  d_i ~ log-Normal( mu= 5.36, sigma= 0.8 )
#   based on mean(sim.bar$d) ~= 253, E[d_i] = exp(5.36) ~= 212.7
# Owen's latent ARD model: d_i ~ Normal( mu= -12, sigma= 1 )
# => scale is very different

# other differences:
# n_subpop = 30
# perc_subpop is defined by known/N

#set.seed(100)
#perc_subpop <- round(rgamma(n = 15, shape = 1, rate = 10), 3)
#summary(sort(perc_subpop))

#summary(sort(with(McCarty, known/N)))
# McCarty subpop prop.s are much smaller


## simulate from model with barrier effects
#nindiv_sim <- 100  # from nsum example
nindiv_sim <- 1000
sim.bar <- with(McCarty, nsum.simulate(nindiv_sim, known, unknown, N, model="barrier", 
                                       mu, sigma, rho))
# I've verified that this simulation code looks right

# simulates nindiv_sim x K matrix of ARData
#           nindiv_sim personal degree sizes (/= rowSums(y))

## estimate unknown population size
dat.bar <- sim.bar$y
start_mcmc <- Sys.time()
mcmc <- with(McCarty, nsum.mcmc( dat.bar, known, N, model="barrier",
                                 iterations=30000, # per Maltiel et al.
                                 burnin=5000 )) # total iter = iterations + burnin
end_mcmc <- Sys.time()
print(end_mcmc - start_mcmc) # 35 k total iter => ~ 1.5min
# when they fit this model they use floor(d_i) in each iteration

## view posterior distribution of subpopulation sizes for the first subpopulation
hist(mcmc$NK.values[1,])
abline(v=McCarty$unknown,col="grey40")
# overestimates 

## posterior distr. for d[i]'s
pick_is <- sample.int(nindiv_sim,5)
par(mfrow=c(1,5))
for (i in pick_is){
  hist( mcmc$d.values[i,] )
  abline( v = sim.bar$d[i], col="grey40" )
}
# typically underestimates?

## posterior distr. of mu
hist( mcmc$mu.values )
abline( v = McCarty$mu, col="grey40" )
# looks great!

## posterior distr. of sigma
hist( mcmc$sigma.values )
abline( v = McCarty$sigma, col="grey40" )
# overestimates

## view posterior distribution of barrier effect parameters for the first subpopulation
hist(mcmc$rho.values[1,])
abline(v=McCarty$rho,col="grey40")
# should be shifted left a bit


## Fit the stan model

# check stan code

stan_data_null <- list( N = nrow(sim.bar$y),
                        K = ncol(sim.bar$y),
                        y = sim.bar$y,
                        n_known = length(McCarty$known),
                        idx = 1:length(McCarty$known),
                        known_prev = sum(McCarty$known)/McCarty$N )

stan_file_null_01 <- here("stan_models")


