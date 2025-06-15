# Simulate data from the NSUM package's
#   subpopulation sizes and parameters for the
#   McCarty data (Killworth et al., 1998)


# Load packages ----------------------------------------------------------------
library(tidyverse)
library(rstan)
library(tidymodels)
library(posterior)
library(bayesplot)
options(mc.cores = parallel::detectCores())

# load prereq's for NSUM (see NSUM/DESCRIPTION file)
library(MASS)
library(MCMCpack)

source(paste0(getwd(),"/helper/helper_model_checking.R"))
source(paste0(getwd(),"/helper/helper_plots.R"))

gg_color_hue <- function(n) { 
  hues = seq(15, 375, length = n + 1) 
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Load NSUM files from .zip  ---------------------------------------------------
# accompanies the Maltiel, Raftery, McCormick (2013) paper
nsum_pkg_path <- paste0(getwd(),"/NSUM/")
nsum_files <- list.files(paste0(nsum_pkg_path,"R/"))

for (file_i in nsum_files) source(paste0(nsum_pkg_path,"R/",file_i))

# Simulate McCarty data  --------------------------------------

## Import McCarty parameters  -----------------------------------------------------------
load(paste0(nsum_pkg_path,"/data/McCarty.Rdata"))

McCarty
# $known: counts for each of 29 known subpopulations
# $unknown: (500,000) (estimated?) size of unknown population(s)
## K = total # of known/unknown populations = 30
# $N: (250,000,000) the (known) total population size
#       ASIDE: Sim. data #1 uses N = 100,000
# $mu: (5.36) location parameter for the log-normal distribution of network 
#               degrees
# $sigma: (0.8) scale parameter for the log-normal distribution of net. degrees
# $rho: dispersion parameters / barrier effects (0<=rho[k]<=1) for each
#         of the 30 total known/unknown subpop.s
# $tauK: (1) multipliers for transmission biases (length = # of unknown
#             subpop.s )

# McCarty model:  d_i ~ log-Normal( mu= 5.36, sigma= 0.8 )
#   based on mean(sim.bar$d) ~= 253, E[d_i] = exp(5.36) ~= 212.7
# Owen's latent ARD model: d_i ~ Normal( mu= -12, sigma= 1 )
# => scale is very different

# other differences:
# n_subpop = 30 (vs. 15)
# perc_subpop is defined by known/N (vs. sim'd from a gamma)

# compare subpop. proportions:
#set.seed(100)
#perc_subpop <- round(rgamma(n = 15, shape = 1, rate = 10), 3)
#summary(sort(perc_subpop))

#summary(sort(with(McCarty, known/N)))
# McCarty subpop prop.s are much smaller

## Simulate from barrier effects model --------------------------------------
nindiv_sim <- 1000
sim.bar <- with(McCarty, nsum.simulate(nindiv_sim, known, unknown, N, 
                                       model="barrier", 
                                       mu, sigma, rho))
# I've verified that this simulation code looks right

# simulates nindiv_sim x K matrix of ARData
#           nindiv_sim personal degree sizes (i.e., est. size of ego-network)


## True degrees ----------------------------------------------------------------

deg_hist <- ggplot(data.frame(x=sim.bar$d),aes(x)) + 
  geom_histogram() + 
  theme_single() +
  ylab("") + xlab("Sample Degree")

deg_hist

# MCMC and Stan setup ----------------------------------------------------------

iter_warmup <- 3000
iter_warmup_nsum <- 34500 # Maltiel et al. use 30k iter, 5k burnin
iter_after_warmup <- 500 # keep this small-ish because we will use all of these
#  iterations in our ppcs
# note:  in cmdstan, total # of iter. = iter_sampling + iter_warmup
#        in rstan, total # of iter = iter ( and warmup < iter)

# MCMC from NSUM ---------------------------------------------------------------
## estimate unknown population size
dat.bar <- sim.bar$y
start_mcmc <- Sys.time()
mcmc <- with(McCarty, nsum.mcmc( dat.bar, known, N, model="barrier",
                                 indices.k = ncol(dat.bar),
                                 iterations = iter_after_warmup, 
                                 burnin = iter_warmup_nsum )) # total iter = iterations + burnin
end_mcmc <- Sys.time()
print(end_mcmc - start_mcmc) # 35 k total iter => (~1.5min desktop; ~17min laptop)
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

## view posterior distribution of barrier effect parameters for the first 
##     subpopulation
hist(mcmc$rho.values[1,])
abline(v=McCarty$rho,col="grey40")
# should be shifted left a bit

## Subpop. sizes ---------------------------------------------------------------
# want to show that we can do as well (e.g., N_k est.s, pppc) as the Stan models
#   which also have entry-level parameters

size_true <- data.frame( parameter = c(names(McCarty$known),"unknown"),
                         true = c(McCarty$known,McCarty$unknown) )

# doesn't appear to estimate N_k for known subpop.s???
size_ests_plotdata_nsum <- mcmc_intervals_data( matrix( 
  mcmc$NK.values,
  ncol = 1,
  dimnames = list(NULL,
#                  c(names(McCarty$known),"unknown")) )*
#    McCarty$N ) %>%
                  "unknown"))) 

size_ests_plotdata_nsum %>%
  left_join(size_true) %>% 
  pivot_longer( cols=c(m,true),
                names_to = "point_type", 
                values_to = "point_value") %>%
  ggplot() +
  geom_segment( aes(x=ll,xend=hh,y=parameter,yend=parameter),
                #color = get_color("mid"),
                linewidth = .5 ) +
  geom_segment( aes(x=l,xend=h,y=parameter,yend=parameter),
                linewidth = 2 ) +
  geom_point( aes( x = point_value, y = parameter,
                   col = point_type), size = 1.5 ) +
  bayesplot_theme_get() +
  labs(y = "Subpopulation", x = "Size" )

## PPC: P(Y_ik = j) ------------------------------------------------------------

### PPDraws --------------------------------------------------------------------

start_nsum_ppds <- Sys.time()
ysim_nsum <- NULL
for (i in 1:iter_after_warmup){
  ysim_nsum_i <- nsum.simulate( nrow(sim.bar$y), 
                                McCarty$known,
                                mcmc$NK.values[i],
                                McCarty$N,
                                model="barrier",
                                mcmc$mu.values[i],
                                mcmc$sigma.values[i],
                                mcmc$rho.values[,i] )$y
  colnames(ysim_nsum_i) <- 1:ncol(sim.bar$y)
  rownames(ysim_nsum_i) <- 1:nrow(sim.bar$y)
  
  ysim_nsum_i_long <- ysim_nsum_i %>% 
    as_tibble() %>% 
    mutate(node = row_number()) %>% 
    pivot_longer( !node, 
                  values_to = "y", 
                  names_to = "subpop" ) %>%
    mutate( par = paste0("y_sim[",node,",",subpop,"]"),
            .draw = i) %>%
    dplyr::select(par, y, .draw)
  
  ysim_nsum <- ysim_nsum %>% 
    bind_rows(ysim_nsum_i_long)
}
end_nsum_ppds <- Sys.time()
end_nsum_ppds - start_nsum_ppds

ysim_nsum_dfdraws <- ysim_nsum %>%
  pivot_wider( names_from = par,
               values_from = y ) %>%
  posterior::as_draws_df()

### PPChecks -------------------------------------------------------------------

ppc_nsum <- construct_ppc( ysim_nsum_dfdraws, sim.bar$y)
ppc_nsum_plot <- with(ppc_nsum, plot_ests_all( ppc_draws,
                                                   y_tibble ))
ppc_nsum_plot$final_plot

# Stan models ------------------------------------------------------------------

stan_data_null <- list( N = nrow(sim.bar$y),
                        K = ncol(sim.bar$y),
                        y = sim.bar$y,
                        n_known = length(McCarty$known),
                        idx = 1:length(McCarty$known),
                        known_prev = sum(McCarty$known)/McCarty$N,
                        N_totalpop = McCarty$N )

## Erdos Renyi -----------------------------------------------------------------
# (assumes a fixed/common degree)

stan_file_null_01 <- paste0(getwd(),"/stan_models/null_model_01_scaled.stan")
mod_null_01 <- stan_model(stan_file_null_01)

stan_fit_null_01 <- sampling( mod_null_01,
                      data = stan_data_null,
                      seed = 123,
                      chains = 4,
                      #iter_sampling = 1000, #iter_warmup = 1000,
                      iter = iter_warmup + iter_after_warmup,
                      warmup = iter_warmup, 
                      #parallel_chains = 4,
                      refresh = 100 ) # how often progress is reported

summary( summary(stan_fit_null_01)$summary[,"Rhat"] )
rstan::traceplot(stan_fit_null_01, inc_warmup = FALSE, pars="scaled_log_d" )

## check the population sizes 
plot( stan_fit_null_01,
      pars = "scaled_log_d" )

### Subpop. sizes --------------------------------------------------------------

size_ests_plotdata_ER <- mcmc_intervals_data( matrix( 
  as.array(rstan::extract(stan_fit_null_01,pars="b")[[1]]),
  ncol = stan_data_null$K,
  dimnames = list(NULL,
                  c(names(McCarty$known),"unknown")) )*
  McCarty$N )

nfacets <- 4
facet_quantile_lims <- seq(0,1,length=nfacets+1)
facet_labs <- paste0( levels(cut(
  seq(0,1,length=nfacets+1)*100, 4,
  include.lowest = TRUE)),
  "th Percentiles")

size_ests_plotdata_ER  %>%
  arrange(m) %>% # sort subpop.s by posterior median
  mutate( parameter_sorted = factor(parameter,
                                    levels = size_ests_plotdata_ER %>% 
                                      arrange(m) %>% 
                                      pull(parameter)) ) %>%
  mutate( facet_factor = cut(m,quantile(m,probs=facet_quantile_lims),
                             include.lowest=TRUE,
                             labels = facet_labs) ) %>%
  left_join(size_true) %>% 
  pivot_longer( cols=c(m,true),
                names_to = "point_type", 
                values_to = "point_value") %>%
  ggplot() +
  geom_segment( aes(x=ll,xend=hh,y=parameter_sorted,yend=parameter_sorted),
                #color = get_color("mid"),
                linewidth = .5 ) +
  facet_wrap(vars(facet_factor),
             scales="free") +
  geom_segment( aes(x=l,xend=h,y=parameter,yend=parameter),
                linewidth = 2 ) +
  #geom_point(aes(x=m,y=parameter),size=2,shape=21,
  #           col="lightblue") +
  #geom_point(aes(x=true,y=parameter),
  #           size = 2, shape=2, col="green") +
  geom_point( aes( x = point_value, y = parameter,
                   col = point_type), size = 1.5 ) +
  #scale_y_discrete()
  bayesplot_theme_get() +
  labs(y = "Subpopulation", x = "Size" )

### PPC: P(Y_ik = j) ------------------------------------------------------------

ppc_null_1 <- construct_ppc(stan_fit_null_01, sim.bar$y)
ppc_null_1_plot <- with(ppc_null_1, plot_ests_all( ppc_draws,
                                               y_tibble ))
ppc_null_1_plot$final_plot

## Varying degree --------------------------------------------------------------

stan_file_null_02 <- paste0(getwd(),"/stan_models/null_model_02_scaled.stan")
mod_null_02 <- stan_model(stan_file_null_02)

stan_fit_null_02 <- sampling( mod_null_02,
                              data = stan_data_null,
                              seed = 123,
                              chains = 4,
                              #iter_sampling = 1000, #iter_warmup = 1000,
                              iter = iter_warmup + iter_after_warmup,
                              warmup = iter_warmup, 
                              #parallel_chains = 4,
                              refresh = 100 ) # how often progress is reported

summary( summary(stan_fit_null_02)$summary[,"Rhat"] )

rstan::traceplot(stan_fit_null_02, inc_warmup = FALSE,
          pars=c("scaled_beta"))

rstan::traceplot(stan_fit_null_02, inc_warmup = FALSE,
          pars=c("scaled_log_d[1]",
                 "scaled_log_d[2]"))

### Subpop. sizes --------------------------------------------------------------

size_ests_plotdata_degree <- mcmc_intervals_data( matrix( 
  as.array(rstan::extract(stan_fit_null_02,pars="b")[[1]]),
  ncol = stan_data_null$K,
  dimnames = list(NULL,
                  c(names(McCarty$known),"unknown")) )*
    McCarty$N )

nfacets <- 4
facet_quantile_lims <- seq(0,1,length=nfacets+1)
facet_labs <- paste0( levels(cut(
  seq(0,1,length=nfacets+1)*100, 4,
  include.lowest = TRUE)),
  "th Percentiles")

size_ests_plotdata_degree  %>%
  arrange(m) %>% # sort subpop.s by posterior median
  mutate( parameter_sorted = factor(parameter,
                                    levels = size_ests_plotdata_degree %>% 
                                      arrange(m) %>% 
                                      pull(parameter)) ) %>%
  mutate( facet_factor = cut(m,quantile(m,probs=facet_quantile_lims),
                             include.lowest=TRUE,
                             labels = facet_labs) ) %>%
  left_join(size_true) %>% 
  pivot_longer( cols=c(m,true),
                names_to = "point_type", 
                values_to = "point_value") %>%
  ggplot() +
  geom_segment( aes(x=ll,xend=hh,y=parameter_sorted,yend=parameter_sorted),
                #color = get_color("mid"),
                linewidth = .5 ) +
  facet_wrap(vars(facet_factor),
             scales="free") +
  geom_segment( aes(x=l,xend=h,y=parameter,yend=parameter),
                linewidth = 2 ) +
  #geom_point(aes(x=m,y=parameter),size=2,shape=21,
  #           col="lightblue") +
  #geom_point(aes(x=true,y=parameter),
  #           size = 2, shape=2, col="green") +
  geom_point( aes( x = point_value, y = parameter,
                   col = point_type), size = 1.5 ) +
  #scale_y_discrete()
  bayesplot_theme_get() +
  labs(y = "Subpopulation", x = "Size" )

### PPC: P(Y_ik = j) -----------------------------------------------------------

ppc_null_2 <- construct_ppc(stan_fit_null_02, sim.bar$y)
ppc_null_2_plot <- with(ppc_null_2, plot_ests_all( ppc_draws,
                                                   y_tibble ))
ppc_null_2_plot$final_plot

## Overdispersed model (Zheng et al., 2006) ------------------------------------

stan_file_zheng <- paste0(getwd(),"/stan_models/zheng_et_al_2006_scaled.stan")
mod_zheng <- stan_model(stan_file_zheng)

stan_fit_zheng <- sampling( mod_zheng,
                            data = stan_data_null,
                            seed = 123,
                            chains = 4,
                            #iter_sampling = 1000, #iter_warmup = 1000,
                            iter = iter_warmup + iter_after_warmup,
                            warmup = iter_warmup, 
                            #parallel_chains = 4,
                            refresh = 100 ) # how often progress is reported

summary( summary(stan_fit_zheng)$summary[,"Rhat"] )

rstan::traceplot(stan_fit_zheng, inc_warmup = FALSE,
                 pars=c("sigma_alpha"))

rstan::traceplot(stan_fit_zheng, inc_warmup = FALSE,
                 pars=c("scaled_beta"))

### Subpop. sizes --------------------------------------------------------------

size_ests_plotdata_zheng <- mcmc_intervals_data( matrix( 
  exp( as.array(rstan::extract(stan_fit_zheng,pars="scaled_beta")[[1]]) ),
  ncol = stan_data_null$K,
  dimnames = list(NULL,
                  c(names(McCarty$known),"unknown")) )*
    McCarty$N )

nfacets <- 4
facet_quantile_lims <- seq(0,1,length=nfacets+1)
facet_labs <- paste0( levels(cut(
  seq(0,1,length=nfacets+1)*100, 4,
  include.lowest = TRUE)),
  "th Percentiles")

size_ests_plotdata_zheng  %>%
  arrange(m) %>% # sort subpop.s by posterior median
  mutate( parameter_sorted = factor(parameter,
                                    levels = size_ests_plotdata_zheng %>% 
                                      arrange(m) %>% 
                                      pull(parameter)) ) %>%
  mutate( facet_factor = cut(m,quantile(m,probs=facet_quantile_lims),
                             include.lowest=TRUE,
                             labels = facet_labs) ) %>%
  left_join(size_true) %>% 
  pivot_longer( cols=c(m,true),
                names_to = "point_type", 
                values_to = "point_value") %>%
  ggplot() +
  geom_segment( aes(x=ll,xend=hh,y=parameter_sorted,yend=parameter_sorted),
                #color = get_color("mid"),
                linewidth = .5 ) +
  facet_wrap(vars(facet_factor),
             scales="free") +
  geom_segment( aes(x=l,xend=h,y=parameter,yend=parameter),
                linewidth = 2 ) +
  #geom_point(aes(x=m,y=parameter),size=2,shape=21,
  #           col="lightblue") +
  #geom_point(aes(x=true,y=parameter),
  #           size = 2, shape=2, col="green") +
  geom_point( aes( x = point_value, y = parameter,
                   col = point_type), size = 1.5 ) +
  #scale_y_discrete()
  bayesplot_theme_get() +
  labs(y = "Subpopulation", x = "Size" )


zheng_subpops_plotdata <- size_ests_plotdata_zheng  %>%
  arrange(m) %>% # sort subpop.s by posterior median
  mutate( parameter_sorted = factor(parameter,
                                    levels = size_ests_plotdata_zheng %>% 
                                      arrange(m) %>% 
                                      pull(parameter)) ) %>%
  mutate( facet_factor = cut(m,quantile(m,probs=facet_quantile_lims),
                             include.lowest=TRUE,
                             labels = facet_labs) ) %>%
  left_join(size_true) %>% 
  pivot_longer( cols=c(m,true),
                names_to = "point_type", 
                values_to = "point_value") %>%
  mutate( Subpopulation = ifelse( parameter == "unknown",
                                  "Unknown",
                                  "Known" ))

zheng_subpops_plot <- zheng_subpops_plotdata %>%
  ggplot() +
  geom_segment( aes(x=ll,xend=hh,y=parameter_sorted,yend=parameter_sorted,col=Subpopulation),
                linewidth = .5 ) +
  facet_wrap(vars(facet_factor),
             scales="free") +
  geom_segment( aes(x=l,xend=h,y=parameter,yend=parameter,col=Subpopulation),
                linewidth = 2 ) +
  geom_point( data = zheng_subpops_plotdata %>%
                filter( point_type == "m" ),
              aes( x = point_value, y = parameter,col=Subpopulation), size = 3, shape = "|" ) +
  geom_point( data = zheng_subpops_plotdata %>%
                filter( point_type != "m" ),
              aes( x = point_value, y = parameter)) + 
  bayesplot_theme_get() +
  theme(legend.position = "bottom") + 
  labs(y = "Subpopulation", x = "Size" )


### PPC: P(Y_ik = j) -----------------------------------------------------------

ppc_zheng <- construct_ppc(stan_fit_zheng, sim.bar$y)
ppc_zheng_plot <- with(ppc_zheng, plot_ests_all( ppc_draws,
                                                   y_tibble ))
ppc_zheng_plot$final_plot

# Stansum: Stan version of NSUM ------------------------------------------------
#   GitHub:  https://github.com/coalesce-lab/stansum/blob/master/src/stan/MaltielBEM_count.stan
#   Reference:  https://coalesce-lab.github.io/stansum/reference/maltiel_bem.html
#  requires cmdstanr package, so instead just take the .stan file and implement
#   as above

tryBem <- FALSE
# this stan code treats ALL N_k / N as observed!!

if (tryBem){

  stan_data_bem <- c(stan_data_null[c("N","K","y")],
                   list( 
                     # fractional sizes of sub-populations
                     m = c(McCarty$known,McCarty$unknown)/
                           stan_data_null$N_totalpop,
                     # min. degree
                     L = rep(25,stan_data_null$N) ))

  nsum_inits <- with(stan_data_null, 
                   killworth.start( y, McCarty$known, N_totalpop ))

  init_bem_killworth <- function(...) list( d_raw = nsum_inits$d.start,
                                          mu = nsum_inits$mu.start,
                                          sigma = nsum_inits$sigma.start )

  stan_file_maltiel <- paste0(getwd(),"/stan_models/MaltielBEM_count.stan")
  mod_bem <- stan_model(stan_file_maltiel)

  stan_fit_bem <- sampling( mod_bem,
                            data = stan_data_bem,
                            seed = 123,
                            chains = 4,
                            init = init_bem_killworth,
                            #iter_sampling = 1000, #iter_warmup = 1000,
                            iter = iter_after_warmup,
                            warmup = 3000, 
                            #parallel_chains = 4,
                            refresh = 100 ) # how often progress is reported

  summary( summary(stan_fit_bem)$summary[,"Rhat"] )
}

# All models -------------------------------------------------------------------

modelnames <- rev(c("Erdos Renyi",
                    "Varying Degree",
                    "Overdispersed",
                    "Barrier Effects\n(NSUM)"))

## Degree distributions --------------------------------------------------------

degs_all <- data.frame( #true = sim.bar$d,
                        nsum = rowMeans(mcmc$d.values),
                        er = rep( exp(summary(stan_fit_null_01)$summary[
                          "scaled_log_d","mean"]), nindiv_sim),
                        deg = exp(summary(stan_fit_null_02)$summary[
                          paste0("scaled_log_d[",1:nindiv_sim,"]"),"mean"]),
                        overd = exp(summary(stan_fit_zheng)$summary[
                          paste0("scaled_alpha[",1:nindiv_sim,"]"),"mean"]) )

deg_hist_all <- degs_all %>% 
  pivot_longer( everything(), values_to = "degree",
                names_to = "model_abbr" ) %>%
  mutate( model = factor( model_abbr,
                          levels = c("nsum","overd","deg","er"),
                          labels = modelnames )) %>%
  filter( model_abbr != "er" ) %>%
  ggplot(aes(degree,col=model)) + 
  geom_histogram( data = data.frame(degree=sim.bar$d),
                  aes(degree,y=after_stat(density)),
                  col="lightgrey",alpha=.1,bins=50) +
  geom_density( data = data.frame(degree=sim.bar$d,model="Observed"),
                aes(degree),
                lty = 2, key_glyph = "path" ) +
  geom_segment(aes(x = mean(degs_all$er), 
                   xend = mean(degs_all$er), 
                   y=0, yend=Inf,
                   col="Erdos-Renyi")) + 
  geom_density(key_glyph = "path") + 
  labs(x="Sample Degree", y="") + 
  theme_bw() + 
  scale_color_manual( values=c(gg_color_hue(4),"black")) + 
  guides(color = guide_legend(title = "Model",
                              override.aes = 
                                list(linetype=c(1,1,1,1,2),
                                     col=c(gg_color_hue(4),"black")) ) ) +
  theme(legend.position = "bottom")

deg_hist_all

## Subpop. sizes ---------------------------------------------------------------
size_ests_plotdata_all <- size_ests_plotdata_nsum %>% 
  mutate( model = "Barrier Effects\n(NSUM)") %>%
  bind_rows( size_ests_plotdata_ER %>% mutate( model = "Erdos Renyi" ),
             size_ests_plotdata_degree %>% mutate( model = "Varying Degree"),
             size_ests_plotdata_zheng %>% mutate( model = "Overdispersed") ) %>%
  mutate( model = factor( model, levels = modelnames ) )

size_ests_all_plot <- size_ests_plotdata_all %>% 
  filter(parameter=="unknown") %>%
  ggplot() +
  geom_segment( aes(x=ll,xend=hh,y=model,yend=model),
                linewidth = .5 ) +
  geom_segment( aes(x=l,xend=h,y=model,yend=model),
                linewidth = 2 ) +
  geom_point( aes( x = m, y = model,
                   col=model), size = 3.5 ) +
  geom_vline(xintercept = size_true %>% 
               filter(parameter=="unknown") %>% 
               pull(true), lty=2) +
  #scale_y_discrete()
  bayesplot_theme_get() + theme(legend.position = "bottom") +
  guides(col=guide_legend(title="Model")) +
  labs(y = "Model", x = "Size" )

size_ests_all_plot

## PPC: P(Y_ik = j) ---------------------------------------------------------------

ppc_all_plotdata <- ppc_nsum_plot$final_plot_data %>%
  mutate( model = "Barrier Effects (NSUM)") %>%
  bind_rows( ppc_null_1_plot$final_plot_data %>% mutate( model = "Erdos Renyi" ),
             ppc_null_2_plot$final_plot_data %>% mutate( model = "Varying Degree" ),
             ppc_zheng_plot$final_plot_data %>% mutate( model = "Overdispersed" )) %>%
  mutate( model = factor( model, levels = rev(modelnames) ) )

ppc_all_plotdata_lims <- ppc_all_plotdata %>%
  group_by(model,count) %>%
  summarize( square_lim = max(true_prop,avg),
             ebar_ht = .075*square_lim )

ppc_all_plot <- ppc_all_plotdata %>%
  filter( count %in% paste0("P(y_ik =",c(0,1,3,5,10),")") ) %>%
  mutate( covers_true = (true_prop >= lower) & (true_prop <= upper)) %>%
  left_join(ppc_all_plotdata_lims) %>%
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

## Save plots ------------------------------------------------------------------

ggsave("Figures/McCartyData_Degrees_061225.pdf",
       deg_hist,
       width = 5, height = 5)

ggsave("Figures/McCartyData_Degrees_Est_061225.pdf",
       deg_hist_all,
       width = 5, height = 5)

ggsave("Figures/McCartyData_Subpop_Known_Zheng_061225.pdf",
       zheng_subpops_plot,
       width = 10, height = 10)


ggsave("Figures/McCartyData_Subpop_Ests_061225.pdf",
       size_ests_all_plot,
       width = 8, height = 6)

ggsave("Figures/McCartyData_PPCs_061225.pdf",
       ppc_all_plot,
       width = 12, height = 10)
