#### Sept 23 2024 ####

### working on simulating and fitting existing models for this data
### aim is to simulate and fit models from 2006, 2015 and maybe 2019 paper




# Some Latent Space Simulation Code 2009 --------------------------------------------------

## code taken from Tian's Dropbox, 
## Project_LatentMix/prgm/simulations/fake_generation.R

dist_xy <- function(x, y, method_used = "euclidean") {
  n_x <- nrow(x)
  n_y <- nrow(y)
  dist_mat <- as.matrix(dist(rbind(x, y), method = method_used))
  dist_mat <- dist_mat[1:n_x, n_x + 1:n_y]
  return(dist_mat)
}

# center (0,0)
n_popu <- 1000000
perc_subpopu <- 0.001
n_subpopu <- round(n_popu * perc_subpopu)
n_sample <- 1000
sub_k <- 5 * 3
dim_p <- 2
rho_sub <- rep(0, sub_k)
center_sub <- matrix(runif(sub_k * 2, 0, 1), sub_k, 2)
var_sub <- 1.7^(1:sub_k + 1)
lamda <- 100
degree_mean <- 500
degree_sd <- 50

# population are distributed uniformly in the unit square (cube if higher dimension)
x_popu <- matrix(NA, n_popu, 2)
for (i in 1:dim_p) {
  x_popu[, i] <- runif(n_popu, 0, 1)
}

x_sample_degree <- rnorm(n_sample, degree_mean, degree_sd)
x_sample <- x_popu[sample(1:n_popu, n_sample), ]

fake_data <- NULL

par(mfrow = c(sub_k / 3, 3 * 2), pty = "s", font_main = 1, cex_main = 1)

for (i in 1:sub_k) {
  print(c(i, var_sub[i], date()))
  
  # subpopulation structure
  cov_sup <- matrix(rep(rho_sub[i], dim_p^2), dim_p, dim_p)
  diag(cov_sup) <- 1
  cov_sup <- cov_sup * var_sub[i]
  eigen_cov_sup <- eigen(solve(cov_sup))
  
  temp_dist <- x_popu
  temp_center <- center_sub[i, ]
  for (j in 1:dim_p) {
    temp_dist[, j] <- temp_dist[, j] - temp_center[j]
  }
  temp_dist <- (temp_dist) %*% t(eigen_cov_sup$vectors) %*% 
    diag(sqrt(eigen_cov_sup$values))
  temp_dist <- sqrt(rowSums(temp_dist^2))
  
  prob_sub <- 2 * (exp(-lamda * temp_dist) / (1 + exp(-lamda * temp_dist)))
  prob_sub <- prob_sub * n_popu * perc_subpopu / sum(prob_sub)
  prob_sub <- ifelse(prob_sub > 1, 1, prob_sub)
  sub_member <- rbinom(n_popu, 1, p = prob_sub)
  print(sum(sub_member))
  plot(x_popu[sub_member == 1, ],
       xlim = c(0, 1), ylim = c(0, 1), pch = ".", 
       cex = 1.2, col = "light blue", xlab = "", ylab = "",
       main = paste("Subpopulation", format(i))
  )
  points(center_sub[i, 1], center_sub[i, 2], pch = 16, col = "blue", cex = 1.3)
  points(x_sample, pch = ".", col = "light salmon", cex = 1.2)
  
  x_sub <- x_popu[sub_member == 1, ]
  
  temp_dist_sub <- x_sub
  temp_dist_sample <- x_sample
  
  dist_sub_sample <- dist_xy(temp_dist_sample, temp_dist_sub)
  dist_sample_sample <- dist_xy(temp_dist_sample, temp_dist_sample)
  est_degree_sample <- rowSums(2 * (exp(-lamda * dist_sample_sample) / 
                                      (1 + exp(-lamda * dist_sample_sample)))) /
    n_sample * n_popu
  
  ## this seems to be used to assign nodes to the sample from
  ## that subpopulation, but not clear
  prob_sub_sample <- diag(x_sample_degree / est_degree_sample) %*% 
    (exp(-lamda * dist_sub_sample) / (1 + exp(-lamda * dist_sub_sample))) * 2
  
  prob_sub_sample <- ifelse(prob_sub_sample > 1, 1, prob_sub_sample)
  
  fake_data <- cbind(fake_data, rowSums(matrix(rbinom(n_sample * n_subpopu,
                                                      1,
                                                      prob_sub_sample),
                                               n_sample, n_subpopu)))
  
  hist(fake_data[, i], 
       main = paste("mean num known=",
                                    format(mean(fake_data[, i]))),
       breaks = seq(0, max(fake_data[, i]) + 1, 1) - 0.5)
}

dim(fake_data)

head(fake_data)

y_sim <- fake_data



# Simulate Data Corresponding to the Spherical Model ----------------------

library(rotasym)

n_pop <- 1000000
perc_subpopu <- 0.001
n_subpopu <- round(n_popu * perc_subpopu)
n_sample <- 1000
num_subpop <- 10


## first need to simulate the subpopulation locations on the sphere



## then assign nodes to subpopulations based on distance to subpopulations




# Fit the Latent ARD Model ------------------------------------------------



source("Summer_2024/latent_surface_model.R")
source("Summer_2024/fit_latent.R")

n <- 1000
K <- 15

dim(y_sim)


## need to figure out how to fit the latent surface model using all of this

ls.dim <- 3
n <- dim(y_sim)[1]
n.iter <- 3000
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
               # n.iter = n.iter,
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




# Compare results ---------------------------------------------------------


## think I should maybe be using x_sample_degree here, not sure

## compare the true to mean degrees here

## y is how many people they know in the entire population from that subpop
## but not clear what the degree estimates should then correspond to

true_degrees <- apply(y_sim, 1, sum)

deg_hat <- apply(est.degrees, 2, mean)
deg_lower <- apply(est.degrees, 2, quantile, p = .025)
deg_upper <- apply(est.degrees, 2, quantile, p = .975)




x_sample_degree * perc_subpopu

true_degrees[1:20]

as.numeric(deg_hat)[1:10] 


summary(deg_hat - true_degrees)


tibble(true = true_degrees,
       est = deg_hat) |> 
  filter(true == 10) |> 
  ggplot(aes(est)) +
  geom_histogram()

## some clear pattern of underestimation, which makes me think 
## I need to be scaling these somehow


### need to make this a bit tidier

plot(true_degrees, 
     deg_hat * 10)
abline(a = 0, b = 1, col = "red")


for (i in 1:n) lines(c(true_degrees[i], true_degrees[i]),
                     c(deg_lower[i], deg_upper[i]),
                     col = "grey")

points(true_degrees, deg_hat, pch = 16)



sorted_order <- order(true_degrees)

sorted_degree <- true_degrees[sorted_order]

sorted_hat <- deg_hat[sorted_order]
sorted_lower <- deg_lower[sorted_order]
sorted_upper <- deg_upper[sorted_order]

plot(sorted_degree, 
     sorted_hat)
for (i in 1:n) lines(c(sorted_deg[i], sorted_deg[i]),
                     c(sorted_lower[i], sorted_upper[i]),
                     col = "grey")

points(sorted_degree, sorted_hat, pch = 16)
