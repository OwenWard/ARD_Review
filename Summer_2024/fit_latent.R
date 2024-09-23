### first initial example, which does manage to fit the 
### 2015 latent space model for simulated data



source("Summer_2024/latent_surface_model.R")


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

## need to figure out these estimated latent positions,
## I guess these are in polar coords or something



## compare the true to mean degrees here

true_degrees <- apply(y_sim, 1, sum)

deg_hat <- apply(est.degrees, 2, mean)
deg_lower <- apply(est.degrees, 2, quantile, p = .025)
deg_upper <- apply(est.degrees, 2, quantile, p = .975)


### need to make this a bit tidier

plot(true_degrees, 
     deg_hat)

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

## I think this is kind of working, just need to figure out the relationship 
## here



## we don't need this because this is inferring something for nodes 
## which are not in the ard sample, I think
## not sure below here will work, but anyway

posteriorAll <- getPosteriorAllnodes(distance.matrix, est.gi, est.latent.pos, Knn.K, ls.dim)
est.gi.all <- posteriorAll$est.gi.all
est.latent.pos.all <- posteriorAll$est.latent.pos.all
g.sims <- simulate.graph.all(est.degrees, est.eta, est.latent.pos, est.gi, est.gi.all, est.latent.pos.all, ls.dim)



generateRandomInitial <- function(n, p) {
  z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  z <- sweep(z, MARGIN = 1, 1 / sqrt(rowSums(z^2)), `*`)
  return(z)
}

getPosteriorAllnodes <- function(distance.matrix, est.gi, est.latent.pos, Knn.K, ls.dim) {
  n.ARD <- dim(distance.matrix)[2]
  n.nonARD <- dim(distance.matrix)[1]
  est.gi.all <- NULL
  est.latent.pos.all <- NULL
  for (ind in 1:dim(est.gi)[1]) {
    g.ARD <- est.gi[ind, ]
    z.ARD <- matrix(est.latent.pos[ind, ], byrow = F, nrow = n.ARD, ncol = ls.dim)

    g.nonARD <- NULL
    z.nonARD <- NULL
    for (i in 1:n.nonARD) {
      if (sort(distance.matrix[i, ])[1] != 0) {
        K.nn <- order(distance.matrix[i, ])[1:Knn.K]
        weights <- (1 / sort(distance.matrix[i, ])[1:Knn.K]) / sum(1 / sort(distance.matrix[i, ])[1:Knn.K])
        g.nonARD <- c(g.nonARD, sum(g.ARD[K.nn] * weights))
        z.tmp <- colSums(sweep(z.ARD[K.nn, ], MARGIN = 1, weights, "*"))
        z.nonARD <- rbind(z.nonARD, z.tmp / sqrt(sum(z.tmp^2)))
      } else {
        g.nonARD <- c(g.nonARD, g.ARD[order(distance.matrix[i, ])[1]])
        z.nonARD <- rbind(z.nonARD, z.ARD[order(distance.matrix[i, ])[1]])
      }
    }
    g <- c(g.ARD, g.nonARD)
    est.gi.all <- rbind(est.gi.all, g)
    z <- rbind(z.ARD, z.nonARD)
    est.latent.pos.all <- rbind(est.latent.pos.all, c(z))
  }
  return(list(est.gi.all = est.gi.all, est.latent.pos.all = est.latent.pos.all))
}

getPosterior <- function(out, n.iter, m.iter, n.thin, n) {
  est.degrees <- NULL
  est.eta <- NULL
  est.latent.pos <- NULL
  for (n.ind in (n.iter / n.thin / 2 + 1):(n.iter / n.thin)) {
    for (m.ind in 1:m.iter) {
      est.degrees <- rbind(est.degrees, out$sims[n.ind, m.ind, 1:n])
      est.eta <- c(est.eta, out$sims[n.ind, m.ind, ][length(out$sims[n.ind, m.ind, ])])
      est.latent.pos <- rbind(est.latent.pos, c(out$sims.latent[n.ind, m.ind, 1:n, ]))
    }
  }
  return(list(est.degrees = est.degrees, est.eta = est.eta, est.latent.pos = est.latent.pos))
}

getGi <- function(est.degrees, est.eta) {
  gi.m <- NULL
  for (ind in 1:length(est.eta)) {
    nexp.gi <- sqrt(sum(exp(est.degrees[ind, ])) * cp.fcn(est.eta[ind]) / cp.fcn(.000001))
    gi <- exp(est.degrees[ind, ]) / nexp.gi * cp.fcn(est.eta[ind]) / cp.fcn(.000001)
    gi.m <- rbind(gi.m, gi)
  }
  return(log(gi.m))
}
