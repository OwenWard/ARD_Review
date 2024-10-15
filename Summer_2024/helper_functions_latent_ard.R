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
