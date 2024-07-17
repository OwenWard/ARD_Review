# Plots a mixing matrix
plot_mixing <- function(M, y_max, ego_labels) {
  barplot(t(M[,1:4]), horiz = T, beside = T, xlim = c(-0.6, 0.6), ylim = c(0, y_max))
  barplot(t(-M[,4 + 1:4]), horiz = T, beside = T, add = T)
  for(i in 1:nrow(M)) {
    text(-0.5, 3 + 5*(i-1), ego_labels[i])
    text(0.5, 2.25 + 5*(i-1) - 1, "1-17")
    text(0.5, 2.25 + 5*(i-1) + .125, "18-24", col = "darkgray")
    text(0.5, 2.25 + 5*(i-1) + 1.125, "25-64", col = "gray")
    text(0.5, 2.25 + 5*(i-1) + 2.25, "65+", col = "lightgray")
    lines(c(-.42, -.42), 3 + 5*(i-1) + c(-1.5,1.5))
  }
  text(-.5, y_max, "Ego Groups")
  text(0.2, y_max, "Male Alters")
  text(-0.2, y_max, "Female Alters")
  text(.5, y_max, "Alter Ages")
}

# Simulates responses from parameters of a mixing stan model
simulate_mixing <- function(d, omega, M, beta, ego) {
  library(MASS)
  y <- matrix(nrow = length(d), ncol = ncol(beta))
  
  for(i in 1:nrow(y)) {
    for(j in 1:ncol(y)) {
      mu_ij <- omega[j] * d[i] * M[ego[i],] %*% beta[,j]
      y[i,j] <- rnegbin(1, mu_ij, omega[j])
    }
  }
  
  return(y)
}

# Plots a comparison between actual and simulated responses 
plot_comparison <- function(data1, data2, ego, ego_names) {
  for(e in 1:length(unique(ego))) {
    mu1 <- colMeans(data1[ego == e,])
    mu2 <- colMeans(data2[ego == e,])
    plot(mu1, mu2, pch = 21, bg = 'grey', main = ego_names[e], xlim = c(min(mu1)-0.5, max(mu1)+0.5), ylim = c(min(mu2)-0.2, max(mu2)+0.2), xlab = "Mean Responses (Actual)", ylab = "Mean Responses (Simulated)")
    abline(a = 0, b = 1, lty = 2)
    text(mu1, mu2, pos = ifelse(mu1 > mu2, 1, 3), colnames(data1))
  }
}