#
# ... Estimate the Degrees Using Responses to Names and Occupations ...
#
# ... 1. Run the Model in Stan ....

rm_names <- c(1,2,7,8)
rm_occ <- c(2)

mcmc_data_comb <- list(E = length(unique(ego_sex_age)), A = nrow(beta_names), J = ncol(beta_names[,-rm_names]), K = ncol(beta_occ_sex_age[,-rm_occ]), 
                       N = nrow(names_data), y = names_data[,-rm_names], z = occ_data[,-rm_occ], ego = ego_sex_age, BetaY = beta_names[,-rm_names], 
                       BetaZ = beta_occ_sex_age[,-rm_occ], theta_d = c(6.2,.5), theta_o = c(3,2), alpha = rep(1,nrow(beta_names))/nrow(beta_names))
fit_degree_comb <- sampling(degree_combined_mixing_fit, data = mcmc_data_comb, iter = 1000, chains = 4)

degree$CombSexAge <- round(colMeans(extract(fit_degree_comb)$d))
omega_comb_names_sex_age <- colMeans(extract(fit_degree_comb)$inv_omegaY)
omega_comb_occ_sex_age <- colMeans(extract(fit_degree_comb)$inv_omegaZ)

# ... 2. Plot Degree Distributions ....

par(mfrow=c(2,3))
plot_degrees(degree$CombSexAge, ego_sex_age, ego_names_verbose)

#
# ... Plot Mixing Matrices ...
#

par(mfrow=c(1,1))
M_comb_sex_age <- apply(extract(fit_degree_comb)$M, c(2,3), mean)
M_comb_weighted <- M_comb_sex_age
for(i in 1:6) { M_comb_weighted[i,] <- M_comb_weighted[i,] * pop_age_sex[i]}

# ... 1. Mixing Matrix for Sex and Age ... 

plot_mixing(M_comb_sex_age, 32, ego_names)

# ... 2. Mixing Matrix for Sex ... 

M_comb_sex <- matrix(0, nrow = 2, ncol = 8)
M_comb_sex[1,] <- colSums(M_comb_weighted[1:3,])/sum(pop_age_sex[1:3])
M_comb_sex[2,] <- colSums(M_comb_weighted[4:6,])/sum(pop_age_sex[4:6])

plot_mixing(M_comb_sex, 11, c("Male", "Female"))

# ... 3. Mixing Matrix for Age ... 

M_comb_age <- matrix(0, nrow = 3, ncol = 8)
M_comb_age[1,] <- colSums(M_comb_weighted[c(1,4),])/sum(pop_age_sex[c(1,4)])
M_comb_age[2,] <- colSums(M_comb_weighted[c(2,5),])/sum(pop_age_sex[c(2,5)])
M_comb_age[3,] <- colSums(M_comb_weighted[c(3,6),])/sum(pop_age_sex[c(3,6)])

plot_mixing(M_comb_age, 16, c("18-24", "25-64", "65+"))

