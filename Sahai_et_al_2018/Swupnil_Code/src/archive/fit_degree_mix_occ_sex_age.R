require("rstan")
options(mc.cores = parallel::detectCores())

#
# ... Prepare the Raw Data ...
#
# ... 1. Import National Occupation Estimates ...

occ_raw <- read.csv("data/occ_sex_age_race.csv")
occ_codes <- c(2200, 4600, 3850, 2100, 2010, 6355, 4510, 4040)

#
# ... Prepare the Data ...
#

# ... Group Occupations Into Alter Categories ...

beta_occ_sex_age <- matrix(0, ncol = ncol(occ_data), nrow = ncol(pop_age_sex))
colnames(beta_occ_sex_age) <- colnames(occ_data)
rownames(beta_occ_sex_age) <- colnames(pop_age_sex)
  
age_list <- list(c(1:17), c(18:24), c(45:64), c(65:99))
for(j in 1:ncol(beta_occ_sex_age)) {
  for(s in 0:1) {
    for(a in 0:3) {
      current_occ_sex_age <- occ_raw$weight[(occ_raw$occ == occ_codes[j]) & (occ_raw$sex == s+1) & (occ_raw$age %in% age_list[[a+1]])]
      current_sex_age <- occ_raw$weight[(occ_raw$sex == s+1) & (occ_raw$age %in% age_list[[a+1]])]
      beta_occ_sex_age[s*4 + a + 1,j] <- sum(current_occ_sex_age) / sum(current_sex_age)
    }
  }
}

#
# ... Estimate the Degrees Using Responses to Occupations ...
#
# ... 1. Run the Models in Stan ....

rm_occ <- c(2)
mcmc_data_occ_sex_age <- list(E = length(unique(ego_sex_age)), 
                              A = nrow(beta_occ_sex_age), 
                              K = ncol(beta_occ_sex_age[,-rm_occ]), 
                              N = nrow(occ_data), 
                              y = occ_data[,-rm_occ], 
                              ego = ego_sex_age, 
                              Beta = beta_occ_sex_age[,-rm_occ], 
                              theta_d = c(6.2,.5), 
                              theta_o = c(3,2), 
                              alpha = rep(1,nrow(beta_occ_sex_age))/nrow(beta_occ_sex_age))
fit_degree_occ_sex_age <- sampling(degree_mixing_fit, data = mcmc_data_occ_sex_age, iter = 1000, chains = 4)

degree$OccSexAge <- round(colMeans(extract(fit_degree_occ_sex_age)$d))
omega_occ_sex_age <- colMeans(extract(fit_degree_occ_sex_age)$inv_omega)

# ... 2. Plot Degree Distributions ....

par(mfrow=c(2,3))
plot_degrees(degree$OccSexAge, ego_sex_age, ego_names_verbose)

#
# ... Plot Mixing Matrix ...
#

par(mfrow=c(1,1))
M_occ_sex_age <- apply(extract(fit_degree_occ_sex_age)$M, c(2,3), mean)
M_occ_sex_age_weighted <- M_occ_sex_age
for(i in 1:6) { M_occ_sex_age_weighted[i,] <- M_occ_sex_age_weighted[i,] * pop_age_sex[i]}

# ... 1. Mixing Matrix for Sex and Age ... 

plot_mixing(M_occ_sex_age, 32, ego_names)

# ... 2. Mixing Matrix for Sex ... 

M_occ_sex <- matrix(0, nrow = 2, ncol = 8)
M_occ_sex[1,] <- colSums(M_occ_sex_age_weighted[1:3,])/sum(pop_age_sex[1:3])
M_occ_sex[2,] <- colSums(M_occ_sex_age_weighted[4:6,])/sum(pop_age_sex[4:6])

plot_mixing(M_occ_sex, 11, c("Male", "Female"))

# ... 3. Mixing Matrix for Age ... 

M_occ_age <- matrix(0, nrow = 3, ncol = 8)
M_occ_age[1,] <- colSums(M_occ_sex_age_weighted[c(1,4),])/sum(pop_age_sex[c(1,4)])
M_occ_age[2,] <- colSums(M_occ_sex_age_weighted[c(2,5),])/sum(pop_age_sex[c(2,5)])
M_occ_age[3,] <- colSums(M_occ_sex_age_weighted[c(3,6),])/sum(pop_age_sex[c(3,6)])

plot_mixing(M_occ_age, 16, c("18-24", "25-64", "65+"))

