require("rstan")
options(mc.cores = parallel::detectCores())

#
# ... Prepare the Raw Data ...
#
# ... 1. Import Historical National Name Birth Proportions ...

beta_names_raw <- beta_names_raw[,-c(1:16)]   #only keep 1916 onwards

#
# ... Prepare the Derived Data ...
#
# ... 1. Separate Respondents into Ego Categories ....

ego_sex_age <- rep(0, length(data))

ego_sex_age[(gender == genders[2]) & (age %in% ages[1:2])] = 1    #Male 18-24
ego_sex_age[(gender == genders[2]) & (age %in% ages[3:10])] = 2   #Male 25-64
ego_sex_age[(gender == genders[2]) & (age %in% ages[11:13])] = 3  #Male 65+

ego_sex_age[(gender == genders[1]) & (age %in% ages[1:2])] = 4    #Female 18-24
ego_sex_age[(gender == genders[1]) & (age %in% ages[3:10])] = 5   #Female 25-64
ego_sex_age[(gender == genders[1]) & (age %in% ages[11:13])] = 6  #Female 65+

ego_names <- c("M_18-24","M_25-64","M_65+","F_18-24","F_25-64","F_65+")

# ... 2. Adjust Beta Matrix for Deaths (Multiply Name Birth Proportions By Number Alive Today) ...

beta_names_alive <- beta_names_raw
for(i in 1:nrow(beta_names_alive)){
  for(j in 1:ncol(beta_names_alive)){
    if(i > 6) {
      beta_names_alive[i,j] <- beta_names_raw[i,j] * rev(pop_raw$Male)[j]
    } else {
      beta_names_alive[i,j] <- beta_names_raw[i,j] * rev(pop_raw$Female)[j]
    }
  }
}

# ... 3. Group Names Into Alter Categories ...

beta_names <- matrix(0, ncol = ncol(names_data), nrow = ncol(pop_age_sex))
colnames(beta_names) <- colnames(names_data)
rownames(beta_names) <- colnames(pop_age_sex)

for(j in 1:ncol(beta_names)){
  if(j > 6) {
    beta_names[1,j] <- sum(beta_names_alive[j,83:99])/pop_age_sex[1]   #F 1-17 -- 1998-2014
    beta_names[2,j] <- sum(beta_names_alive[j,76:82])/pop_age_sex[2]   #M 18-24 -- 1991-1997
    beta_names[3,j] <- sum(beta_names_alive[j,36:75])/pop_age_sex[3]   #M 25-64 -- 1951-1990
    beta_names[4,j] <- sum(beta_names_alive[j,1:35])/pop_age_sex[4]    #M 65+ -- 1916-1950
  } else {
    beta_names[5,j] <- sum(beta_names_alive[j,83:99])/pop_age_sex[5]   #F 1-17 -- 1998-2014
    beta_names[6,j] <- sum(beta_names_alive[j,76:82])/pop_age_sex[6]   #F 18-24 -- 1991-1997
    beta_names[7,j] <- sum(beta_names_alive[j,36:75])/pop_age_sex[7]   #F 25-64 -- 1951-1990
    beta_names[8,j] <- sum(beta_names_alive[j,1:35])/pop_age_sex[8]    #F 65+   -- 1916-1950
  }
}

#
# ... Estimate the Degrees Using Responses to Names ...
#
# ... 1. Run the Model in Stan ....

rm_names <- c(1,2,7,8)
mcmc_data_names <- list(E = length(unique(ego_sex_age)), A = nrow(beta_names), K = ncol(beta_names[,-rm_names]), 
                        N = nrow(names_data), y = names_data[,-rm_names], ego = ego_sex_age, 
                        Beta = beta_names[,-rm_names], theta_d = c(6.2,.5), theta_o = c(3,2), 
                        alpha = rep(1,nrow(beta_names))/nrow(beta_names), p = 1)
fit_degree_names <- sampling(degree_mixing_fit, data = mcmc_data_names, iter = 1000, chains = 4)

degree$NameSexAge <- round(colMeans(extract(fit_degree_names)$d))
omega_names_sex_age <- colMeans(extract(fit_degree_names)$inv_omega)

# ... 2. Plot Degree Distributions ....

par(mfrow=c(2,3))
for(i in 1:6) {
  subset_ego_i <- which(ego_sex_age == i)
  hist(degree$NameSexAge[subset_ego_i], xlab = "Degree", xlim = c(0,5000), main = rownames(beta_names)[i], breaks = c(0:20)*250)
  abline(v = mean(degree$NameSexAge[subset_ego_i]), lty = 2, col = "red")
  text(mean(degree$NameSexAge[subset_ego_i]), length(subset_ego_i)/5, paste("Mean = ",round(mean(degree$NameSexAge[subset_ego_i])), sep =""), col = "red", pos = 4)
  abline(v = median(degree$NameSexAge[subset_ego_i]), lty = 2, col = "blue")
  text(median(degree$NameSexAge[subset_ego_i]), length(subset_ego_i)/4, paste("Median = ",round(median(degree$NameSexAge[subset_ego_i])), sep =""), col = "blue", pos = 4)
}

#
# ... Plot Mixing Matrices ...
#

par(mfrow=c(1,1))
M_names_sex_age <- apply(extract(fit_degree_names)$M, c(2,3), mean)
M_names_weighted <- M_names_sex_age
for(i in 1:6) { M_names_weighted[i,] <- M_names_weighted[i,] * pop_age_sex[i]}

# ... 1. Mixing Matrix for Sex and Age ... 

plot_mixing(M_names_sex_age, 32, ego_names)

# ... 2. Mixing Matrix for Sex ... 

M_names_sex <- matrix(0, nrow = 2, ncol = 8)
M_names_sex[1,] <- colSums(M_names_weighted[1:3,])/sum(pop_age_sex[1:3])
M_names_sex[2,] <- colSums(M_names_weighted[4:6,])/sum(pop_age_sex[4:6])

plot_mixing(M_names_sex, 11, c("Male", "Female"))

# ... 3. Mixing Matrix for Age ... 

M_names_age <- matrix(0, nrow = 3, ncol = 8)
M_names_age[1,] <- colSums(M_names_weighted[c(1,4),])/sum(pop_age_sex[c(1,4)])
M_names_age[2,] <- colSums(M_names_weighted[c(2,5),])/sum(pop_age_sex[c(2,5)])
M_names_age[3,] <- colSums(M_names_weighted[c(3,6),])/sum(pop_age_sex[c(3,6)])

plot_mixing(M_names_age, 16, c("18-24", "25-64", "65+"))

