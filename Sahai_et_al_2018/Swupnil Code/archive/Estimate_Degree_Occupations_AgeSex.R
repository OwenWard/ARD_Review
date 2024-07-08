require("rstan")
options(mc.cores = parallel::detectCores())

#
# ... Import National Data ...
#

# ... 1. Import National Population Estimates ....
pop_age_sex <- c(sum(pop_raw$Male[1:7]), sum(pop_raw$Male[8:47]), sum(pop_raw$Male[48:82]))
pop_age_sex <- c(pop_age_sex, sum(pop_raw$Female[1:7]), sum(pop_raw$Female[8:47]), sum(pop_raw$Female[48:82]))
pop_age_sex <- t(as.matrix(pop_age_sex))
colnames(pop_age_sex) <- c("M_18-24","M_25-64","M_65+","F_18-24","F_25-64","F_65+")

# ... 2. Import National Occupation Estimates ...
occ_raw <- read.csv("data/occ_sex_age_race.csv")
occ_codes <- c(2200, 4600, 3850, 2100, 2010, 6355, 4510, 4040)

#
# ... Prepare the Data ...
#

# ... 1. Separate Respondents into Ego Categories ....
ego_sex_age <- rep(0, length(data))

ego_sex_age[(gender == genders[2]) & (age %in% ages[1:2])] = 1    #Male 18-24
ego_sex_age[(gender == genders[2]) & (age %in% ages[3:10])] = 2   #Male 25-64
ego_sex_age[(gender == genders[2]) & (age %in% ages[11:13])] = 3  #Male 65+

ego_sex_age[(gender == genders[1]) & (age %in% ages[1:2])] = 4    #Female 18-24
ego_sex_age[(gender == genders[1]) & (age %in% ages[3:10])] = 5   #Female 25-64
ego_sex_age[(gender == genders[1]) & (age %in% ages[11:13])] = 6  #Female 65+

# ... 2. Group Occupations Into Alter Categories ...
beta_occ_sex_age <- matrix(ncol = ncol(occ_data), nrow = length(unique(ego_sex_age)))
colnames(beta_occ_sex_age) <- colnames(occ_data)
rownames(beta_occ_sex_age) <- colnames(pop_age_sex)

age_list <- list(c(18:24), c(45:64), c(65:100))
for(j in 1:ncol(beta_occ_sex_age)) {
  for(s in 0:1) {
    for(a in 0:2) {
      current_occ_sex_age <- occ_raw$weight[(occ_raw$occ == occ_codes[j]) & (occ_raw$sex == s+1) & (occ_raw$age %in% age_list[[a+1]])]
      current_sex_age <- occ_raw$weight[(occ_raw$sex == s+1) & (occ_raw$age %in% age_list[[a+1]])]
      beta_occ_sex_age[s*3 + a + 1,j] <- sum(current_occ_sex_age) / sum(current_sex_age)
    }
  }
}

#
# ... Estimate the Degrees Using Responses to Occupations ...
#

# ... 1. Run the Models in Stan ....
mcmc_data_occ_sex_age <- list(E = length(unique(ego_sex_age)), A = nrow(beta_occ_sex_age), K = ncol(beta_occ_sex_age), N = nrow(occ_data), y = occ_data, ego = ego_sex_age, Beta = beta_occ_sex_age, theta_d = c(6.2,.5), theta_o = c(3,2), alpha = rep(1,6)/6, p = 1)
fit_degree_occ_sex_age <- sampling(degree_mixing_fit, data = mcmc_data_occ_sex_age, iter = 1000, chains = 4)

degree$OccSexAge <- round(colMeans(extract(fit_degree_occ_sex_age)$d))

# ... 2. Plot Degree Distributions ....
par(mfrow=c(2,3))
for(i in 1:6) {
  subset_ego_i <- which(ego_sex_age == i)
  hist(degree$OccSexAge[subset_ego_i], xlab = "Degree", xlim = c(0,5000), main = rownames(beta_occ_sex_age)[i], breaks = c(0:20)*250)
  abline(v = mean(degree$OccSexAge[subset_ego_i]), lty = 2, col = "red")
  text(mean(degree$OccSexAge[subset_ego_i]), length(subset_ego_i)/5, paste("Mean = ",round(mean(degree$OccSexAge[subset_ego_i])), sep =""), col = "red", pos = 4)
  abline(v = median(degree$OccSexAge[subset_ego_i]), lty = 2, col = "blue")
  text(median(degree$OccSexAge[subset_ego_i]), length(subset_ego_i)/4, paste("Median = ",round(median(degree$OccSexAge[subset_ego_i])), sep =""), col = "blue", pos = 4)
}

#
# ... Plot Mixing Matrix ...
#

M_occ_sex_age <- apply(extract(fit_degree_occ_sex_age)$M, c(2,3), mean)
M_occ_sex_age_weighted <- M_occ_sex_age
for(i in 1:6) { M_occ_sex_age_weighted[i,] <- M_occ_sex_age_weighted[i,] * pop_age_sex[i]}

# ... 1. Mixing Matrix for Sex and Age ... 

barplot(t(M_occ_sex_age[,1:3]), horiz = T, beside = T, xlim = c(-0.6, 0.6), ylim = c(0, 26))
barplot(t(-M_occ_sex_age[,3 + 1:3]), horiz = T, beside = T, add = T)
for(i in 1:6) {
  text(-0.5, 2.5 + 4*(i-1), rownames(beta_occ_sex_age)[i])
  text(0.5, 2.5 + 4*(i-1) - 1.25, "18-24")
  text(0.5, 2.5 + 4*(i-1), "25-64", col = "darkgray")
  text(0.5, 2.5 + 4*(i-1) + 1.25, "65+", col = "lightgray")
  lines(c(-.42, -.42), 2.5 + 4*(i-1) + c(-1.5,1.5))
}
text(-.5, 26, "Ego Groups")
text(0.2, 26, "Male Alters")
text(-0.2, 26, "Female Alters")
text(.5, 26, "Alter Ages")

# ... 2. Mixing Matrix for Sex ... 

M_occ_sex <- matrix(0, nrow = 2, ncol = 6)
M_occ_sex[1,] <- colSums(M_occ_sex_age_weighted[1:3,])/sum(pop_age_sex[1:3])
M_occ_sex[2,] <- colSums(M_occ_sex_age_weighted[4:6,])/sum(pop_age_sex[4:6])

barplot(t(M_occ_sex[,1:3]), horiz = T, beside = T, xlim = c(-0.6, 0.6), ylim = c(0, 9))
barplot(t(-M_occ_sex[,3 + 1:3]), horiz = T, beside = T, add = T)
for(i in 1:2) {
  text(-0.5, 2.5 + 4*(i-1), c("Male", "Female")[i])
  text(0.5, 2.5 + 4*(i-1) - 1.25, "18-24")
  text(0.5, 2.5 + 4*(i-1), "25-64", col = "darkgray")
  text(0.5, 2.5 + 4*(i-1) + 1.25, "65+", col = "lightgray")
  lines(c(-.42, -.42), 2.5 + 4*(i-1) + c(-1.5,1.5))
}
text(-.5, 9, "Ego Groups")
text(0.2, 9, "Male Alters")
text(-0.2, 9, "Female Alters")
text(.5, 9, "Alter Ages")

# ... 3. Mixing Matrix for Age ... 

M_occ_age <- matrix(0, nrow = 3, ncol = 6)
M_occ_age[1,] <- colSums(M_occ_sex_age_weighted[c(1,4),])/sum(pop_age_sex[c(1,4)])
M_occ_age[2,] <- colSums(M_occ_sex_age_weighted[c(2,5),])/sum(pop_age_sex[c(2,5)])
M_occ_age[3,] <- colSums(M_occ_sex_age_weighted[c(3,6),])/sum(pop_age_sex[c(3,6)])

barplot(t(M_occ_age[,1:3]), horiz = T, beside = T, xlim = c(-0.6, 0.6), ylim = c(0, 13))
barplot(t(-M_occ_age[,3 + 1:3]), horiz = T, beside = T, add = T)
for(i in 1:3) {
  text(-0.5, 2.5 + 4*(i-1), c("18-24", "25-64", "65+")[i])
  text(0.5, 2.5 + 4*(i-1) - 1.25, "18-24")
  text(0.5, 2.5 + 4*(i-1), "25-64", col = "darkgray")
  text(0.5, 2.5 + 4*(i-1) + 1.25, "65+", col = "lightgray")
  lines(c(-.42, -.42), 2.5 + 4*(i-1) + c(-1.5,1.5))
}
text(-.5, 13, "Ego Groups")
text(0.2, 13, "Male Alters")
text(-0.2, 13, "Female Alters")
text(.5, 13, "Alter Ages")
