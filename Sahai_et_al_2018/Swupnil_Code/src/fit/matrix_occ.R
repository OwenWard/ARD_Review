# ---------------------------------------------------------------------------
# ---------- Estimating Degree from Occupations Using Mixing Model ----------
# ---------------------------------------------------------------------------

require("rstan");
options(mc.cores = parallel::detectCores());

# -----------------------------
# ---------- Fitting ----------
# -----------------------------

# 1. Run the Models in Stan ....
mcmc_data_occ <- list(E = length(unique(ego_sex_age_race)), 
                      A = nrow(beta_occ), 
                      K = ncol(beta_occ), 
                      N = nrow(occ_data), 
                      y = occ_data, 
                      ego = ego_sex_age_race, 
                      Beta = beta_occ, 
                      theta_d = c(6.2,.5), 
                      theta_o = c(3,2), 
                      alpha = rep(1,24)/24);
fit_degree_occ <- sampling(degree_mixing_fit, data = mcmc_data_occ, iter = 100, chains = 4);

# --------------------------------
# ---------- Assessment ----------
# --------------------------------

degree$OccSexAgeRace <- round(colMeans(extract(fit_degree_occ)$d));

# ... Plot Degree Distributions ....

par(mfrow=c(2,3))
plot_degrees(degree$OccSexAgeRace, ego_sex_age, ego_names_verbose);
plot_degrees(degree$OccSexAgeRace, ego_sex_age_race, ego_names_verbose_race);

# ... Plot Mixing Matrix ...
M_occ_sex_age_race <- apply(extract(fit_degree_occ)$M, c(2,3), mean)
colnames(M_occ_sex_age_race) <- rownames(M_occ_sex_age_race) <- rownames(beta_occ)
M_occ_weighted <- M_occ_sex_age_race
for(i in 1:12) { M_occ_weighted[i,] <- M_occ_weighted[i,] * pop_age_sex_race[i]}

# 1. Mixing Matrix for Sex and Age and Race ... 
barplot(t(M_occ_sex_age_race[,1:12]), horiz = T, beside = T, xlim = c(-0.6, 0.6), ylim = c(0, 26))
barplot(t(-M_occ_sex_age_race[,12 + 1:12]), horiz = T, beside = T, add = T)
for(i in 1:6) {
  text(-0.5, 2.5 + 4*(i-1), rownames(beta_occ)[i])
  text(0.5, 2.5 + 4*(i-1) - 1.25, "18-24")
  text(0.5, 2.5 + 4*(i-1), "25-64", col = "darkgray")
  text(0.5, 2.5 + 4*(i-1) + 1.25, "65+", col = "lightgray")
  lines(c(-.45, -.45), 2.5 + 4*(i-1) + c(-1.5,1.5))
}
text(-.5, 26, "Ego Groups")
text(0.2, 26, "Male Alters")
text(-0.2, 26, "Female Alters")
text(.5, 26, "Alter Ages")

# 2. Mixing Matrix for Sex ... 
M_occ_sex <- matrix(0, nrow = 2, ncol = 12)
M_occ_sex[1,] <- colSums(M_occ_weighted[1:12,])/sum(pop_age_sex_race[1:12])
M_occ_sex[2,] <- colSums(M_occ_weighted[13:24,])/sum(pop_age_sex_race[13:24])

barplot(t(M_occ_sex[,1:12]), horiz = T, beside = T, xlim = c(-0.6, 0.6), ylim = c(0, 9))
barplot(t(-M_occ_sex[,12 + 1:12]), horiz = T, beside = T, add = T)
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

M_occ_age <- matrix(0, nrow = 3, ncol = 12)
M_occ_age[1,] <- colSums(M_occ_weighted[c(1,4,7,10,13,16,19,22),])/sum(pop_age_sex[c(1,4,7,10,13,16,19,22)])
M_occ_age[2,] <- colSums(M_occ_weighted[c(1,4,7,10,13,16,19,22) + 1,])/sum(pop_age_sex[c(1,4,7,10,13,16,19,22) + 1])
M_occ_age[3,] <- colSums(M_occ_weighted[c(1,4,7,10,13,16,19,22) + 2,])/sum(pop_age_sex[c(1,4,7,10,13,16,19,22) + 2])

barplot(t(M_occ_age[,1:12]), horiz = T, beside = T, xlim = c(-0.6, 0.6), ylim = c(0, 13))
barplot(t(-M_occ_age[,12 + 1:12]), horiz = T, beside = T, add = T)
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

# ... 4. Mixing Matrix for Race ... 

M_occ_race <- matrix(0, nrow = 4, ncol = 12)
M_occ_age[1,] <- colSums(M_occ_weighted[c(1:3, 13:15),])/sum(pop_age_sex[c(1:3, 13:15)])
M_occ_age[2,] <- colSums(M_occ_weighted[c(1:3, 13:15) + 3,])/sum(pop_age_sex[c(1:3, 13:15) + 3])
M_occ_age[3,] <- colSums(M_occ_weighted[c(1:3, 13:15) + 6,])/sum(pop_age_sex[c(1:3, 13:15) + 6])
M_occ_age[4,] <- colSums(M_occ_weighted[c(1:3, 13:15) + 9,])/sum(pop_age_sex[c(1:3, 13:15) + 9])

barplot(t(M_occ_age[,1:12]), horiz = T, beside = T, xlim = c(-0.6, 0.6), ylim = c(0, 13))
barplot(t(-M_occ_age[,12 + 1:12]), horiz = T, beside = T, add = T)
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