# -------------------------------------
# ---------- Occupation Data ----------
# -------------------------------------

# ... Converts the occupation survey data into a nicely formatted (occ x birthyear) matrix ... 

occs <- colnames(occ_data);
years <- 1915:2014;

# ... Read in the Survey Data ...

occ_raw <- read.csv(paste0(code_path, "data/omni/occ_sex_age_race.csv"),
                    stringsAsFactors = TRUE)
occ_codes <- c(2200, 4600, 3850, 2100, 2010, 6355, 4510, 4040);
g_k <- c(rep(1, 8), rep(2, 8));

# ... Create the Raw Beta Matrix ...
# This contains p(occupation | gender, age)

beta_occ_raw <- matrix(0, nrow = length(occs) * 2, ncol = length(years));
rownames(beta_occ_raw) <- c(paste("Female",occs),paste("Male",occs));
colnames(beta_occ_raw) <- years;

I <- length(occ_codes);
for(j in 1:ncol(beta_occ_raw)){
  for(i in 1:nrow(beta_occ_raw)){
    age_ <- 2014 - years[j];
    sex_ <- 3 - g_k[i];
    occ_ <- occ_codes[(i-1)%%I + 1];
    occ_sex_age_ <- occ_raw$weight[(occ_raw$occ == occ_) & (occ_raw$sex == sex_) & (occ_raw$age == age_)];
    sex_age_ <- occ_raw$weight[(occ_raw$sex == sex_) & (occ_raw$age == age_)];
    if(sum(sex_age_) > 0) {
      beta_occ_raw[i,j] <- sum(occ_sex_age_) / sum(sex_age_);
    }
  }
  rm(age_, sex_, occ_, occ_sex_age_, sex_age_);
}

# -------------------------------------
# ---------- Beta Adjustment ----------
# -------------------------------------

# ... Adjust Beta Matrix for Deaths (Multiply Name Birth Proportions By Number Alive Today) ...
beta_occ_alive <- beta_occ_raw;
for(i in 1:nrow(beta_occ_alive)){
  for(j in 1:ncol(beta_occ_alive)){
    if(i > 8) {
      beta_occ_alive[i,j] <- beta_occ_raw[i,j] * rev(pop_raw$Male)[j];
    } 
    else {
      beta_occ_alive[i,j] <- beta_occ_raw[i,j] * rev(pop_raw$Female)[j];
    }
  }
}

# -----------------------------------
# ---------- Beta Smoothed ----------
# -----------------------------------

# ... Adjust Smooth Beta Matrix To Adjust Zero Values ...
beta_occ_smooth <- beta_occ_raw;
ages_ <- 2014 - years;
par(mfrow=c(1,2));
for(i in 1:nrow(beta_occ_raw)){
  sex_ <- 3 - g_k[i];
  occ_ <- occ_codes[(i-1)%%I + 1];
  dens <- density(occ_raw$age[occ_raw$occ==occ_ & occ_raw$sex==sex_], 
                  weights = occ_raw$weight[occ_raw$occ==occ_ & occ_raw$sex==sex_]);
  dens_binned <- rep(0, ncol(beta_occ_raw));
  for(j in 1:length(dens_binned)) {
    current_ <- which((dens$x > ages_[j] - 0.5) & (dens$x <= ages_[j] + 0.5));
    if(length(current_) > 0) {
      dens_binned[j] <- mean(dens$y[current_]);
    }
  }
  sum_dens <- sum(dens_binned[beta_occ_raw[i,] > 0]);
  sum_beta <- sum(beta_occ_raw[i,]);
  for(j in 1:ncol(beta_occ_raw)){
    if(i > 8) {
      beta_occ_smooth[i,j] <- dens_binned[j] * sum_beta / sum_dens;
    } 
    else {
      beta_occ_smooth[i,j] <- dens_binned[j] * sum_beta / sum_dens;
    }
  }
  rm(dens, sex_, dens_binned, current_, sum_dens, sum_beta);
}

# -------------------------------------------------------
# ---------- Derived Data: Name Age Categories ----------
# -------------------------------------------------------

# ... Group Occupations Into Alter Categories ...

beta_occ <- matrix(ncol = ncol(occ_data), nrow = length(unique(ego_sex_age_race)));
colnames(beta_occ) <- colnames(occ_data);
rownames(beta_occ) <- colnames(pop_age_sex_race);

age_list <- list(c(18:24), c(45:64), c(65:100))
for(j in 1:ncol(beta_occ)) {
  for(s in 0:1) {
    for(a in 0:2) {
      current_occ_sex_race_age <- occ_raw$weight[(occ_raw$occ == occ_codes[j]) & (occ_raw$sex == s+1) & (occ_raw$age %in% age_list[[a+1]]) & (occ_raw$hisp != 0)]
      current_sex_race_age <- occ_raw$weight[(occ_raw$sex == 1) & (occ_raw$age %in% age_list[[a+1]]) & (occ_raw$hisp != 0)]
      beta_occ[s*12 + a + 1,j] <- sum(current_occ_sex_race_age) / sum(current_sex_race_age)
      
      current_occ_sex_race_age <- occ_raw$weight[(occ_raw$occ == occ_codes[j]) & (occ_raw$sex == s+1) & (occ_raw$age %in% age_list[[a+1]]) & (occ_raw$race == 100) & (occ_raw$hisp == 0)]
      current_sex_race_age <- occ_raw$weight[(occ_raw$sex == 1) & (occ_raw$age %in% age_list[[a+1]]) & (occ_raw$race == 100) & (occ_raw$hisp == 0)]
      beta_occ[s*12 + a + 4,j] <- sum(current_occ_sex_race_age) / sum(current_sex_race_age)
      
      current_occ_sex_race_age <- occ_raw$weight[(occ_raw$occ == occ_codes[j]) & (occ_raw$sex == s+1) & (occ_raw$age %in% age_list[[a+1]]) & (occ_raw$race == 200) & (occ_raw$hisp == 0)]
      current_sex_race_age <- occ_raw$weight[(occ_raw$sex == 1) & (occ_raw$age %in% age_list[[a+1]]) & (occ_raw$race == 200) & (occ_raw$hisp == 0)]
      beta_occ[s*12 + a + 7,j] <- sum(current_occ_sex_race_age) / sum(current_sex_race_age)
      
      current_occ_sex_race_age <- occ_raw$weight[(occ_raw$occ == occ_codes[j]) & (occ_raw$sex == s+1) & (occ_raw$age %in% age_list[[a+1]]) & (occ_raw$race != 100) & (occ_raw$race != 200) & (occ_raw$hisp == 0)]
      current_sex_race_age <- occ_raw$weight[(occ_raw$sex == 1) & (occ_raw$age %in% age_list[[a+1]]) & (occ_raw$race != 100) & (occ_raw$race != 200) & (occ_raw$hisp == 0)]
      beta_occ[s*12 + a + 10,j] <- sum(current_occ_sex_race_age) / sum(current_sex_race_age)
    }
  }
}

# ----------------------------------------------------
# ---------- Derived Data: Occs Age Moments ----------
# ----------------------------------------------------

age_max_2014 <- 99;
age_mean_2014 <- sum(c(0:99) * (pop_raw$Total)/(sum(pop_raw$Total)));

g_k <- g_k_occ <- c(rep(1, 8), rep(2, 8));
x_axis <- 2016-c(1915:2014);
mu_k_occ <- sigma_k_occ <- sum_prob_k_occ <- c();

for(k in 1:nrow(beta_occ_alive)) {
  if(g_k[k] == 1) {
    # Female Alters: sum_(a_j) p(j in G_k | a_j, g_j)
    sum_prob_k_occ[k] <- sum(beta_occ_smooth[k,]);
    y_axis <- (beta_occ_smooth[k,])/sum_prob_k_occ[k];
  }
  else {
    # Male Alters: sum_(a_j) p(j in G_k | a_j, g_j)
    sum_prob_k_occ[k] <- sum(beta_occ_smooth[k,]);
    y_axis <- (beta_occ_smooth[k,])/sum_prob_k_occ[k];
  }
  # Ego Perspective: mu and sigma for p(j in G_k | a_j, g_j)/[sum_(a_j) p(j in G_k | a_j, g_j)]
  mu_k_occ[k] <- x_axis%*%y_axis/sum(y_axis);
  sigma_k_occ[k] <- sqrt((x_axis-mu_k_occ[k])^2 %*% y_axis/sum(y_axis));
}

# Plot the age distributions of the occupations using weighted density plots

par(mfrow=c(4,2));
for(i in 1:length(occs)) {
  male_ <- beta_occ_smooth[i+length(occs),];
  female_ <- beta_occ_smooth[i,];
  if(sum(male_) > sum(female_)) {
    large_color <- 1; small_color <- 2;
    large_ <- male_; small_ <- female_;
  }
  else {
    large_color <- 2; small_color <- 1;
    large_ <- female_; small_ <- male_;
  }
  plot(2014-years, large_,
       main = occs[i], ylab = "p(G_k | age, gender)",
       xlab = "Age", col = large_color,
       type = 'l', bty = "l");
  lines(2014-years, small_, col = small_color);
  rm(small_,large_,small_color,large_color,male_,female_)
}
