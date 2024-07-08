# -----------------------------------------
# ---------- SSA Birth Names Data ---------
# -----------------------------------------

# ... Converts the SSA scraped names birth proportions into a nicely formatted (name x birthyear) matrix ... 

names <- colnames(names_data);
years <- 1915:2014;

# ... Read in the Raw Scraped Data ...

raw_names <- read.csv(paste0(code_path,"data/ssa/ssa_names.csv"), header = FALSE)
colnames(raw_names) <- c("Year", "Rank", "Male_Name", "Male_Prop", "Female_Name", "Female_Prop");

# ... Create the Raw Beta Matrix ...
# This contains p(name | gender, year born)

beta_names_raw <- matrix(0, nrow = length(names), ncol = length(years));
rownames(beta_names_raw) <- names;
colnames(beta_names_raw) <- years;

for(j in 1:length(years)) {
  for(i in 1:length(names)) {
    current_prop <- 0;
    current_male <- raw_names$Male_Prop[(raw_names$Male_Name == names[i]) & (raw_names$Year == years[j])];
    current_female <- raw_names$Female_Prop[(raw_names$Female_Name == names[i]) & (raw_names$Year == years[j])];
    
    if(length(current_male) > 0){
      current_prop <- current_prop + current_male;
    }
    if(length(current_female) > 0){
      current_prop <- current_prop + current_female;
    }
    
    beta_names_raw[i,j] <- current_prop/100;
  }
}

rm(raw_names, names, current_prop, current_male, current_female);

# -------------------------------------
# ---------- Beta Adjustment ----------
# -------------------------------------

# ... Adjust Beta Matrix for Deaths (Multiply Name Birth Proportions By Number Alive Today) ...
beta_names_alive <- beta_names_raw;
for(i in 1:nrow(beta_names_alive)){
  for(j in 1:ncol(beta_names_alive)){
    if(i > 6) {
      beta_names_alive[i,j] <- beta_names_raw[i,j] * rev(pop_raw$Male)[j];
    } 
    else {
      beta_names_alive[i,j] <- beta_names_raw[i,j] * rev(pop_raw$Female)[j];
    }
  }
}

# -------------------------------------------------------
# ---------- Derived Data: Name Age Categories ----------
# -------------------------------------------------------

beta_names <- matrix(0, ncol = ncol(names_data), nrow = ncol(pop_age_sex));
colnames(beta_names) <- colnames(names_data);
rownames(beta_names) <- colnames(pop_age_sex);

for(j in 1:ncol(beta_names)){
  if(j > 6) {
    beta_names[1,j] <- sum(beta_names_alive[j,83:99])/pop_age_sex[1];   #F 1-17 -- 1998-2014
    beta_names[2,j] <- sum(beta_names_alive[j,76:82])/pop_age_sex[2];   #M 18-24 -- 1991-1997
    beta_names[3,j] <- sum(beta_names_alive[j,36:75])/pop_age_sex[3];   #M 25-64 -- 1951-1990
    beta_names[4,j] <- sum(beta_names_alive[j,1:35])/pop_age_sex[4];    #M 65+ -- 1916-1950
  } 
  else {
    beta_names[5,j] <- sum(beta_names_alive[j,83:99])/pop_age_sex[5];   #F 1-17 -- 1998-2014
    beta_names[6,j] <- sum(beta_names_alive[j,76:82])/pop_age_sex[6];   #F 18-24 -- 1991-1997
    beta_names[7,j] <- sum(beta_names_alive[j,36:75])/pop_age_sex[7];   #F 25-64 -- 1951-1990
    beta_names[8,j] <- sum(beta_names_alive[j,1:35])/pop_age_sex[8];    #F 65+   -- 1916-1950
  }
}

# ----------------------------------------------------
# ---------- Derived Data: Name Age Moments ----------
# ----------------------------------------------------

age_max_2014 <- 99;
age_mean_2014 <- sum(c(0:99) * (pop_raw$Total)/(sum(pop_raw$Total)));

g_k <- g_k_names <- c(rep(1, 6), rep(2, 6));
x_axis <- 2016-c(1915:2014);
mu_k_name <- sigma_k_name <- sum_prob_k_name <- c();

for(k in 1:nrow(beta_names_alive)) {
  if(g_k[k] == 1) {
    # Female Alters: sum_(a_j) p(j in G_k | a_j, g_j)
    sum_prob_k_name[k] <- sum(beta_names_raw[k,]);
    y_axis <- (beta_names_raw[k,])/sum_prob_k_name[k];
  }
  else {
    # Male Alters: sum_(a_j) p(j in G_k | a_j, g_j)
    sum_prob_k_name[k] <- sum(beta_names_raw[k,]);
    y_axis <- (beta_names_raw[k,])/sum_prob_k_name[k];
  }
  # Ego Perspective: mu and sigma for p(j in G_k | a_j, g_j)/[sum_(a_j) p(j in G_k | a_j, g_j)]
  mu_k_name[k] <- x_axis%*%y_axis/sum(y_axis);
  sigma_k_name[k] <- sqrt((x_axis-mu_k_name[k])^2 %*% y_axis/sum(y_axis));
}

# Plot the age distributions of the names

par(mfrow=c(3,4))
for(i in 1:nrow(beta_names_alive)) {
  if(i > 8) {
    dens <- beta_names_alive[i,]/rev(pop_raw$Male);
  }
  else {
    dens <- beta_names_alive[i,]/rev(pop_raw$Female);
  }
  plot(2014-years,
       dens,
       ylab = "p(name|age,gender)",
       xlab = "Age",
       main = rownames(beta_names_alive)[i],
       bty = "l",
       type = 'l');
}
