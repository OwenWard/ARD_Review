# -----------------------------------------
# ---------- Cleaned Survey Data ----------
# -----------------------------------------

mccarty <- read.csv("Data/mccarty.csv");
mccarty <- mccarty[,c(1:4,6,8,10,12,14,5,7,9,11,13,15)];
mccarty <- mccarty[-which(is.na(mccarty$Sex) | is.na(mccarty$Age)),];

for(i in 1:nrow(mccarty)) {
  current <- mccarty[i,4:15];
  if (length(which(current==0)) == 0 && length(which(is.na(current))) != length(current)) {
    current[which(is.na(current))] <- 0;
  } 
  else {
    current[which(is.na(current))] <- -1;
  } 
  mccarty[i,4:15] <- current;
  
  rm(current);
}

# ----------------------------------
# ---------- Derived Data ----------
# ----------------------------------

names_data_mccarty <- mccarty[,4:15];

gender_mccarty <- mccarty$Sex;
age_mccarty <- mccarty$Age;

# ... Separate Respondents into Ego Categories ....

ego_sex_age_mccarty <- rep(0, length(mccarty));

ego_sex_age_mccarty[(gender_mccarty == 1) & (age_mccarty >= 18) & (age_mccarty <= 24)] = 1;    #Male 18-24
ego_sex_age_mccarty[(gender_mccarty == 1) & (age_mccarty >= 25) & (age_mccarty <= 64)] = 2;   #Male 25-64
ego_sex_age_mccarty[(gender_mccarty == 1) & (age_mccarty >= 65)] = 3;  #Male 65+

ego_sex_age_mccarty[(gender_mccarty == 2) & (age_mccarty >= 18) & (age_mccarty <= 24)] = 4;    #Female 18-24
ego_sex_age_mccarty[(gender_mccarty == 2) & (age_mccarty >= 25) & (age_mccarty <= 64)] = 5;   #Female 25-64
ego_sex_age_mccarty[(gender_mccarty == 2) & (age_mccarty >= 65)] = 6;  #Female 65+

# ------------------------------------------
# ---------- Population Estimates ----------
# ------------------------------------------

pop_raw_mccarty <- read.csv("data/pop_age_2000.csv")[-c(1,87),];

# -------------------------------------
# ---------- Beta Adjustment ----------
# -------------------------------------

# ... Adjust Beta Matrix for Deaths (Multiply Name Birth Proportions By Number Alive Today) ...
beta_names_alive_mccarty <- beta_names_raw_mccarty
for(i in 1:nrow(beta_names_alive_mccarty)){
  for(j in 1:ncol(beta_names_alive_mccarty)){
    if(i <= 6) {
      beta_names_alive_mccarty[i,j] <- beta_names_raw_mccarty[i,j] * rev(pop_raw_mccarty$Male)[j]
    } else {
      beta_names_alive_mccarty[i,j] <- beta_names_raw_mccarty[i,j] * rev(pop_raw_mccarty$Female)[j]
    }
  }
}