# -----------------------------------------
# ---------- Cleaned Survey Data ----------
# -----------------------------------------

data_raw <- data <- read.csv(paste0(code_path,"data/omni/omni.csv"),stringsAsFactors = TRUE)
# Replace NAs for names, occupations, and cancers
for (i in 1:nrow(data)) {
  segments <- list(2:13, 14:21, 45:49, 50:54)

  for (s in 1:length(segments)) {
    seg <- segments[[s]]
    current <- data[i, seg]

    if (length(which(current == 0)) == 0 && length(which(is.na(current))) != length(current)) {
      current[which(is.na(current))] <- 0
    } else {
      current[which(is.na(current))] <- -1
    }

    data[i, seg] <- current
  }
}

# ----------------------------------
# ---------- Derived Data ----------
# ----------------------------------

names_data <- data[, 2:13]
occ_data <- data[, 14:21]
allergy_data <- data[, 45:49]
cancer_data <- data[, 50:54]
gender <- data$Gender
genders <- levels(gender)
age <- data$Age
ages <- levels(age)
race <- data$Race
races <- levels(race)
rm_na_names <- which(rowSums(names_data) < 0)
rm_na_occs <- which(rowSums(occ_data) < 0)
rm_na_comb <- which(rowSums(occ_data) + rowSums(names_data) < 0)
rm_old <- which(data$age > 70)
rm_old_na_names <- which((data$age > 70) | (rowSums(names_data) < 0))
rm_old_na_occs <- which((data$age > 70) | (rowSums(occ_data) < 0))
rm_old_na_comb <- which((data$age > 70) | (rowSums(occ_data) + rowSums(names_data) < 0))
names_data_woOld <- names_data
names_data_woOld[rm_old, ] <- -1
occ_data_woOld <- occ_data
occ_data_woOld[rm_old, ] <- -1

## check this, does this need to be changed?

degree <- matrix(0, nrow = nrow(data), ncol = 0)
degree$ID <- data[, 1]
# ... Separate Respondents into Ego Categories ....

ego_sex_age <- rep(0, nrow(data))
ego_sex_age[(gender == genders[2]) & (age %in% ages[1:2])] <- 1
# Male 18-24
ego_sex_age[(gender == genders[2]) & (age %in% ages[3:10])] <- 2
# Male 25-64
ego_sex_age[(gender == genders[2]) & (age %in% ages[11:13])] <- 3
# Male 65+
ego_sex_age[(gender == genders[1]) & (age %in% ages[1:2])] <- 4
# Female 18-24
ego_sex_age[(gender == genders[1]) & (age %in% ages[3:10])] <- 5
# Female 25-64
ego_sex_age[(gender == genders[1]) & (age %in% ages[11:13])] <- 6
# Female 65+

ego_sex_age_race <- rep(0, nrow(data))
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[1:2]) & (race == races[3])] <- 1 # Male Hispanic 18-24
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[3:10]) & (race == races[3])] <- 2 # Male Hispanic 25-64
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[11:13]) & (race == races[3])] <- 3 # Male Hispanic 65+
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[1:2]) & (race == races[5])] <- 4 # Male White 18-24
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[3:10]) & (race == races[5])] <- 5 # Male White 25-64
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[11:13]) & (race == races[5])] <- 6 # Male White 65+
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[1:2]) & (race == races[2])] <- 7 # Male Black 18-24
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[3:10]) & (race == races[2])] <- 8 # Male Black 25-64
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[11:13]) & (race == races[2])] <- 9 # Male Black 65+
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[1:2]) & (race %in% races[c(1, 4)])] <- 10 # Male Other 18-24
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[3:10]) & (race %in% races[c(1, 4)])] <- 11 # Male Other 25-64
ego_sex_age_race[(gender == genders[2]) & (age %in% ages[11:13]) & (race %in% races[c(1, 4)])] <- 12 # Male Other 65+
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[1:2]) & (race == races[3])] <- 13 # Female Hispanic 18-24
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[3:10]) & (race == races[3])] <- 14 # Female Hispanic 25-64
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[11:13]) & (race == races[3])] <- 15 # Female Hispanic 65+
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[1:2]) & (race == races[5])] <- 16 # Female White 18-24
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[3:10]) & (race == races[5])] <- 17 # Female White 25-64
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[11:13]) & (race == races[5])] <- 18 # Female White 65+
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[1:2]) & (race == races[2])] <- 19 # Female Black 18-24
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[3:10]) & (race == races[2])] <- 20 # Female Black 25-64
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[11:13]) & (race == races[2])] <- 21 # Female Black 65+
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[1:2]) & (race %in% races[c(1, 4)])] <- 22 # Female Other 18-24
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[3:10]) & (race %in% races[c(1, 4)])] <- 23 # Female Other 25-64
ego_sex_age_race[(gender == genders[1]) & (age %in% ages[11:13]) & (race %in% races[c(1, 4)])] <- 24 # Female Other 65+

ego_names <- c("M_18-24", "M_25-64", "M_65+", "F_18-24", "F_25-64", "F_65+")
ego_names_verbose <- c("Male 18-24", "Male 25-64", "Male 65+", "Female 18-24", "Female 25-64", "Female 65+")
ego_names_verbose_race <- c(
  "Male Hispanic 18-24", "Male Hispanic 25-64", "Male Hispanic 65+", "Male White 18-24", "Male White 25-64", "Male White 65+",
  "Male Black 18-24", "Male Black 25-64", "Male Black 65+", "Male Other 18-24", "Male Other 25-64", "Male Other 65+",
  "Female Hispanic 18-24", "Female Hispanic 25-64", "Female Hispanic 65+", "Female White 18-24", "Female White 25-64", "Female White 65+",
  "Female Black 18-24", "Female Black 25-64", "Female Black 65+", "Female Other 18-24", "Female Other 25-64", "Female Other 65+"
)
# ------------------------------------------
# ---------- Population Estimates ----------
# ------------------------------------------

pop_raw <- read.csv(paste0(code_path, "data/pop/pop_age_2013.csv"), 
                    stringsAsFactors = TRUE)[-c(1, 102), ]
# Keep Only 0-99 Year Olds

pop_age_sex <- c(sum(pop_raw$Male[2:18]), sum(pop_raw$Male[19:25]), sum(pop_raw$Male[26:65]), sum(pop_raw$Male[66:100]))
pop_age_sex <- c(pop_age_sex, sum(pop_raw$Female[2:18]), sum(pop_raw$Female[19:25]), sum(pop_raw$Female[26:65]), sum(pop_raw$Female[66:100]))
pop_age_sex <- t(as.matrix(pop_age_sex))
colnames(pop_age_sex) <- c("M_1-17", "M_18-24", "M_25-64", "M_65+", "F_1-17", "F_18-24", "F_25-64", "F_65+")
pop_age_sex_race <- c(sum(pop_raw[1:7, 2]), sum(pop_raw[8:47, 2]), sum(pop_raw[48:82, 2]))
for (i in 3:9) {
  pop_age_sex_race <- c(pop_age_sex_race, sum(pop_raw[1:7, i]), sum(pop_raw[8:47, i]), sum(pop_raw[48:82, i]))
}
pop_age_sex_race <- t(as.matrix(pop_age_sex_race))
colnames(pop_age_sex_race) <- c(
  "M_H_18-24", "M_H_25-64", "M_H_65+", "M_W_18-24", "M_W_25-64", "M_W_65+",
  "M_B_18-24", "M_B_25-64", "M_B_65+", "M_O_18-24", "M_O_25-64", "M_O_65+",
  "F_H_18-24", "F_H_25-64", "F_H_65+", "F_W_18-24", "F_W_25-64", "F_W_65+",
  "F_B_18-24", "F_B_25-64", "F_B_65+", "F_O_18-24", "F_O_25-64", "F_O_65+"
)
