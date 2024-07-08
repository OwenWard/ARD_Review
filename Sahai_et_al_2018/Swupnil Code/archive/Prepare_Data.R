# ... 1. Import Cleaned Survey Data ....
data_raw <- data <- read.csv("Data/Omni.csv")

for(i in 1:nrow(data)) {
  current <- data[i,2:13]
  if (length(which(current==0)) == 0 && length(which(is.na(current))) != length(current)) {
    current[which(is.na(current))] <- 0
  } else {
    current[which(is.na(current))] <- -1
  } 
  data[i,2:13] <- current
  
  current <- data[i,14:21]
  if (length(which(current==0)) == 0 && length(which(is.na(current))) != length(current)) {
    current[which(is.na(current))] <- 0
  } else {
    current[which(is.na(current))] <- -1
  } 
  data[i,14:21] <- current
  
  rm(current)
}

# ... Derived Data ...

names_data <- data[,2:13]
occ_data <- data[,14:21]

gender <- data$Gender
genders <- levels(gender)
age <- data$Age
ages <- levels(age)
race <- data$Race
races <- levels(race)

degree <- matrix(0, nrow = nrow(data), ncol = 0)
degree$ID <- data[,1]

# ... 2. Import National Population Estimates ....

pop_raw <- read.csv("data/pop_age_race.csv")
pop_raw <- pop_raw[-c(1:19, 102),-2]          #Keep Only 18-99 Year Olds
pop_raw$Male <- read.csv("data/pop_age.csv")[-c(1:19, 102), 3]
pop_raw$Female <- read.csv("data/pop_age.csv")[-c(1:19, 102), 4]

# ... 3. Compile the Model in Stan ....
degree_mixing_fit <- stan_model(model_code = degree_mixing_code, model_name = "Degree")

# ... 4. Exploratory Plots ....

# ... Age Histograms ....

par(mfrow=c(2,3))
for(r in races){
  hist(as.numeric(age[race == r])*5 + 15, col = "gray", xlab = "Age", main = paste(r, "Age Distribution"), breaks = 14)
}

# ... Response Histograms ....

par(mfrow=c(2,3))
for(i in 1:ncol(names_data)){
  hist(names_data[-which(names_data[,i]<0),i], col = "gray", xlab = "Number Known", xlim = c(0,40), main = colnames(names_data)[i], breaks = c(0:150)*2)
}

par(mfrow=c(2,4))
for(i in 1:ncol(occ_data)){
  hist(occ_data[-which(occ_data[,i]<0),i], col = "gray", xlab = "Number Known", xlim = c(0,50), main = colnames(occ_data)[i], breaks = c(0:375)*2)
}

# ... Missing Values Histograms ....

par(mfrow=c(2,1))
barplot(apply(data[,2:13], 2, function(x) length(which(x == -1))/length(x)*100), ylim = c(0,7), ylab = "Responses Missing (%)", main = "Missing Name Responses")
barplot(apply(data[,14:21], 2, function(x) length(which(x == -1))/length(x)*100), ylim = c(0,9), ylab = "Responses Missing (%)", main = "Missing Occupation Responses")
