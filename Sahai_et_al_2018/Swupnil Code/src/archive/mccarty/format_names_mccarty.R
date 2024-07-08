# -----------------------------------------
# ---------- SSA Birth Names Data ---------
# -----------------------------------------

# ... Converts the SSA scraped names birth proportions into a nicely formatted (name x birthyear) matrix ... 

names <- colnames(mccarty)[4:15];
years <- 1917:2001;

# ... Read in the Raw Scraped Data ...

raw_names <- read.csv("Data/SSA-scraper/output.csv", header=FALSE);
colnames(raw_names) <- c("Year", "Rank", "Male_Name", "Male_Prop", "Female_Name", "Female_Prop");

# ... Create the Raw Beta Matrix ...
# This contains p(name | gender, age)

beta_names_raw_mccarty <- matrix(0, nrow = length(names), ncol = length(years));
rownames(beta_names_raw_mccarty) <- names;
colnames(beta_names_raw_mccarty) <- years;

for(j in 1:length(years)) {
  for(i in 1:length(names)) {
    current_prop <- 0
    current_male <- raw_names$Male_Prop[(raw_names$Male_Name == names[i]) & (raw_names$Year == years[j])];
    current_female <- raw_names$Female_Prop[(raw_names$Female_Name == names[i]) & (raw_names$Year == years[j])];
    
    if(length(current_male) > 0){
      current_prop <- current_prop + current_male;
    }
    if(length(current_female) > 0){
      current_prop <- current_prop + current_female;
    }
    
    beta_names_raw_mccarty[i,j] <- current_prop/100;
  }
}

rm(raw_names, names, years);
