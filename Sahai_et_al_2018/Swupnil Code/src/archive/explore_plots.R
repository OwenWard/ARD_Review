# ---------------------------------------
# ---------- Exploratory Plots ----------
# ---------------------------------------

# ... Age Histograms ....

par(mfrow=c(2,3))
for(r in races){
  hist(as.numeric(age[race == r])*5 + 15, col = "gray", xlab = "Age", main = paste(r, "Age Distribution"), breaks = 14)
}

#
#
# ... Response Histograms ....

par(mfrow=c(2,3))
for(i in 1:ncol(names_data)){
  hist(names_data[-which(names_data[,i]<0),i], col = "gray", xlab = "Number Known", xlim = c(0,40), main = colnames(names_data)[i], breaks = c(0:150)*2)
}

par(mfrow=c(2,4))
for(i in 1:ncol(occ_data)){
  hist(occ_data[-which(occ_data[,i]<0),i], col = "gray", xlab = "Number Known", xlim = c(0,50), main = colnames(occ_data)[i], breaks = c(0:375)*2)
}

#
#
# ... Missing Values Histograms ....

par(mfrow=c(2,1))
barplot(apply(data[,2:13], 2, function(x) length(which(x == -1))/length(x)*100), ylim = c(0,7), ylab = "Responses Missing (%)", main = "Missing Name Responses")
barplot(apply(data[,14:21], 2, function(x) length(which(x == -1))/length(x)*100), ylim = c(0,9), ylab = "Responses Missing (%)", main = "Missing Occupation Responses")
