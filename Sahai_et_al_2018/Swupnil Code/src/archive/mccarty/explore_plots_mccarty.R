# ---------------------------------------
# ---------- Exploratory Plots ----------
# ---------------------------------------

# ... Response Histograms ....

par(mfrow=c(2,3))
for(i in 1:ncol(names_data_mccarty)){
  hist(names_data[-which(names_data_mccarty[,i]<0),i], col = "gray", xlab = "Number Known", xlim = c(0,40), main = colnames(names_data_mccarty)[i], breaks = c(0:150)*2)
}

#
#
# ... Missing Values Histograms ....

par(mfrow=c(2,1))
barplot(apply(data[,2:13], 2, function(x) length(which(x == -1))/length(x)*100), ylim = c(0,7), ylab = "Responses Missing (%)", main = "Missing Name Responses")
barplot(apply(data[,14:21], 2, function(x) length(which(x == -1))/length(x)*100), ylim = c(0,9), ylab = "Responses Missing (%)", main = "Missing Occupation Responses")
