par(mfrow=c(2,3))
plot(degree$NameSexAge,degree$CombSexAge,
     xlab="Matrix Name Sex Age",  
     ylab="Matrix Combined Sex Age",
     xlim = c(0,2500), ylim = c(0,2500), bty = 'l');
abline(a=0,b=1)
plot(degree$NameSexAge,degree$NameSexAgeKernel,
     xlab="Matrix Name Sex Age", 
     ylab="Kernel Name Sex Age",
     main="Degree Estimates",
     xlim = c(0,2500), ylim = c(0,2500), bty = 'l');
abline(a=0,b=1)
plot(degree$CombSexAge,degree$NameSexAgeKernel,
     xlab="Matrix Combined Sex Age", 
     ylab="Kernel Name Sex Age",
     xlim = c(0,2500), ylim = c(0,2500), bty = 'l');
abline(a=0,b=1)

plot(degree$OccSexAgeRace,degree$CombSexAge,
     xlab="Matrix Occupation Sex Age Race",
     ylab="Matrix Combined Sex Age",
     xlim = c(0,2500), ylim = c(0,2500), bty = 'l');
abline(a=0,b=1)
plot(degree$OccSexAgeRace,degree$OccSexAgeKernel,
     xlab="Matrix Occupation Sex Age Race",
     ylab="Kernel Occupation Sex Age",
     xlim = c(0,2500), ylim = c(0,2500), bty = 'l');
abline(a=0,b=1)
plot(degree$OccSexAgeKernel,degree$NameSexAgeKernel,
     xlab="Kernel Occupation Sex Age", 
     ylab="Kernel Name Sex Age",
     xlim = c(0,2500), ylim = c(0,2500), bty = 'l');
abline(a=0,b=1)