slip <- read.csv("~/slip-model.csv")

chain <- slip[,1:3]

# sets the burnin size, removes all rows from the chain that are associated with the burnin 
burnin <- ceiling(0.1*(nrow(chain)-1))
acceptance <- sum(slip[-(1:burnin),"Accept"]) / nrow(slip)

med1 <- round(median(chain[-(1:burnin),1]),4)
med2 <- round(median(chain[-(1:burnin),2]),4)
med3 <- round(median(chain[-(1:burnin),3]),4)
par(mfrow=c(1,3), mar=c(5,5,4,1))
# PLOTTING 
hist(chain[-(1:burnin),1], main="Posterior of Enter", xlab="Prob(Enter)",ylab="Frequency",col="lightskyblue")
abline(v = med1, col='red',lwd=2)
#text(0.000168, 70000, paste0("Median = ", med1))
#text(0.000168, 60000, paste0("Acceptance = ", as.character(acceptance)))
hist(chain[-(1:burnin),2], main="Posterior of Stay", xlab="Prob(Stay)", ylab="Frequency",col="lightskyblue")
#text(0.905, 70000, paste0("Median = ", med2))
#abline(v = med2, col='red',lwd=2)
hist(chain[-(1:burnin),3], main="Posterior of Rate", xlab="Rate", ylab="Frequency",col="lightskyblue")
#text(0.00017, 400, paste0("Acceptance = ",as.character(acceptance)))
dev.off()

len <- length(chain[,1])
thinned <- seq(burnin, len, 100)

par(mfrow=c(1,2), mar=c(5,5,4,1))
plot(chain[thinned,1], type = "l", xlab="MCMC Steps" , ylab="Prob(Enter)",main = "Chain values of Enter")
abline(h = med1, col="red")
plot(chain[thinned,2], type = "l", xlab="MCMC Steps" , ylab="Prob(Stay)",main = "Chain values of Stay")
abline(h = med2, col="red")