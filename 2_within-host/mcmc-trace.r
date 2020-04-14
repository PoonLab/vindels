slip <- read.csv("~/PycharmProjects/hiv-withinhost/slip-model.csv")

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

# ----- SLIPPAGE MODEL -----
csv <- read.csv("~/PycharmProjects/hiv-withinhost/slip-model-simple-indel.csv", stringsAsFactors = F, skip=1, header=F)

# csv <- as.data.frame(sapply(1:ncol(csv), function(x){
#   sapply(1:nrow(csv), function(y){
#     as.numeric(strsplit(csv[y,x]," ")[[1]][1])
#   })
# }))
tru <- c(0.001, 0.70, 0.0001)
colnames(csv) <- c('p.enter', 'p.stay', 'rate', 'slip.changed', 'accept', 'time')
# colnames(csv) <- c('p.enter', 'slope','int', 'rate', 'slip.changed', 'accept', 'time')
#png(file="~/vindels/Figures/within-host/finalized/slippage-trace2.png", width=600, height=800)
par(mar=c(2.5,5,3,1), mfrow=c(3,1))
plot(csv$p.enter, type = "l", xlab="MCMC Steps (x10)" , ylab="Prob(Enter)",
     main = "Chain values of Enter", cex.axis=1.3, cex.lab=1.4, cex.main=1.7)
abline(h=tru[1],col='red',lwd=2)
plot(csv$p.stay, type = "l", xlab="MCMC Steps (x10)" , ylab="P(Stay)",
     main = "Chain values of Stay",  cex.axis=1.3, cex.lab=1.4, cex.main=1.7)
abline(h=tru[2],col='red',lwd=2)
plot(csv$rate, type = "l", xlab="MCMC Steps (x10)" , ylab="Rate",
     main = "Chain values of Rate",  cex.axis=1.3, cex.lab=1.4, cex.main=1.7)
abline(h=tru[3],col='red',lwd=2)
tru <- c(0.001, 25, 0.5, 0.0005)

med1 <- median(csv$p.enter)
med2 <- median(csv$p.stay)
med3 <- median(csv$rate)
#dev.off()