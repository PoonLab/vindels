# ----- SLIPPAGE MODEL -----
csv <- read.csv("~/PycharmProjects/hiv-withinhost/15_modeling/slip-15-fix-perfect.csv", stringsAsFactors = F, comment.char="#")
csv2 <- read.csv("~/PycharmProjects/hiv-withinhost/15_modeling/slip-16-fix-perfect.csv", stringsAsFactors = F, comment.char="#")

csv 

# csv <- as.data.frame(sapply(1:ncol(csv), function(x){
#   sapply(1:nrow(csv), function(y){
#     as.numeric(strsplit(csv[y,x]," ")[[1]][1])
#   })
# }))
burnin <- ceiling(0.1*nrow(csv))
tru <- c(0.00052, 0.80, 0.00001, 0.09)
#colnames(csv) <- c('p.enter', 'p.stay', 'rate', 'fix','likelihood', 'slip.changed', 'accept', 'time')
# colnames(csv) <- c('p.enter', 'slope','int', 'rate', 'slip.changed', 'accept', 'time')
#png(file="~/vindels/Figures/within-host/finalized/slippage-trace2.png", width=600, height=800)
par(mar=c(2.5,5,3,1), mfrow=c(2,2))
plot(csv[-(1:burnin)
         ,1], type = "l", xlab="MCMC Steps (x10)" , ylab="P(Enter)",
     main = "Chain values of Enter", cex.axis=1.3, cex.lab=1.4, cex.main=1.7)#, ylim=c(0.0005, 0.0012))
abline(h=tru[1],col='red',lwd=2)
plot(csv[-(1:burnin)
         ,2], type = "l", xlab="MCMC Steps (x10)" , ylab="P(Stay)",
     main = "Chain values of Stay",  cex.axis=1.3, cex.lab=1.4, cex.main=1.7)#, ylim=c(0.62, 0.77))
abline(h=tru[2],col='red',lwd=2)
plot(csv[-(1:burnin)
         ,3], type = "l", xlab="MCMC Steps (x10)" , ylab="Rate",
     main = "Chain values of Rate",  cex.axis=1.3, cex.lab=1.4, cex.main=1.7)#, ylim=c(0.00001,0.000018))
abline(h=tru[3],col='red',lwd=2)
plot(csv[-(1:burnin)
         ,4], type = "l", xlab="MCMC Steps (x10)" , ylab="Fix",
     main = "Chain values of Fix",  cex.axis=1.3, cex.lab=1.4, cex.main=1.)#,ylim=c(0.1,0.25))
abline(h=tru[4],col='red',lwd=2)
plot(csv[-(1:burnin)
         ,'fix-sd'], type = "l", xlab="MCMC Steps (x10)" , ylab="Fix - SD",
     main = "Chain values of Fix SD",  cex.axis=1.3, cex.lab=1.4, cex.main=1.7)#, ylim=c(0.0001,0.00015))

par(mar=c(6,5.5,6,1),mfrow=c(1,3))
hist(csv[-(1:burnin),1], 
     col="dodgerblue", 
     breaks=15, 
     main="Posterior of 'Enter'", 
     xlab="P(Enter)", cex.axis=1.3, cex.lab=1.5, cex.main=1.8)
abline(v=tru[1],col='red',lwd=2)
hist(csv[-(1:burnin),2], 
     col="dodgerblue", 
     breaks=15, 
     main="Posterior of 'Stay'", 
     xlab="P(Stay)", cex.axis=1.3, cex.lab=1.5, cex.main=1.8)
abline(v=tru[2],col='red',lwd=2)
hist(csv[-(1:burnin),3], 
     col="dodgerblue", 
     breaks=15, 
     main="Posterior of 'Rate'", 
     xlab="Rate", cex.axis=1.3, cex.lab=1.5, cex.main=1.8)
abline(v=tru[3],col='red',lwd=2)
med1 <- median(csv$p.enter)
med2 <- median(csv$p.stay)
med3 <- median(csv$rate)
#dev.off()