# ----- SLIPPAGE MODEL -----
csv <- read.csv("~/PycharmProjects/hiv-withinhost/15_modeling/slip-21-full-stacks.csv", stringsAsFactors = F, comment.char="#")
csv2 <- read.csv("~/PycharmProjects/hiv-withinhost/15_modeling/slip-18-fix-perfect.csv", stringsAsFactors = F, comment.char="#")

# csv <- as.data.frame(sapply(1:ncol(csv), function(x){
#   sapply(1:nrow(csv), function(y){
#     as.numeric(strsplit(csv[y,x]," ")[[1]][1])
#   })
# }))
data <- csv
burnin <- ceiling(0.1*nrow(data))
tru <- c(0.00052, 0.85, 0.00001, 0.09)
#colnames(data) <- c('p.enter', 'p.stay', 'rate', 'fix','likelihood', 'slip.changed', 'accept', 'time')
# colnames(data) <- c('p.enter', 'slope','int', 'rate', 'slip.changed', 'accept', 'time')
#png(file="~/vindels/Figures/within-host/finalized/slippage-trace2.png", width=600, height=800)
par(mar=c(2.5,5,3,1), mfrow=c(2,2))
plot(data[-(1:burnin)
         ,1], type = "l", xlab="MCMC Steps (x10)" , ylab="P(Enter)",
     main = "Chain values of Enter", cex.axis=1.3, cex.lab=1.4, cex.main=1.7)#, ylim=c(0.0005, 0.0012))
abline(h=tru[1],col='red',lwd=2)
plot(data[-(1:burnin)
         ,2], type = "l", xlab="MCMC Steps (x10)" , ylab="P(Stay)",
     main = "Chain values of Stay",  cex.axis=1.3, cex.lab=1.4, cex.main=1.7)#, ylim=c(0.62, 0.77))
abline(h=tru[2],col='red',lwd=2)
plot(data[-(1:burnin)
         ,3], type = "l", xlab="MCMC Steps (x10)" , ylab="Rate",
     main = "Chain values of Rate",  cex.axis=1.3, cex.lab=1.4, cex.main=1.7, ylim=c(0.00001,0.000013))
abline(h=tru[3],col='red',lwd=2)
plot(data[-(1:burnin)
         ,4], type = "l", xlab="MCMC Steps (x10)" , ylab="Fix",
     main = "Chain values of Fix",  cex.axis=1.3, cex.lab=1.4, cex.main=1.)#,ylim=c(0.1,0.25))
abline(h=tru[4],col='red',lwd=2)


plot(csv[-(1:burnin)
         ,'fix-sd'], type = "l", xlab="MCMC Steps (x10)" , ylab="Fix - SD",
     main = "Chain values of Fix SD",  cex.axis=1.3, cex.lab=1.4, cex.main=1.7)#, ylim=c(0.0001,0.00015))


# ---- Pretty Blue Histograms ----
main <- 2
lab <- 1.8
ax <- 1.5

par(mar=c(6,5,5,0.5),las=1,mfrow=c(1,4))
hist(data[-(1:burnin),1], 
     freq=F,
     col="dodgerblue", 
     breaks=15, 
     main="Posterior of 'Enter'", 
     ylab="",
     yaxt="n",
     xlab="P(Enter)", cex.axis=ax, cex.lab=lab, cex.main=main)
abline(v=tru[1],col='black',lwd=2,lty=2)
lines(xy.coords(x=c(0,1), y=c(1,1)),col="red",lwd=2.5)
#title(ylab="Density", line=3, cex.lab=lab)
axis(2,
     labels=c("0","5e3", "1e4", "1.5e4"),
     at=c(0,5000,10000,15000),
     cex.axis=ax)
hist(rep(csv[-(1:burnin),2],3), 
     freq=F,
     col="dodgerblue", 
     breaks=15, 
     main="Posterior of 'Stay'", 
     xlab="P(Stay)", 
     ylab="",
     cex.axis=ax, cex.lab=lab, cex.main=main)
abline(v=tru[2],col='black',lwd=2, lty=2)
lines(xy.coords(x=c(0.4,0.9), y=c(2,2)),col="red",lwd=2.5)
options(scipen=100)
hist(csv[-(1:burnin),3], 
     freq=F,
     col="dodgerblue", 
     breaks=15, 
     ylab="",
     yaxt="n",
     main="Posterior of 'Rate'", 
     xlim=c(0.00001,0.000013),
     xlab="Rate", cex.axis=ax, cex.lab=lab, cex.main=main)
abline(v=tru[3],col='black',lwd=2, lty=2)
lines(xy.coords(x=c(1e-3,1e-7), y=c(1000,1000)),col="red",lwd=2.5)
axis(2, 
     labels=c("0","4e5", "8e5", "1.2e6"), 
     at=c(0,4e5,8e5,1.2e6),
     cex.axis=ax)
hist(csv[-(1:burnin),4], 
     freq=F,
     col="dodgerblue", 
     breaks=15, 
     ylab="",
     main="Posterior of 'Fixation'", 
     xlab="P(Fixation)", cex.axis=ax, cex.lab=lab, cex.main=main)
abline(v=tru[4],col='black',lwd=2, lty=2)
lines(density(rbeta(10000,3,20)),col="red",lwd=2.5)
med1 <- median(csv$p.enter)
med2 <- median(csv$p.stay)
med3 <- median(csv$rate)


# single parameter
par(mar=c(6,7,5,2),las=1)#,mfrow=c(2,2))
hist(chain[-(1:burnin),1], 
     freq=F,
     col="dodgerblue", 
     breaks=15, 
     main="Posterior of 'Slip'", 
     ylab="",
     #yaxt="n",
     xlab="P(Slip)", cex.axis=ax, cex.lab=lab, cex.main=main)
abline(v=0.001,col='black',lwd=2,lty=2)
lines(xy.coords(x=c(0,1), y=c(1,1)),col="red",lwd=2.5)
title(ylab="Density", line=4.5, cex.lab=lab)




#dev.off()