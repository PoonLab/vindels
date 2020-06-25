# random figures to generate for thesis 
wid <- 4
bw <- 0.2
data1 <- rnorm(10000, mean=5, sd=2)
data2 <- rnorm(10000, mean=5.8, sd=0.7)
data3 <- rnorm(10000, mean=6.4, sd=0.3)
data4 <- rnorm(10000, mean=6.6, sd=0.07)
plot(density(data1, bw=0.8), xlim=c(2,8), ylim=c(0,2), lwd=wid, col="black")
lines(density(data2, bw=bw), lwd=wid, col="red3")
lines(density(data3, bw=bw), lwd=wid, col="red2")
lines(density(data4, bw=bw), lwd=wid, col="red1")
