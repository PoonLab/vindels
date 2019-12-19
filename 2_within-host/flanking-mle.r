# MAXIMUM LIKELIHOOD ESTIMATION 
# ------------------------------------------------------

# POISSON DISTRIBUTION


flanking <- read.csv("~/PycharmProjects/hiv-withinhost/14_flanking/flanking.csv",row.names=1, stringsAsFactors = F)
all <- read.csv("~/PycharmProjects/hiv-withinhost/14_flanking/flanking-all.csv",row.names=1, stringsAsFactors = F)

pll <- function(rate, count, len){
  lam <- rate * len 
  res <- -lam + count*log(lam) 
  res <- res[!is.na(res)]
  sum(res)
}
obj.f <- function(rate) -pll(rate, all$new.count, all$Vlength)
mle.result <- bbmle::mle2(obj.f, start=list(rate=1), method = "Brent", lower=1e-12, upper = 1)


# BINOMIAL DISTRIBUTION 
binomll <- function(prob, count, len){
  N <- len
  k <- count
  p <- prob
  chs <- factorial(N) / (factorial(k) * factorial(N - k))
  sum(log(chs) +  k*log(p) +  (N - k)*log(1-p))
}
obj.f2 <- function(prob) -binomll(prob, all$new.count, all$Vlength)
mle.result2 <- bbmle::mle2(obj.f2, start=list(prob=1), method = "Brent", lower = 1e-12, upper=1)


# GEOMETRIC LIKELIHOOD FUNCTION 

geomll <- function(forward, count, N){
  N * log(forward) + (sum(count)-N) * log(1-forward)
}
slips.whole <- c(all$count.flanking, rep(0,sum(nchar(all$Vseq))-nrow(all)))
slips.nt <- c(all$new.count, rep(0,sum(nchar(all$Vseq))-nrow(all)))
slips.nt2 <- slips.nt + 1
objf3 <- function(forward) -geomll(forward, slips.nt2, length(slips.nt2))
mle3 <- bbmle::mle2(objf3, start=list(forward=1), method="Brent" , lower=1e-8, upper=1)


count.subs <- function()
  
  
  # CUSTOM GEOMETRIC LIKELIHOOD FUNCTION 
  custom1 <- function(slip, mut, len, subs, N){
    N * log(slip) + sum(count) * log(1-slip) + log(choose(count, subs)) + subs * log(mut) + (counts - subs) * log(1-mut)
  }

obj4 <- function(slip, mut) -custom1(slip, mut, slips.nt, length(slips.nt))
mle4 <- bbmle::mle2(obj4, start=list(slip=1, mut=1), method="L-BFGS-B",lower=1e-12, upper=1)

x <- runif(5000,min=1e-7, max=1)
y <- unname(sapply(x, objf3))
plot(x=x, y=y)

obs <- all$new.count
lens <- c()
dists <- list()
for (i in 1:25){
  sim <- unname(sapply(all$Vseq, simulateData))
  diff <- nchar(sim) - nchar(all$Vseq)
  simdiff <- diff[diff!=0]
  lens[i] <- length(simdiff)
  dists[[i]] <- simdiff
}
