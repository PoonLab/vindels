# stan model for indel rates 
library(rstan)

num.pat <- 100
num.data <- 200
rstan_options(auto_write=T)
options(mc.cores = parallel::detectCores()-2)

t.rate <- 10
sub.rate <- rlnorm(num.pat, meanlog=log(t.rate), sdlog=0.3)

hist(sub.rate, breaks=5, col="grey")

mat <- matrix(nrow=num.data, ncol=num.pat)

for (i in 1:num.pat){
  mat[,i] <- rnorm(num.data, mean=sub.rate[i], sd=1.2)
}


data.stan <- list(npat = num.pat,
                  ndata = num.data,
                  mat = mat)

stan.fit <- stan("~/vindels/2_within-host/rates.stan",
                 data= data.stan, 
                 chains=1,
                 iter=1000000)
rstan::summary(stan.fit)

traceplot(stan.fit, alpha=0.7)
w <- extract(stan.fit)

plot(density(w[[2]]))
abline(v=mm)


plot(density(w$lambda))
