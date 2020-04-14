# stan model for indel rates 
library(rstan)


# will receive data in from 9_5 indel rates 
num.pat <- 19
num.data <- 200
rstan_options(auto_write=T)
options(mc.cores = parallel::detectCores()-2)


# ---- Simulated Data ----
t.rate <- 2.5
sub.rate <- rlnorm(num.pat, meanlog=log(t.rate), sdlog=0.3)

hist(sub.rate, breaks=5, col="grey")

mat <- matrix(nrow=num.data, ncol=num.pat)

for (i in 1:num.pat){
  mat[,i] <- rnorm(num.data, mean=sub.rate[i], sd=0.5)
}

# ----- Real Data ----
mat <- matrix(nrow=num.data, ncol=num.pat)
for (x in 1:length(V1)){
  mat[,x] <- V1[[x]]
}

# ----- Stan Model ----
data.stan <- list(npat = num.pat,
                  ndata = num.data,
                  mat = mat)

stan.fit <- stan("~/vindels/2_within-host/rates.stan",
                 data= data.stan, 
                 chains=1,
                 iter=100000)
pdf()
traceplot()
dev.off()
rstan::summary(stan.fit)

traceplot(stan.fit, alpha=0.7)
x <- extract(stan.fit)

par(mfrow=c(1,2))
plot(density(x[[1]]))
abline(v=log(t.rate), col="red")
plot(density(x[[2]]))
#abline(v=median(mat), col="red")


plot(density(w$lambda))
