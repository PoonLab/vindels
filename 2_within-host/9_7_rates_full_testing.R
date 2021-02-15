# stan model testing 


num.pat <- 13
num.data <- 200

sub.sd <- 2
sub.rate <- 2.5

pat.sd <- 0.8
pat.rate <- rnorm(num.pat, mean=sub.rate, sd=sub.sd)


mat <- matrix(ncol=num.pat, nrow=num.data)
for (n in 1:num.pat){
  mat[,n] <- rnorm(num.data, pat.rate[n], pat.sd)
}


data.stan <- list(npat = num.pat,
                  ndata = num.data,
                  mat = mat)



stan.fit <- stan("~/vindels/2_within-host/rates.stan",
                 data= data.stan, 
                 chains=1,
                 iter=1000000,
                 control = list(adapt_delta = 0.99))

# RLNORM TESTING
N = 10000
data.stan <- list(N=N,input = rlnorm(N,meanlog=log(0.005),sdlog=1))

stan.fit <- stan("~/vindels/2_within-host/stan_modeling/lnorm.stan",
                 data= data.stan, 
                 chains=1,
                 iter=3000,
                 control = list(adapt_delta = 0.99))

# POISSON TESTING
N = 10000
rate = 0.005
times <- rlnorm(N,4.5,0.3) 
data.stan <- list(N=N,input = rpois(N,lambda=rate*times), times=times)

stan.fit <- stan("~/vindels/2_within-host/stan_modeling/lnorm.stan",
                 data= data.stan, 
                 chains=1,
                 iter=3000,
                 control = list(adapt_delta = 0.99))