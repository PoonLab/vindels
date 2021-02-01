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