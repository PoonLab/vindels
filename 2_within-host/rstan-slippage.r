# rstan script for modeling the 


insertions <- "
data {  
  int N;                   // number of possible insertion locations 
  real<lower=0> counts[N]; // slippage events
  

}
parameters {
  real<lower=0, upper=1> slip;
  real<lower=0, upper=1> mut;
}

model { 
   counts ~ neg_binomial(1, slip/(1-slip))
} 
"

slip_data <- list(N=c(100), counts=slips)
fit <- stan(model_code="slip", data=slip_data, iter=10000, warmup=100, chains=8)

affll <- function(p.enter, p.stay, slips){
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[slips != 0] - 1)
  
  print(x)
  print(y)
  print(z)
  
  (1 - p.enter)^x * p.enter^y * (1 - p.stay)^y * p.stay^z
}

obj.f <- function(p.enter, p.stay){
  -affine(p.enter, p.stay)
}
mle3 <- bbmle::mle2(objf, start=list(p.enter=1, p.stay=1), method="Brent" , lower=1e-12, upper=1)


mcmc <- function(tipseq, ancseq, slip, ){
  
}
  

slip_start <- runif(1)
  
  
runif