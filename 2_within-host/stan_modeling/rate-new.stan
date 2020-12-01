data {
  int<lower=0> N;                     // Number of columns    
  int<lower=0> counts[N];         // Count data
  real<lower=0> branches[N];     // Time data
}


parameters {
  real<lower=0, upper=1> rate;
}

model {
  rate ~ uniform(0,1);
  for (i in 1:N){
    counts[i] ~ poisson( exp(rate + log(branches[i])));
  }
  
}

