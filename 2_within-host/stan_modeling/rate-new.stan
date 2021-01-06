data {
  int<lower=0> N;                     // Number of columns    
  int<lower=0> counts[N];         // Count data
  real<lower=0> branches[N];     // Time data
}


parameters {
  real<lower=-20, upper=5> rate;
}

model {
  rate ~ normal(0,20);
  for (i in 1:N){
    counts[i] ~ poisson(exp(rate) * branches[i]);
  }
}

