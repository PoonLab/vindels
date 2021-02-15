data {
  int<lower=0> N;                     // Number of columns    
  int<lower=0> counts[N];         // Count data
  real<lower=0> branches[N];     // Time data
}


parameters {
  real<lower=-30, upper=30> rate;
}

model {
  rate ~ uniform(-30,30);
  for (i in 1:N){
    counts[i] ~ poisson_log_lpmf(rate + log(branches[i]));
  }
}

