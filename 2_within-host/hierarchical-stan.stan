//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//


data {
  int<lower=0> n;
  int<lower=0> np;
  vector[np] M[n];  // https://mc-stan.org/docs/2_22/stan-users-guide/basic-motivation.html
}

parameters {
  //real<lower=0.01, upper = 100> lambda;
  real<lower=0.01, upper=1> mlog;
  real<lower=0.01, upper = 100> lambda;
}

model {
  
  lambda ~ lognormal(mlog, 0.3);
  
  for (p in 1:np){
    // M[p] ~ exponential(lambda);
    M[p] ~ exponential(lambda);
  }
}

