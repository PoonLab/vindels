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

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> npat;
  int<lower=0> ndata;
  vector[npat] mat[ndata];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real <lower=0.01, upper=100> sub_rate;
  real<lower=0.01, upper=100> pat_rate;
  real <lower=0,upper=10> sub_sd;
  real <lower=0,upper=10> pat_sd;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  pat_rate ~ lognormal(sub_rate, sub_sd);
  
  for (pat in 1:npat){
    mat[pat] ~ normal(pat_rate, pat_sd);
  }
}
