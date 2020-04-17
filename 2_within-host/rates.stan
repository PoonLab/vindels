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
  //real sub_rate;
  real pat_rate;
  real<lower=0, upper=5> pat_sd;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  pat_rate ~ uniform(0, 7);  // --> after running, posterior is normally distributed; 
  pat_sd ~ uniform(0,10);           // this shows the prior distribution and the
  //pat_rate ~ lognormal(sub_rate, sub_sd);
  
  for (i in 1:npat){
    mat[i] ~ normal(pat_rate, pat_sd); // within patient distribution 
  }
}

