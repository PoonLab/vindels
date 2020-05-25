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
  int<lower=0> ntree;
  int<lower=0> nbranch;
  matrix[nbranch,ntree] counts;
  matrix[nbranch,ntree] offset;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real pat_rate;
  real<lower=0> pat_sd;

  vector[ntree] tre_rate;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  tre_rate ~ normal(pat_rate, pat_sd);
  
  for (i in 1:ntree){
    counts[i] ~ poisson_log(tre_rate[i] + offset[i]); //+ offset[i]); // within patient distribution 
  }
}
