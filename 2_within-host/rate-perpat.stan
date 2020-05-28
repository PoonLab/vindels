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
  int<lower=0> ntree;                // Number of columns
  int sizes[npat];                   // Number of rows
  int counts[sum(sizes),ntree];         // Count data
  matrix[sum(sizes),ntree] lengths;     // Time data
}

parameters {
  real<lower=0> sub_rate;
  real<lower=0> sub_sd;
  real<lower=0> pat_rates[npat];
  real<lower=0> pat_sd;
  real<lower=0> tre_rates[ntree];
}

model {
  int pos;
  
  pat_rates ~ normal(sub_rate, sub_sd);
  
  pos = 1;
  for (i in 1:npat){
    tre_rates ~ normal(pat_rates[i], pat_sd);
    
    for (j in 1:ntree){
      segment(counts[j], pos, sizes[i]) ~ poisson(exp(tre_rates[j] * segment(lengths[j], pos, sizes[i])));  
    }
    pos = pos + sizes[i];
    
  }
  

}

