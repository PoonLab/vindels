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
  int<lower=0> ntree;                  // Number of columns           // Number of branches 
  int sizes[npat];                   // Number of rows
  int<lower=0> counts[sum(sizes),ntree];         // Count data
  matrix<lower=0>[sum(sizes),ntree] lengths;     // Time data
}

parameters {
  real<lower=0> sub_rate;
  //real<lower=0> sub_sd;
  real<lower=0> pat_rates[npat];
  //<lower=0> pat_sd;
  real<lower=0> tre_rates[ntree];
}

model {
  int pos;
  
  pat_rates ~ normal(sub_rate, 0.03);
  
  pos = 1;
  for (i in 1:npat){
    tre_rates ~ normal(pat_rates[i], 0.06);
    for (j in 1:ntree){
      //print(pos);
      //print(pos+sizes[i]-1);
      print(counts[pos:(pos+sizes[i]-1), j]) // ~ poisson(exp(tre_rates[j] * lengths[pos:(pos+sizes[i]-1), j]));  
    }
    pos = pos + sizes[i];
  }
}

