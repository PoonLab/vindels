data {
  int<lower=0> N;
  real<lower=0> input[N];
  //real<lower=0> times[N];
}

parameters {
  //real<lower=-20> rate;
  real<lower=-20> mu;
  real<lower=0> sigma;
}

// transformed parameters {
//   real lambda[N];
//   for (i in 1:N){
//     lambda[i] = log(times[i]) + rate
//   }
// }

model {
  // for (i in 1:N){
  //   input[i] ~ poisson_log_lpmf(rate + log(times[i]));
  // }
  
  input ~ lognormal_lpdf(mu, sigma);  
}

