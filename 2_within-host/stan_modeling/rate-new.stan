data {
  int<lower=0> npat;                        // Number of patients
  int<lower=0> ntree;                         // Number of columns    
  int<lower=0> mx;
  int <lower=0> sizes[npat];                            // Vector describing the number of rows belonging to each patient
  int<lower=0> counts[ntree, mx, npat];         // Count data
  real<lower=0> branches[ntree, mx, npat] ;     // Time data
}

parameters {
  real<lower=0.000001, upper=0.1> sub_rate;
  real<lower=0.000001, upper=1> sub_sd;
  real<lower=0.000001, upper=0.1> pat_rates[npat];
  real<lower=0.000001,upper=1> pat_sd;
  real<lower=0.000001, upper=0.1> tre_rates[ntree] ;
}
// transformed parameters {
//   real<lower=0> lambda;
//   lambda=tre_rates * branches
// }


model {
  int pos;
  int lim;
  pat_rates ~ normal(sub_rate, sub_sd);  // 0.2 is arbitrary
  
  pos = 1;
  for (i in 1:npat){
    tre_rates ~ normal(pat_rates[i], pat_sd);  // 0.2 is arbitrary

    counts[1:ntree, 1:sizes[i], i] ~ poisson(tre_rates * branches[1:ntree, 1:sizes[i], i]  );
  }

}

