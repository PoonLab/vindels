data {
  int<lower=0> npat;                        // Number of patients
  int<lower=0> ntree;                         // Number of columns    
  int <lower=0> sizes[npat];                            // Vector describing the number of rows belonging to each patient
  int<lower=0> counts[sum(sizes),ntree];         // Count data
  matrix<lower=0>[sum(sizes),ntree] branches;     // Time data
}

parameters {
  real<lower=-20,upper=10> sub_rate;
  real<lower=0, upper=5> sub_sd;
  real<lower=-20,upper=10> pat_rates[npat];
  real<lower=0,upper=5> pat_sd;
  real<lower=-20, upper=0> tre_rates[ntree, npat];
}

model {
  int pos;
  int lim;
  pat_rates ~ normal_lpdf(sub_rate, sub_sd);  // 0.2 is arbitrary
  
  pos = 1;           
  for (i in 1:npat){
    tre_rates[,i] ~ normal_lpdf(pat_rates[i], pat_sd);  // 0.2 is arbitrary
    lim = pos+sizes[i]-1;
    
    for (j in 1:ntree){      
      counts[pos:lim, j] ~ poisson_log_lpmf(tre_rates[j,i] + log(branches[pos:lim, j])); 
    }
    pos = pos + sizes[i];
  }

}

