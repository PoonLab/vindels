data {
  int<lower=0> npat;                        // Number of patients
  int<lower=0> ntree;                         // Number of columns    
  int <lower=0> sizes[npat];                            // Vector describing the number of rows belonging to each patient
  int<lower=0> counts[sum(sizes),ntree];         // Count data
  matrix<lower=0>[sum(sizes),ntree] lengths;     // Time data
}

parameters {
  real<lower=0.000001, upper=0.1> sub_rate;
  real<lower=0.00001, upper=1> sub_sd;
  real<lower=0.00001, upper=0.1> pat_rates[npat];
  real<lower=0.00001, upper=1> pat_sd;
  real<lower=0.00001, upper=0.1> tre_rates[ntree];
}

model {
  int pos;
  int lim;
  
  pat_rates ~ normal(sub_rate, sub_sd);  // 0.2 is arbitrary
  
  pos = 1;
  for (i in 1:npat){
    tre_rates ~ normal(pat_rates[i], pat_sd);  // 0.2 is arbitrary
    
    lim = pos+sizes[i]-1;
    for (j in 1:ntree){      
      counts[pos:lim, j] ~ poisson(tre_rates[j] * lengths[pos:lim, j]); 
      print(target())
    }
    pos = pos + sizes[i];
  }
}

