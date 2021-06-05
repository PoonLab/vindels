data {
  int<lower=0> npat;                        // Number of patients
  int<lower=0> ntree;                         // Number of columns    
  int <lower=0> sizes[npat];                            // Vector describing the number of rows belonging to each patient
  int<lower=0> counts[sum(sizes),ntree];         // Count data
  matrix[sum(sizes),ntree] l_branches;     // Time data
}

parameters {
  real<lower=-15,upper=0> sub_rate;
  real<lower=0, upper=5> sub_sd;

}

transformed parameters {
  real<lower=-15,upper=0> pat_rates[npat];
  real<lower=0,upper=5> pat_sd;
  real<lower=-15, upper=0> tre_rates[npat, ntree];
  
  vector[npat] alpha; // The patient rates
  real theta[npat, ntree];  // The tree rates
 
  alpha = sub_rate + sub_sd * pat_rate_unitized;
  
  for (i in 1:npat){
    for (j in 1:ntree)
    {
      theta[i,j] = alpha[i] + sigma * tree_rate_unitized[i,j];
    }
  }

}


model {
  int pos;
  int lim;
  pat_rate_unitized ~ normal(0, 1);
  tree_rate_unitized ~ normal(0, 1);
  
  pat_rates ~ normal_lpdf(sub_rate, sub_sd);  // 0.2 is arbitrary
  
  pos = 1;           
  for (i in 1:npat){
    tre_rates[i,] ~ normal_lpdf(pat_rates[i], pat_sd);  // 0.2 is arbitrary
    lim = pos+sizes[i]-1;
    
    for (j in 1:ntree){      
      counts[pos:lim, j] ~ poisson_log_lpmf(tre_rates[i,j] + l_branches[pos:lim,j]); 
    }
    pos = pos + sizes[i];
  }

}

