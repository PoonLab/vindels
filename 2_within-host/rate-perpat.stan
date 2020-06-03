data {
  int<lower=0> npat;                        // Number of patients
  int<lower=0> ntree;                         // Number of columns    
  int sizes[npat];                            // Vector describing the number of rows belonging to each patient
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
  int lim;
  
  pat_rates ~ normal(sub_rate, 0.03);
  
  pos = 1;
  for (i in 1:npat){
    tre_rates ~ normal(pat_rates[i], 0.06);
    
    lim = pos+sizes[i]-1;
    for (j in 1:ntree){      //  1 - #Columns
      counts[pos:lim, j]  ~ poisson(exp(tre_rates[j] * lengths[pos:lim, j]));  
      // Alternate format:  segment(counts[j], pos, sizes[i]) ~ poisson(exp(tre_rates[j] * segment(lengths, pos, sizes[i]))); 
    }
    pos = pos + sizes[i];
  }
}

