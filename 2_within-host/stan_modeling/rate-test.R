
data {                      // Number of patients
  int<lower=0> ntree;                         // Number of columns    
  int row;                            // Vector describing the number of rows belonging to each patient
  int<lower=0> counts[row,ntree];         // Count data
  matrix<lower=0>[row,ntree] lengths;     // Time data
}

parameters {
  real<lower=0> sub_rate;
  real<lower=0> tre_rates[ntree];
}

model {
  
  tre_rates ~ normal(sub_rate, 0.05)
  
  for (j in 1:ntree){      //  1 - #Columns
      counts[j]  ~ poisson_log(tre_rates[j] * lengths[j]);  
  }
  
  print(target())
}
