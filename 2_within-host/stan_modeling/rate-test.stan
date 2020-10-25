
data {                      // Number of patients
  int<lower=0> ntree;                         // Number of columns    
  int<lower=0> nrow;                            // Vector describing the number of rows belonging to each patient
  int<lower=0> counts[nrow,ntree];         // Count data
  matrix<lower=0>[nrow,ntree] lengths;     // Time data
}

parameters {
  real<lower=0> sub_rate;
  real<lower=0,upper=100> tre_rates[ntree];
}



model {
  
  tre_rates ~ normal(sub_rate, 0.05);
  //print(target())
  for (j in 1:ntree){      //  1 - #Columns
    //print(j)
    //print(lengths[j] * tre_rates[j])
    counts[j] ~  poisson_log(tre_rates[j] * lengths[j]);  
  }
  
}
