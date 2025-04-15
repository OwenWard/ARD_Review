// October 15th 2024
// Stan Model for the Simplest Null Model


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> K;
  array[N, K] int y;
}


parameters {
  vector[N] log_d; // the individual degree
  vector<lower=0, upper=1>[K] beta;
}

model {
  
  log_d ~ normal(0, 25);
  beta ~ beta(1, 1); 
  
  for (n in 1:N) {
    real exp_log_d = exp(log_d[n]);
    y[n] ~ poisson(exp_log_d * beta);
  }
}


// skip these for now

generated quantities {
  array[N] int y_sum;
  array[N, K] int y_sim;
  vector[N] log_lik;
  for (n in 1:N) {
    real exp_log_d = exp(log_d[n]);
    y_sim[n] = poisson_rng(exp_log_d * beta);
    y_sum[n] = sum( y_sim[n] );
    log_lik[n] = poisson_log_lpmf(y[n] | exp_log_d * beta);
  }

}
