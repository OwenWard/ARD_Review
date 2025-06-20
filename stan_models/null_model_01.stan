// October 15th 2024
// Stan Model for the Simplest Null Model


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> K;
  array[N, K] int y;
}


parameters {
  real log_d; // the common population degree
  vector[K] beta;
}


transformed parameters {
  // vector<lower=0, upper =1>[K] b = inv_logit(beta);  // element-wise
  vector[K] b = exp(beta);
}


model {
  
  log_d ~ normal(0, 25);
  beta ~ normal(0, 5); 
  real exp_log_d = exp(log_d);
  for (n in 1:N) {
    y[n] ~ poisson(exp_log_d * b);
  }
}


// skip these for now

generated quantities {
  array[N] int y_sum;
  array[N, K] int y_sim;
  vector[N] log_lik;
  real exp_log_d = exp(log_d);
  for (n in 1:N) {
    y_sim[n] = poisson_rng(exp_log_d * b);
    y_sum[n] = sum( y_sim[n] );
    log_lik[n] = poisson_log_lpmf(y[n] | exp_log_d * b);
  }

}
