// October 15th 2024
// Stan Model for the Simplest Null Model



// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> K;
  array[N, K] int y;
  int<lower=0, upper=K> n_known;
  array[n_known] int<lower=1, upper=K> idx;
  real<lower=0, upper=1> known_prev;
  //int<lower=0> N_totalpop;
}


parameters {
  real log_d; // the common population degree
  vector[K] beta;
}


transformed parameters {
  vector[K] scaled_beta;
  real scaled_log_d;
  real C;
  C = log(sum(exp(beta[idx])/known_prev) );  // sum over known indices
  scaled_log_d = log_d + C;
  scaled_beta = beta - C;
  vector[K] b = exp(scaled_beta);
}


model {
  
  log_d ~ normal(0, 25);
  beta ~ normal(0, 5); 
  real exp_log_d = exp(scaled_log_d);
  for (n in 1:N) {
    y[n] ~ poisson(exp_log_d .* b);
  }
}



generated quantities {
  //vector[K] size_ests = b*N_totalpop;
  
  array[N] int y_sum;
  array[N, K] int y_sim;
  array[N, K] real log_lik;
  real curr_log_d = scaled_log_d;
  for (n in 1:N) {
    y_sim[n] = poisson_log_rng(curr_log_d + scaled_beta);
    y_sum[n] = sum( y_sim[n] );
    log_lik[n] = poisson_log_lpmf(y[n] | curr_log_d + scaled_beta);
  }

}
