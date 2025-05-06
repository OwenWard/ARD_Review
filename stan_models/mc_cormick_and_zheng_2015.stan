// May 6th 2025
// Stan Model for the Overdispersed Binomial
// Model of McCormick and Zheng et al 2015
// modified to take varying priors for each comp of beta



// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=1> p; // the latent dimension
  array[N, K] int y;
  int<lower=0, upper=K> n_known;
  array[n_known] int<lower=1, upper=K> idx;
  real<lower=0, upper=1> known_prev;
}


parameters {
  real mu_alpha;
  real mu_beta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector[N] alpha;
  vector[K] beta;
  vector[K] eta;
  real xi;
  array[N] unit_vector[p] z;
  array[K] unit_vector[p] nu;
}


transformed parameters {
  vector[K] scaled_beta;
  vector[N] scaled_alpha;
  real C;
  C = log(sum(exp(beta[idx])/known_prev) );  // sum over known indices
  // TO DO, modify for different notation here (b vs beta)
  scaled_alpha = alpha + C;
  scaled_beta = beta - C;
}

// The model to be estimated..
model {
  xi ~ gamma(2, 1); // setting this for now
  eta ~ gamma(2, 1);
  // sigma_beta ~ normal(0, 5);
  // sigma_alpha ~ normal(0, 5);
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta); 
  
  for(n in 1:N){
    for(k in 1:K){
      real gamma_ik = ...//inv_omega[k]/(1 - inv_omega[k]) ; //the hypersphere part
      real scaled_deg = exp(scaled_alpha[n]);
      y[n,k] ~ poisson(scaled_deg * scaled_beta[k] * gamma_ik);
    }
  }
}



generated quantities {
  array[N] int y_sum;
  array[N, K] int y_sim;
  array[N,K] real log_lik;
  for (n in 1:N) {
    for (k in 1:K) {
      real gamma = inv_omega[k]/(1 - inv_omega[k]) ; //the beta par of neg_bin
      real xi_i_k = gamma * exp(scaled_alpha[n] + scaled_beta[k])  ;
      y_sim[n, k] = neg_binomial_rng(xi_i_k, gamma) ;
      log_lik[n, k] = neg_binomial_lpmf(y[n, k] | xi_i_k, gamma);
    }
    y_sum[n] = sum( y_sim[n] );
  }

}
