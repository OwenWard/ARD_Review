// July 5th 2024
// Stan Model for the Overdispersed Binomial
// Model of Zheng et al 2006



// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> K;
  real mu_beta; // secify these as known
  real<lower=0> sigma_beta; // specify these as known
  array[N, K] int y;
}


parameters {
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[N] alpha;
  vector[K] beta;
  vector<lower=0, upper=1>[K] inv_omega;
}

// The model to be estimated..
model {
  
  mu_alpha ~ normal(0, 5);
  sigma_alpha ~ normal(0, 5);

  
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta); 
  
  for(n in 1:N){
    for(k in 1:K){
      real gamma = inv_omega[k]/(1 - inv_omega[k]) ; //the beta par of neg_bin
      real xi_i_k = gamma * exp(alpha[n] + beta[k])  ;
      y[n,k] ~ neg_binomial(xi_i_k, gamma);
    }
  }
}



generated quantities {
  array[N] int out_degree;
  array[N, K] int y_sim;
  for (n in 1:N) {
    for (k in 1:K) {
      real gamma = inv_omega[k]/(1 - inv_omega[k]) ; //the beta par of neg_bin
      real xi_i_k = gamma * exp(alpha[n] + beta[k])  ;
      y_sim[n, k] = neg_binomial_rng(xi_i_k, gamma) ;
    }
    out_degree[n] = sum( y_sim[n] );
  }
  
}
