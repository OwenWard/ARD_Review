// October 9th 2024
// Stan Model for the Overdispersed Binomial
// Model of Zheng et al 2006
// modified to take varying priors for each comp of beta



// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> K;
  // vector[K] mu_beta; // secify these as known
  // vector<lower=0>[K] sigma_beta; // specify these as known
  array[N, K] int y;
  int<lower=0, upper=K> n_known;
  array[n_known] int<lower=1, upper=K> idx;
  real<lower=0, upper=1> known_prev;
  array[N, K] int<lower=0,upper=1> obs_mask; 
}


parameters {
  real mu_alpha;
  real mu_beta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector[N] alpha;
  vector[K] beta;
  vector<lower=0, upper=1>[K] inv_omega;
}


transformed parameters {
  vector[K] scaled_beta;
  vector[N] scaled_alpha;
  real C;
  C = log(sum(exp(beta[idx])/known_prev) );  // sum over known indices
  scaled_alpha = alpha + C;
  scaled_beta = beta - C;
}

// The model to be estimated..
model {
  
  // mu_alpha ~ normal(0, 25);
  sigma_alpha ~ normal(0, 5);
  alpha ~ normal(0, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta); 
  
  for(n in 1:N){
    for(k in 1:K){
      if(obs_mask[n, k]){
        real gamma = inv_omega[k]/(1 - inv_omega[k]) ; //the beta par of neg_bin
        real xi_i_k = gamma * exp(scaled_alpha[n] + scaled_beta[k])  ;
        y[n,k] ~ neg_binomial(xi_i_k, gamma); 
      }
    }
  }
}



// generated quantities {
//   array[N] int y_sum;
//   array[N, K] int y_sim;
//   array[N,K] real log_lik;
//   for (n in 1:N) {
//     for (k in 1:K) {
//       real gamma = inv_omega[k]/(1 - inv_omega[k]) ; //the beta par of neg_bin
//       real xi_i_k = gamma * exp(scaled_alpha[n] + scaled_beta[k])  ;
//       y_sim[n, k] = neg_binomial_rng(xi_i_k, gamma) ;
//       log_lik[n, k] = neg_binomial_lpmf(y[n, k] | xi_i_k, gamma);
//     }
//     y_sum[n] = sum( y_sim[n] );
//   }
// 
// }
