// May 6th 2025
// Stan Model for the Overdispersed Binomial
// Model of McCormick and Zheng et al 2015
// modified to take varying priors for each comp of beta

functions {
  /**
   *  Log normalizing constant of the von Mises–Fisher distribution on S^p.
   *  The sphere lives in R^{p+1}, so (data) vectors have length p+1.
   *
   *  @param p       integer dimension of the sphere  (p ≥ 2)
   *  @param kappa   concentration parameter κ ≥ 0
   *  @return        log C_p(κ)
   *
   *  C_p(κ) = κ^{p/2 - 1} / [ (2π)^{p/2}  I_{p/2 - 1}(κ) ]
   *
   *  The implementation uses log_modified_bessel_first_kind()
   *  for numerical stability when κ is large.
   */
  real log_vmf_norm(int p, real kappa) {
    if (kappa == 0) {
    return   lgamma(0.5 * (p + 1))
           - log(2)
           - 0.5 * (p + 1) * log(pi());
    }
    real v          = 0.5 * p - 1.0;              // order of the Bessel fn
    real log_bessel = log_modified_bessel_first_kind(v, kappa);
    return v * log(kappa)
           - 0.5 * p * log(2 * pi())
           - log_bessel;
  }
  // real vmf_norm(int p, real kappa) {
  //   return exp(log_vmf_norm(p, kappa));
  // }
}
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
  vector<lower=0>[K] eta;
  real<lower=0> xi;
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
  vector[K] num_const_eta;
  vector[K] den_const_eta;
  for(k in 1:K){
    num_const_eta[k] = log_vmf_norm(p, eta[k]);
    // den_const_eta[k] = 
  }
  matrix[N, p] Z;
  for (n in 1:N)
    Z[n] = z[n]';            // transpose: vector  → row_vector → matrix row

  matrix[K, p] NU;
  for (k in 1:K)
    NU[k] = nu[k]'; 
  matrix[N,K] dot   = Z * NU';                 // BLAS call
  // keep numerical safety near ±1
  matrix[N,K] dclip = fmin(fmax(dot, -1 + 1e-12), 1 - 1e-12);
  matrix[N,K] theta = acos(dclip);  
  real num_part1 = log_vmf_norm(p, xi);
  real den_part1 = log_vmf_norm(p, 0);
}

// The model to be estimated..
model {
profile("priors") {
  xi ~ gamma(2, 1); // setting this for now
  eta ~ gamma(2, 1);
  sigma_beta ~ normal(0, 5);
  sigma_alpha ~ normal(0, 5);
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta); 
}
profile("likelihood"){

  for(n in 1:N){
    for(k in 1:K){
      real den_term = sqrt(xi ^2 +  eta[k]^2 + 2 * xi * eta[k] *
                            theta[n, k]);
      // real den_term = 1;
      real log_num = num_part1 + num_const_eta[k];//log_vmf_norm(p, eta[k]);
      real log_den = den_part1 + log_vmf_norm(p, den_term);
      // print("myvar = ", log_den);
      real gamma_ik = exp(log_num - log_den); //num/den ; //the hypersphere part
      real scaled_deg = exp(scaled_alpha[n]);
      y[n,k] ~ poisson(scaled_deg * exp(scaled_beta[k]) * gamma_ik);
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
