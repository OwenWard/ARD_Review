data {
  int<lower=1> E;   // number of ego groups
  int<lower=1> A;   // number of alter groups
  int<lower=1> K;   // number of names
  int<lower=1> N;   // number of respondents
  
  array[N, K] int y; //int y[N,K];
  array[N] int ego; //int ego[N];
  matrix[A,K] Beta;
  
  vector[2] theta_d;  // prior mu and sigma for log_d
  vector[2] theta_o;  // prior alpha and beta for inv_omega
  vector[A] alpha;    // prior for mixing matrix rows
  real<lower=0> p;    // power transformation: 0.5 or 1
}

parameters {
  array[N] real<upper=8.5> log_d; 
  //real<upper=8.5> log_d[N]; // respondent degrees
  array[K] real<lower=0, upper=1> inv_omega;
  //real<lower=0,upper=1> inv_omega[K];     // inverse over-dispersion
  array[E] simplex[A] M;
  //simplex[A] M[E];                        // mixing matrix
}

transformed parameters {
  array[N] real<lower=0> d;
  //real<lower=0> d[N];
  array[K] real<lower=0> omega;
  //real<lower=0> omega[K];
  
  for(n in 1:N){
    d[n] = exp(log_d[n]);
  }
  for(k in 1:K){
    omega[k] = inv(inv(inv_omega[k]) - 1);
  }
}

model {
  real mu_ik;

  log_d ~ normal(theta_d[1], theta_d[2]);
  inv_omega ~ beta(theta_o[1], theta_o[2]);
  for(e in 1:E){
    M[e] ~ dirichlet(alpha);
  }
  
  for(k in 1:K) {
    for(n in 1:N) {
      if(y[n,k]>=0) {
        mu_ik = omega[k] * d[n] * dot_product(M[ego[n]], sub_col(Beta, 1, k, A));
        mu_ik = pow(mu_ik, p);
        y[n,k] ~ neg_binomial(mu_ik, omega[k]);
      }
    }
  }
}


generated quantities {
  array[N] int out_degree;
  array[N, K] int y_sim;
  for (n in 1:N) {
    for (k in 1:K) {
      if(y[n,k]>=0) {
        real mu_ik = omega[k] * d[n] * dot_product(M[ego[n]], sub_col(Beta, 1, k, A));
        mu_ik = pow(mu_ik, p);
        y_sim[n,k] = neg_binomial_rng(mu_ik, omega[k]);
      }
    }
    out_degree[n] = sum( y_sim[n] );
  }
  
}
