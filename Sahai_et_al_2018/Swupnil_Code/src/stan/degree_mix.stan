# stan model code to calculate degree if responses are "how many X do you know?" (i.e. not binary)
data {
    int<lower=1> E;   // number of ego groups
    int<lower=1> A;   // number of alter groups
    int<lower=1> K;   // number of names
    int<lower=1> N;   // number of respondents
    
    int y[N,K];
    int ego[N];
    matrix[A,K] Beta;
    
    vector[2] theta_d;  // prior mu and sigma for log_d
    vector[2] theta_o;  // prior alpha and beta for inv_omega
    vector<lower=0>[A] alpha;    // prior for mixing matrix rows
}

parameters {
    real<upper=8.5> log_d[N];               // log respondent degrees
    real<lower=0,upper=1> inv_omega[K];     // inverse over-dispersion
    simplex[A] M[E];                        // mixing matrix
}

transformed parameters {
    real<lower=0> d[N];
    real<lower=0> omega[K];

    for(n in 1:N){
      d[n] = exp(log_d[n]);
    }
    for(k in 1:K){
      omega[k] = inv(inv(inv_omega[k]) - 1);
    }
}

model {
    real mu_nk;

    log_d ~ normal(theta_d[1], theta_d[2]);
    inv_omega ~ beta(theta_o[1], theta_o[2]);
    for(e in 1:E){
      M[e] ~ dirichlet(alpha);
    }
    
    for(k in 1:K) {
      for(n in 1:N) {
        if(y[n,k]>=0) {
          mu_nk = d[n] * dot_product(M[ego[n]], sub_col(Beta, 1, k, A));
          y[n,k] ~ neg_binomial(mu_nk * omega[k], omega[k]);
        }
      }
    }
}