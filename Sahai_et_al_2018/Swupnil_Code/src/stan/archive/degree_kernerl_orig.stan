# stan model code to calculate degree if responses are "how many X do you know?" (i.e. not binary)

degree_kernel_continuous_code = '
functions {
  real norm_prod_int(real a, real mu, real lambda, real sigma) {
    real y;
    y <- 1/sqrt(2 * 3.141593 * (lambda + square(sigma))) * exp(-square(a - mu)/(2 * (lambda + square(sigma))));
    return y;
  }
}

data {
  int<lower=1> K;                         // number of names
  int<lower=1> N;                         // number of respondents
  
  int y[N,K];                             // responses to how many X do you know
  real age[N];                            // age of respondents
  int<lower=1> g_n[N];                    // gender of respondents: 1=female 2=male
  int<lower=1> g_k[K];                    // genders of names: 1=female 2=male

  real<lower=0> mu_k[K];                  // age mean for each name k
  real<lower=0> sigma_k[K];               // age standard deviation for each name k
  real<lower=0> sum_prop_k[K];            // summation of age-wise population of each name k divided by age-wise population
  
  real mu_d;                              // prior mu for the degrees
  real<lower=0> sigma_d;                  // prior sigma for the degrees
  real mu_lambda;                         // prior mu for log_lambda
  real<lower=0> sigma_lambda;             // prior sigma for log_lambda
  real<lower=0> alpha_omega;              // prior alpha for inv_omega
  real<lower=0> beta_omega;               // prior beta for inv_omega
  vector<lower=0>[2] alpha_rho;           // prior alpha for gender mixing rows
  
  real recallPower;                       // 0 = no recall adjustment, 0.5 = recall adjustment
}

parameters {
  real<upper=8.5> log_d[N];               // log respondent degrees
  real<lower=0,upper=1> inv_omega[K];     // inverse over-dispersion
  matrix[2,2] log_lambda;                 // log kernel scales
  simplex[2] rho[2];                      // gender mixing rates
}

transformed parameters {
  real<lower=0> d[N];
  real<lower=0> omega[K];
  real<lower=0> lambda[2,2];
  
  for(n in 1:N){
    d[n] <- exp(log_d[n]);
  }
  for(k in 1:K){
    omega[k] <- inv(inv(inv_omega[k]) - 1);
  }
  for(i in 1:2){
    for(j in 1:2){
      lambda[i,j] <- exp(log_lambda[i,j]);
    }
  }
}

model {
  real mu_nk;
  
  // Priors
  log_d ~ normal(mu_d, sigma_d);
  inv_omega ~ beta(alpha_omega, beta_omega);
  for(i in 1:2) {
    rho[i] ~ dirichlet(alpha_rho);
    for(j in 1:2) {
      log_lambda[i,j] ~ normal(mu_lambda, sigma_lambda);
    }
  }
  
  // Likelihood
  for(k in 1:K) {
    for(n in 1:N) {
      if(y[n,k] >= 0) {

        mu_nk <- sum_prop_k[k];
        mu_nk <- mu_nk * rho[g_n[n], g_k[k]];
        mu_nk <- mu_nk * norm_prod_int(age[n], mu_k[k], lambda[g_n[n], g_k[k]], sigma_k[k]);
        
        if(recallPower > 0) {
          mu_nk <- pow(mu_nk, 1 - recallPower) * exp((-6 - exp(-7)/mu_nk) * recallPower);
        }
        
        mu_nk <- mu_nk * d[n];
        y[n,k] ~ neg_binomial(omega[k] * mu_nk, omega[k]);
      }
    }
  }
}
'