# stan model code to calculate degree if responses are "how many X do you know?" (i.e. not binary)

functions {
  // Returns the integration of the product of two normals
  real norm_prod_int(real mu1, real mu2, real var1, real var2) {
    real K;
    K = 1/sqrt(2 * 3.141593 * (var1 + var2)) * exp(-(mu1 - mu2)^2/(2 * (var1 + var2)));
    return K;
  }

  // Returns the expectation of the product of two normals
  real norm_prod_ex(real mu1, real mu2, real var1, real var2) {
    real ex_;
    ex_ = ((mu1 * var2) + (mu2 * var1))/(var1 + var2);
    return(ex_);
  }

  // Returns the variance of the product of two normals
  real norm_prod_var(real mu1, real mu2, real var1, real var2) {
    real var_;
    var_ = (var1 * var2)/(var1 + var2);
    return(var_);
  }

}

data {
  int<lower=1> K;                         // number of names
  int<lower=1> N;                         // number of respondents
  real age_mean;                          // average age
  real degree_mean;                       // average degree

  int y[N,K];                             // responses to how many X do you know
  real age[N];                            // age of respondents
  int<lower=1> g_n[N];                    // gender of respondents: 1=female 2=male
  int<lower=1> g_k[K];                    // genders of names: 1=female 2=male

  real<lower=0> mu_k[K];                  // age mean for each group k
  real<lower=0> sigma_k[K];               // age standard deviation for each group k
  real<lower=0> prop_k[K];                // population of each group k divided by total population

  real mu_d;                              // prior mu for the log_degrees
  real<lower=0> sigma_d;                  // prior sigma for the log_degrees
  real<lower=0> alpha_omega;              // prior alpha for inv_omega
  real<lower=0> beta_omega;               // prior beta for inv_omega
  real mu_lambda;                         // prior mu for log_lambda
  real<lower=0> sigma_lambda;             // prior sigma for log_lambda
  vector<lower=0>[2] alpha_rho;           // prior alpha for gender mixing rows

  real recallPower;                       // 0 = no recall adjustment, 0.5 = recall adjustment
  real isDegreeRegression;                // 0 = constant degree assumption, 1 = non constant degree
}

parameters {
  vector[N] log_d;                        // log respondent degrees
  real<lower=0,upper=1> inv_omega[K];     // inverse over-dispersion
  matrix[2,2] log_lambda;                 // log kernel scales
  simplex[2] rho[2];                      // gender mixing rates
  vector[3] beta;                         // degree regression coefficients
  real log_eta_d;                         // log degree standard deviation
}

transformed parameters {
  real<lower=0> d[N];
  real<lower=0> omega[K];
  real<lower=0> lambda[2,2];
  real<lower=0> beta_3;  
  real<lower=0> eta_d;

  for(n in 1:N) {
    d[n] = exp(log_d[n]);
  }
  for(k in 1:K) {
    omega[k] = inv(inv(inv_omega[k]) - 1);
  }
  for(i in 1:2) {
    for(j in 1:2){
      lambda[i,j] = exp(log_lambda[i,j]);
    }
  }
  beta_3 = exp(beta[3]);
  eta_d = exp(log_eta_d);
}

model {
  real mu_nk;
  real ex_alter_degree;

  // Regression Coefficient Priors
  beta ~ normal(0,2);
  log_eta_d ~ normal(-0.7,0.1);

  // Priors
  inv_omega ~ beta(alpha_omega, beta_omega);
  for(i in 1:2) {
    rho[i] ~ dirichlet(alpha_rho);
    for(j in 1:2) {
      log_lambda[i,j] ~ normal(mu_lambda, sigma_lambda);
    }
  }

  // Degree Prior
  if(isDegreeRegression > 0) {
    for(n in 1:N) {
      log_d[n] ~ normal(beta[1] + beta[2] * (g_n[n]-1) - beta_3 * ((age[n] - age_mean)/age_mean)^2, eta_d);
    }
  }
  else {
    log_d ~ normal(mu_d, sigma_d);
  }

  // Likelihood
  for(k in 1:K) {
    for(n in 1:N) {
      if(y[n,k] >= 0) {

        mu_nk = prop_k[k];
        mu_nk = mu_nk * rho[g_n[n], g_k[k]];
        mu_nk = mu_nk * norm_prod_int(age[n], mu_k[k], lambda[g_n[n],g_k[k]], sigma_k[k]^2);

        //if(isDegreeRegression > 0) { //TODO: normalize by age_mean^2
        //  ex_alter_degree <- exp(beta[1] + beta[2] * (g_k[k] - 1) + (eta_d)^2/2);
        //  ex_alter_degree <- ex_alter_degree / sqrt(beta_3 / 3.141593);
        //  ex_alter_degree <- ex_alter_degree * norm_prod_int(age_mean, 
        //                                                     norm_prod_ex(age[n], mu_k[k], lambda[g_n[n],g_k[k]], sigma_k[k]^2), 
        //                                                     pow(2 * beta_3, -1), 
        //                                                     norm_prod_var(age[n], mu_k[k], lambda[g_n[n],g_k[k]], sigma_k[k]^2));
        //  mu_nk <- mu_nk * ex_alter_degree / degree_mean;
        //} 
        if(recallPower > 0) {
          mu_nk = pow(mu_nk, 1 - recallPower) * exp((-6 - exp(-7)/mu_nk) * recallPower);
        }

        mu_nk = mu_nk * d[n];
        y[n,k] ~ neg_binomial(omega[k] * mu_nk, omega[k]);
      }
    }
  }
}