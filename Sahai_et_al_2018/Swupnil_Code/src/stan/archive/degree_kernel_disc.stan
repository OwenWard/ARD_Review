# stan model code to calculate degree if responses are "how many X do you know?" (i.e. not binary)
degree_kernel_discrete_code = '
  functions {
    real kernel(real ai, real aj, real lambda) {
      real y;
      y <- exp(-0.5 * square(ai-aj)/lambda);
      return y;
    }
    vector kernel_vec(real ai, vector s, real lambda, int J) {
      vector[J] y;
      for(j in 1:J) {
        y[j] <- kernel(ai, s[j], lambda);
      }
      return y;
    }
  }

  data {
    int<lower=1> K;   // number of names
    int<lower=1> N;   // number of respondents
    int<lower=1> J;   // number of alter ages
    
    int y[N,K];       // responses to how many X do you know
    real age[N];      // age of respondents
    int sex[N];       // gender of respondents: 2=male, 1=female
    matrix[K,J] Beta; // name age proportions
    vector[J] s;        // alter ages
    
    vector[2] theta_d;  // prior mu and sigma for log_d
    vector[2] theta_o;  // prior alpha and beta for inv_omega
    vector[2] theta_l;  // prior mu and sigma for log_lambda
    vector[2] theta_m;  // prior alpha for gender mixing rows
  }

  parameters {
    real<upper=8.5> log_d[N];               // log respondent degrees
    real<lower=0,upper=1> inv_omega[K];     // name inverse over-dispersions
    matrix[2,2] log_lambda;                 // log kernel scales
    simplex[2] M[2];                        // gender mixing matrix
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
    real mu_ik;
    int g_k;
    int g_i;

    log_d ~ normal(theta_d[1], theta_d[2]);
    inv_omega ~ beta(theta_o[1], theta_o[2]);
    for(i in 1:2) {
      M[i] ~ dirichlet(theta_m);
      for(j in 1:2) {
        log_lambda[i,j] ~ normal(theta_l[1], theta_l[2]);
      }
    }
    
    for(k in 1:K) {
      if(k<=6) {
        g_k <- 1;
      } else {
        g_k <- 2;
      }

      for(i in 1:N) {
        if(y[i,k]>=0) {
          g_i <- sex[i];
          mu_ik <- dot_product(sub_row(Beta,k,1,J), kernel_vec(age[i], s, lambda[g_i,g_k], J));
          mu_ik <- mu_ik * M[g_i][g_k] * d[i] / sum(kernel_vec(age[i], s, lambda[g_i,g_k], J));
          y[i,k] ~ neg_binomial(omega[k] * mu_ik, omega[k]);
        }
      }
    }
  }
'