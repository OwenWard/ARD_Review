# stan model code to calculate degree if responses are "how many X do you know?" (i.e. not binary)
degree_code = '
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
    vector[A] alpha;    // fixed random mixing rows
    real<lower=0> p;    // power transformation: 0.5 or 1
  }

  parameters {
    real<upper=8.5> log_d[N];               // respondent degrees
    vector<lower=0,upper=1>[K] inv_omega;   // inverse over-dispersion
  }

  transformed parameters {
    real<lower=0> d[N];
    real<lower=0> omega[K];

    for(n in 1:N){
      d[n] <- exp(log_d[n]);
    }
    for(k in 1:K){
      omega[k] <- inv(inv(inv_omega[k]) - 1);
    }
  }

  model {
    real mu_ik;

    log_d ~ normal(theta_d[1], theta_d[2]);
    inv_omega ~ beta(theta_o[1], theta_o[2]);

    for(k in 1:K) {
      for(n in 1:N) {
        mu_ik <- omega[k] * d[n] * dot_product(alpha, sub_col(Beta,1,k,A));
        mu_ik <- pow(mu_ik, p);
        y[n,k] ~ neg_binomial(mu_ik, omega[k]);
      }
    }
  }
'