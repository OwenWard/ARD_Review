# stan model code to calculate degree if responses are "how many X do you know?" (i.e. not binary)
degree_combined_mixing_code = '
  data {
    int<lower=1> E;     // number of ego groups
    int<lower=1> A;     // number of alter groups
    int<lower=1> N;     // number of respondents

    int<lower=1> J;     // number of questions for survery 1
    int<lower=2> K;     // number of questions for survery 2
    
    int y[N,J];         // responses for survery 1
    int z[N,K];         // responses for survery 2

    int ego[N];         // ego labels for respondents
    matrix[A,J] BetaY;  // alter proportions for survey 1
    matrix[A,K] BetaZ;  // alter proportions for survey 2
    
    vector[2] theta_d;  // prior mu and sigma for log_d
    vector[2] theta_o;  // prior alpha and beta for inv_omega
    vector[A] alpha;    // prior for mixing matrix rows
  }

  parameters {
    real<upper=8.5> log_d[N];               // log respondent degrees
    real<lower=0,upper=1> inv_omegaY[J];    // inverse over-dispersion for survey 1
    real<lower=0,upper=1> inv_omegaZ[K];    // inverse over-dispersion for survey 2
    simplex[A] M[E];                        // mixing matrix
  }

  transformed parameters {
    real<lower=0> d[N];
    real<lower=0> omegaY[J];
    real<lower=0> omegaZ[K];

    for(n in 1:N) {
      d[n] <- exp(log_d[n]);
    }
    for(j in 1:J) {
      omegaY[j] <- inv(inv(inv_omegaY[j]) - 1);
    }
    for(k in 1:K) {
      omegaZ[k] <- inv(inv(inv_omegaZ[k]) - 1);
    }
  }

  model {
    real mu_nj;
    real mu_nk;

    log_d ~ normal(theta_d[1], theta_d[2]);
    inv_omegaY ~ beta(theta_o[1], theta_o[2]);
    inv_omegaZ ~ beta(theta_o[1], theta_o[2]);
    for(e in 1:E){
      M[e] ~ dirichlet(alpha);
    }
    
    for(n in 1:N) {
      for(j in 1:J) {
        if(y[n,j]>=0) {
          mu_nj <- d[n] * dot_product(M[ego[n]], sub_col(BetaY, 1, j, A));
          y[n,j] ~ neg_binomial(mu_nj * omegaY[j], omegaY[j]);
        }
      }
      for(k in 1:K) {
        if(z[n,k]>=0) {
          mu_nk <- d[n] * dot_product(M[ego[n]], sub_col(BetaZ, 1, k, A));
          z[n,k] ~ neg_binomial(mu_nk * omegaZ[k], omegaZ[k]);
        }
      }
    }
  }
'