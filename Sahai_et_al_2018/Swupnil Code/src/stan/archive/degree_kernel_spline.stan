# stan model code to calculate degree if responses are "how many X do you know?" (i.e. not binary)

functions {
  // Returns the integration of the product of two normals
  real norm_prod_int(int mu1, real mu2, real var1, real var2) {
    real K;
    K = 1/sqrt(2 * 3.141593 * (var1 + var2)) * exp(-(mu1 - mu2)^2/(2 * (var1 + var2)));
    return K;
  }
  
  // Builds the B spline with specified parameters
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1;
    vector[size(t)] w2;
    w1 = rep_vector(0, size(t));
    w2 = rep_vector(0, size(t));
    if (order==1) {
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    }
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int<lower=1> K;                         // number of groups
  int<lower=1> N;                         // number of respondents
  int<lower=1> KG;                        // number of groups * genders
  int<lower=1> N_K;                       // number of knots for the bandwidth spline
  int<lower=1> N_S;                       // number of grid points to evaluate spline on
  int<lower=1> D;                         // degree of bandwidth spline (order - 1)
  real age_mean;                          // average age

  int Y[N,K];                             // responses to how many X do you know
  real w[N];                              // response survey weights
  int age[N];                             // age of respondents
  int<lower=1> g_n[N];                    // gender of respondents: 1=female 2=male
  int<lower=1> g_k[KG];                   // genders of names: 1=female 2=male

  real<lower=0> mu_k[KG];                 // age mean for each group k
  real<lower=0> sigma_k[KG];              // age standard deviation for each group k
  real<lower=0> sum_prob_k[KG];           // sum population proportion of each group k
  
  vector[N_K] knots;                      // sequence of knots for the bandwidth spline
  real X[N_S];                            // age grid over which to evaluate spline

  real mu_d;                              // prior mean for the log_degrees
  real<lower=0> sigma_d;                  // prior sd for the log_degrees
  real<lower=0> alpha_omega;              // prior alpha for inv_omega
  real<lower=0> beta_omega;               // prior beta for inv_omega
  vector<lower=0>[2] alpha_rho;           // prior alpha for gender mixing rows
  vector[3] mu_beta;                      // prior mean for beta
  vector<lower=0>[3] sigma_beta;          // prior sd for beta

  real recall_power;                      // 0 = no recall adjustment, 0.5 = recall adjustment
  real degree_regression;                 // 0 = simple prior, 1 = agesex regression prior
}

transformed data {
  int N_B;                                // total number of B-splines              
  matrix[N_K + D - 1, N_S] B;             // matrix of B-splines
  vector[D + N_K] ext_knots_temp;
  vector[2*D + N_K] ext_knots;            // set of extended knots
  N_B = N_K + D - 1;    
  ext_knots_temp = append_row(rep_vector(knots[1], D), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[N_K], D));
  for (ind in 1:N_B)
    B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, D + 1));
  B[N_B, N_S] = 1;
}

parameters {
  vector[N] log_d;                        // log respondent degrees
  real<lower=0,upper=1> inv_omega[K];     // inverse group overdispersions
  simplex[2] rho[2];                      // gender mixing rates
  vector[3] beta;                         // degree regression coefs (third is logged)
  real log_eta;                           // log sd of degree
  
  row_vector[N_B] a_raw[2,2];             // raw spline coefficients
  real a0[2,2];                           // spline intercept
  real<lower=0> tau[2,2];                 // spline coefficient prior sd
}

transformed parameters {
  real<lower=0> d[N];                     // respondent degrees
  real<lower=0> omega[K];                 // group overdispersions
  real<lower=0> eta;                      // sd of degree
  row_vector[N_B] a[2,2];                 // spline coefficients
  vector[N_S] lambda[2,2];                // spline evaluated at the age grid
  for(n in 1:N) {
    d[n] = exp(log_d[n]);
  }
  for(k in 1:K) {
    omega[k] = inv(inv(inv_omega[k]) - 1);
  }
  eta = exp(log_eta);
  for(i in 1:2) {
    for(j in 1:2) {
      //a[i,j][1] = a_raw[i,j][1];
      //for(b in 2:N_B)
      //  a[i,j][b] = a[i,j][b-1] + a_raw[i,j][b] * tau[i,j];
      a[i,j] = a_raw[i,j] * tau[i,j];
      lambda[i,j] = exp(a0[i,j] + to_vector(a[i,j] * B));
    }
  }
}

model {
  vector[N] ex_log_d;
  real mu_nk;
  int age_;

  // Parameter Priors
  for(p in 1:3) {
    beta[p] ~ normal(mu_beta[p], sigma_beta[p]);
  }
  log_eta ~ normal(-0.7, 0.1);
  inv_omega ~ beta(alpha_omega, beta_omega);
  for(i in 1:2) {
    rho[i] ~ dirichlet(alpha_rho);
    for(j in 1:2) {
      a_raw[i,j] ~ normal(0, 2);
      a0[i,j] ~ normal(0, 2);
      tau[i,j] ~ normal(0, 2);
    }
  }

  // Degree Priors
  if(degree_regression == 1) {
    for(n in 1:N) {
      ex_log_d[n] = beta[1] + 
                    beta[2] * (g_n[n]-1) - 
                    exp(beta[3]) * ((age[n] - age_mean)/age_mean)^2;
    }
    log_d ~ normal(ex_log_d, eta);
  }
  else {
    log_d ~ normal(mu_d, eta);
  }

  // Likelihood
  for(k in 1:K) {
    for(n in 1:N) {
      if(Y[n,k] >= 0) {
        mu_nk = rho[g_n[n], g_k[k]] *  
                sum_prob_k[k] * 
                norm_prod_int(age[n], 
                              mu_k[k], 
                              lambda[g_n[n], g_k[k]][age[n]-17], 
                              sigma_k[k]^2);
        if(KG == K*2) {
          mu_nk = mu_nk + rho[g_n[n], g_k[k+K]] *  
                          sum_prob_k[k+K] *
                          norm_prod_int(age[n], 
                                        mu_k[k+K], 
                                        lambda[g_n[n], g_k[k+K]][age[n]-17], 
                                        sigma_k[k+K]^2);
        }
        if(recall_power > 0) {
          mu_nk = pow(mu_nk, 1 - recall_power) * exp((-6 - exp(-7)/mu_nk) * recall_power);
        }
        mu_nk = mu_nk * d[n];
        target += w[n] * neg_binomial_lpmf(Y[n,k] | omega[k] * mu_nk, omega[k]);
      }
    }
  }
}
