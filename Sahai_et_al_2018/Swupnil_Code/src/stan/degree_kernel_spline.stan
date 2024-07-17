// stan model code to calculate degree if responses are "how many X do you know?" (i.e. not binary)

functions {
  // Returns the integration of the product of two normals
  real norm_prod_int(int mu1, real mu2, real var1, real var2) {
    real K;
    K = 1 / sqrt(2 * 3.141593 * (var1 + var2))
        * exp(-(mu1 - mu2) ^ 2 / (2 * (var1 + var2)));
    return K;
  }
  
  // Builds the B spline with specified parameters
  
  vector build_b_spline(array[] real t, array[] real ext_knots, int ind,
                        int order) {
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
    if (order == 1) {
      for (i in 1 : size(t))  // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind + 1]);
    } else {
      if (ext_knots[ind] != ext_knots[ind + order - 1]) 
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t)))
             / (ext_knots[ind + order - 1] - ext_knots[ind]);
      if (ext_knots[ind + 1] != ext_knots[ind + order]) 
        w2 = 1
             - (to_vector(t) - rep_vector(ext_knots[ind + 1], size(t)))
               / (ext_knots[ind + order] - ext_knots[ind + 1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order - 1)
                 + w2 .* build_b_spline(t, ext_knots, ind + 1, order - 1);
    }
    return b_spline;
  }
}
data {
  int<lower=0> K_1; // number of names
  int<lower=0> K_2; // number of occupations
  int<lower=1> N; // number of respondents
  int<lower=1> N_K; // number of knots for the bandwidth spline
  int<lower=1> N_S; // number of grid points to evaluate spline on
  int<lower=1> D; // degree of bandwidth spline (order - 1)
  real age_mean; // average age
  
  array[N, K_1 + K_2] int Y; // responses to how many X do you know
  array[N] real w; // response survey weights
  array[N] int age; // age of respondents
  array[N] int<lower=1> g_n; // gender of respondents: 1=female 2=male
  array[K_1 + K_2 * 2] int<lower=1> g_k; // genders of names: 1=female 2=male (names have one each, but occs have two genders each)
  
  array[K_1 + K_2 * 2] real<lower=0> mu_k; // age mean for each gender-group k
  array[K_1 + K_2 * 2] real<lower=0> sigma_k; // age standard deviation for each gender-group k
  array[K_1 + K_2 * 2] real<lower=0> sum_prob_k; // sum population proportion of each gender-group k
  
  vector[N_K] knots; // sequence of knots for the bandwidth spline
  array[N_S] real X; // age grid over which to evaluate spline
  
  real mu_d; // prior mean for the log_degrees
  real<lower=0> sigma_d; // prior sd for the log_degrees
  real<lower=0> alpha_omega; // prior alpha for inv_omega
  real<lower=0> beta_omega; // prior beta for inv_omega
  vector<lower=0>[2] alpha_rho; // prior alpha for gender mixing rows
  vector[3] mu_beta; // prior mean for beta
  vector<lower=0>[3] sigma_beta; // prior sd for beta
  
  real recall_power; // 0 = no recall adjustment, 0.5 = recall adjustment
  real degree_regression; // 0 = simple prior, 1 = agesex regression prior
}
transformed data {
  int K; // total number of groups
  int N_B; // total number of B-splines              
  matrix[N_K + D - 1, N_S] B; // matrix of B-splines
  vector[D + N_K] ext_knots_temp;
  vector[2 * D + N_K] ext_knots; // set of extended knots
  K = K_1 + K_2;
  N_B = N_K + D - 1;
  ext_knots_temp = append_row(rep_vector(knots[1], D), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[N_K], D));
  for (ind in 1 : N_B) 
    B[ind,  : ] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots),
                                               ind, D + 1));
  B[N_B, N_S] = 1;
}
parameters {
  vector[N] log_d; // log respondent degrees
  array[K] real<lower=0, upper=1> inv_omega; // inverse group overdispersions
  array[2] simplex[2] rho; // gender mixing rates
  vector[3] beta; // degree regression coefs (third is logged)
  real log_eta; // log sd of degree
  
  array[2, 2] row_vector[N_B] a_raw; // raw spline coefficients
  array[2, 2] real a0; // spline intercept
  array[2, 2] real<lower=0> tau; // spline coefficient prior sd
}
transformed parameters {
  array[N] real<lower=0> d; // respondent degrees
  array[K] real<lower=0> omega; // group overdispersions
  real<lower=0> eta; // sd of degree
  array[2, 2] row_vector[N_B] a; // spline coefficients
  array[2, 2] vector[N_S] lambda; // spline evaluated at the age grid
  for (n in 1 : N) {
    d[n] = exp(log_d[n]);
  }
  for (k in 1 : K) {
    omega[k] = inv(inv(inv_omega[k]) - 1);
  }
  eta = exp(log_eta);
  for (i in 1 : 2) {
    for (j in 1 : 2) {
      a[i, j] = a_raw[i, j] * tau[i, j];
      lambda[i, j] = exp(a0[i, j] + to_vector(a[i, j] * B));
    }
  }
}
model {
  vector[N] ex_log_d;
  real mu_nk;
  int age_;
  
  // Parameter Priors
  for (p in 1 : 3) {
    beta[p] ~ normal(mu_beta[p], sigma_beta[p]);
  }
  log_eta ~ normal(-0.7, 0.1);
  inv_omega ~ beta(alpha_omega, beta_omega);
  for (i in 1 : 2) {
    rho[i] ~ dirichlet(alpha_rho);
    for (j in 1 : 2) {
      a_raw[i, j] ~ normal(0, 2);
      a0[i, j] ~ normal(0, 2);
      tau[i, j] ~ normal(0, 2);
    }
  }
  
  // Degree Priors
  if (degree_regression == 1) {
    for (n in 1 : N) {
      ex_log_d[n] = beta[1] + beta[2] * (g_n[n] - 1)
                    - exp(beta[3]) * ((age[n] - age_mean) / age_mean) ^ 2;
    }
    log_d ~ normal(ex_log_d, eta);
  } else {
    log_d ~ normal(mu_d, sigma_d);
  }
  
  // Likelihood
  for (k in 1 : K) {
    for (n in 1 : N) {
      if (Y[n, k] >= 0) {
        mu_nk = rho[g_n[n], g_k[k]] * sum_prob_k[k]
                * norm_prod_int(age[n], mu_k[k],
                                lambda[g_n[n], g_k[k]][age[n] - 17],
                                sigma_k[k] ^ 2);
        // For occupations, we add a term for other gender
        if (k > K_1) {
          mu_nk = mu_nk
                  + rho[g_n[n], g_k[k + K_2]] * sum_prob_k[k + K_2]
                    * norm_prod_int(age[n], mu_k[k + K_2],
                                    lambda[g_n[n], g_k[k + K_2]][age[n] - 17],
                                    sigma_k[k + K_2] ^ 2);
        }
        if (recall_power > 0) {
          mu_nk = pow(mu_nk, 1 - recall_power)
                  * exp((-6 - exp(-7) / mu_nk) * recall_power);
        }
        mu_nk = mu_nk * d[n];
        target += w[n]
                  * neg_binomial_lpmf(Y[n, k] | omega[k] * mu_nk, omega[k]);
      }
    }
  }
}


