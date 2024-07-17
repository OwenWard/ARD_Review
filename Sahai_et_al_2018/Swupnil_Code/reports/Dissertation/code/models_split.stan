functions {
  // Flat luminosity integral. 
  real[] luminosity_flat(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    // theta[1] = Omega_m ; theta[2] = w
    real dydt[1];
    dydt[1] = (
                theta[1] * (1 + t)^3 +
                (1 - theta[1]) * (1 + t)^(3 * (1 + theta[2]))
              )^(-0.5);
    return dydt;
  }
  
  // Curved luminosity integral. 
  real[] luminosity_curved(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    // theta[1] = Omega_m ; theta[2] = Omega_L
    real dydt[1];
    dydt[1] = (
                theta[1] * (1 + t)^3 +
                theta[2] +
                (1 - theta[1] - theta[2]) * (1 + t)^2
              )^(-0.5);
    return dydt;
  }
  
  // Sinn function.
  real sinn(real Omega_k_, real x_) {
    real sqrt_Omega_k_;
    sqrt_Omega_k_ = sqrt(fabs(Omega_k_));
    if(Omega_k_ == 0) {
      return(x_);
    }
    else if(Omega_k_ > 0) {
      return(sinh(sqrt_Omega_k_ * x_)/sqrt_Omega_k_);
    }
    else {
      return(sin(sqrt_Omega_k_ * x_)/sqrt_Omega_k_);
    }
  }
  
  // Luminosity Distance.
  real lum_dist(real c_, real H_0_, real z_, real Omega_k_, real int_) {
    real dL_;
    dL_ = c_/H_0_ * (1 + z_) * sinn(Omega_k_, int_);
    return dL_;
  }
  
  // Distance modulus.
  real dist_mod(real c_, real H_0_, real z_, real Omega_k_, real int_) {
    real mu_;
    mu_ = 25 + 5 * log10(lum_dist(c_, H_0_, z_, Omega_k_, int_));
    return mu_;
  }
}

data {
  int<lower=1> N;                     # number of data points
  int<lower=1> N_sub;                 # number of data points for the third variable
  
  // Data
  vector[N] z_cmb;                    # observed CMB frame redshift
  vector[N] z_hel;                    # observed heliocntric redshift
  vector[N] mB_hat;                   # observed B-band peak magnitude
  vector[N] x1_hat;                   # observed stretch correction
  vector[N] cl_hat;                   # observed color correction
  matrix[3,3] C_hat[N];               # observed covariance for mag, stretch, color
  vector[N_sub] meth_hat;             # observed metallicity
  vector[N_sub] rate_hat;             # observed log star formation rate
  vector[N_sub] age_hat;              # observed log age
  
  // Constants
  real c;                             # speed of light in km/s
  real H_0;                           # dimensionless Hubble parameter today
  
  // Universe Assumption
  int<lower=0, upper=1> is_curved;    # whether or not the universe is curved
  int<lower=0, upper=4> third_var;    # 2 = metallicity; 3 = formation rate; 4 = age
  
  // Prior Centers and Scales
  vector[12] pri_mu;                  # means for the priors
  vector[12] pri_sd;                  # sds for the priors
}

transformed data {
  real x_r[0];                        # real inputs to the luminosity integral
  int x_i[0];                         # integer inputs to the luminosity integral
  real ts[N,2];                       # limits of the luminosity integral
  real t0;                            # initial limit value for the ode integration
  real y0[1];                         # value of integral at initial limit value t0
  for(n in 1:N) {
    ts[n,1] = 0;
    ts[n,2] = z_cmb[n];
  }
  t0 = -0.5;  
  y0[1] = 0;
}

parameters {
  // Cosmological Parameters
  real<lower=0,upper=2> Omega_m;           # total matter density (baryonic + dark)
  real<lower=0,upper=2> Omega_L;           # dark energy density
  real<lower=-2,upper=0> w;                # dark energy constant
  
  // Covariates
  real<lower=0> alpha;                     # coefficient of stretch covariate 
  real<lower=0> beta_0;                    # coefficient of color covariate  
  
  // Covariates for Subset
  real<lower=0> alpha_1;                   # coefficient of stretch covariate 
  real<lower=0> beta_1;                    # coefficient of color covariate
  real gamma;                              # coefficient of third covariate
  
  // Population-Level
  real M_0;                                # mean absolute magnitude (M_0^epsilon)
  real<lower=0> sigma_res;                 # sd of absolute magnitude
  real x1_0;                               # mean stretch correction (x1*)
  real log_R_x1;                           # log sd of stretch correction (R_x1)     
  real cl_0;                               # mean color correction (c*)
  real log_R_cl;                           # log sd of color correction (R_c)  
  
  // Local-Level
  vector[N] M;                             # true absolute magnitude
  vector[N] x1;                            # true stretch correction
  vector[N] cl;                            # true color correction
}

transformed parameters {
  real theta[2];                           # vectorized cosmological parameters
  real<lower=0> R_x1;                      # sd of stretch correction (R_x1)    
  real<lower=0> R_cl;                      # sd of color correction (R_c) 
  // Cosmological Parameters
  theta[1] = Omega_m;
  if(is_curved) {
    theta[2] = Omega_L;
  }
  else {
    theta[2] = w;
  }
  // Population Level Parameters
  R_x1 = exp(log_R_x1);
  R_cl = exp(log_R_cl);
}

model {
  real lum_limits[2,1];    # integrated luminosity limits
  real int_lum;            # integrated luminosity
  vector[N] mu;            # distance modulus
  vector[N] mB;            # true B-band peak magnitude
  vector[3] salt2_hat[N];  # observed salt2 (script D_hat)
  vector[3] salt2[N];      # true salt2 (script D)
  
  ## True B-band Peak Magnitude
  // Curved luminosity integration
  if(is_curved) {
    for(n in 1:N) {
      lum_limits = integrate_ode_rk45(luminosity_curved, y0, t0, ts[n,], theta, x_r, x_i);
      int_lum = lum_limits[2,1] - lum_limits[1,1];
      mu[n] = dist_mod(c, H_0, z_hel[n], 1 - Omega_m - Omega_L, int_lum);
    }
  }
  // Flat luminosity integration
  else {
    for(n in 1:N) {
      lum_limits = integrate_ode_rk45(luminosity_flat, y0, t0, ts[n,], theta, x_r, x_i);
      int_lum = lum_limits[2,1] - lum_limits[1,1];
      mu[n] = dist_mod(c, H_0, z_hel[n], 0, int_lum);
    }
  }
  
  // Baseline Model
  for(n in 1:N) {
    if (n > N_sub) {
      mB[n] = mu[n] - alpha * x1[n] + beta_0 * cl[n] + (M[n] - 19.3);
    }
    else {
      mB[n] = mu[n] - alpha_1 * x1[n] + beta_1 * cl[n] + (M[n] - 19.3);
      if (third_var == 2) {
        mB[n] = mB[n] + gamma * meth_hat[n];
      }
      else if (third_var == 3) {
        mB[n] = mB[n] + gamma * rate_hat[n];
      }
      else if (third_var == 4) {
        mB[n] = mB[n] + gamma * age_hat[n];
      }
    }
  }
  
  ## Priors
  // Cosmological Parameters
  Omega_m   ~ normal(pri_mu[1], pri_sd[1]);
  Omega_L   ~ normal(pri_mu[2], pri_sd[2]);
  w         ~ normal(pri_mu[3], pri_sd[3]);      
  // Covariates
  alpha     ~ normal(pri_mu[4], pri_sd[4]);
  beta_0    ~ normal(pri_mu[5], pri_sd[5]);
  alpha_1   ~ normal(pri_mu[4], pri_sd[4]);
  beta_1    ~ normal(pri_mu[5], pri_sd[5]);
  gamma     ~ normal(pri_mu[6], pri_sd[6]);
  // Population-Level
  M_0       ~ normal(pri_mu[7], pri_sd[7]);
  sigma_res ~ normal(pri_mu[8], pri_sd[8]);
  x1_0      ~ normal(pri_mu[9], pri_sd[9]);      
  log_R_x1  ~ normal(pri_mu[10], pri_sd[10]);
  cl_0      ~ normal(pri_mu[11], pri_sd[11]);
  log_R_cl  ~ normal(pri_mu[12], pri_sd[12]);
  
  // Local-Level
  M ~ normal(M_0, sigma_res);
  x1 ~ normal(x1_0, R_x1);
  cl ~ normal(cl_0, R_cl);
  
  ## Likelihood
  for(n in 1:N) {
    salt2_hat[n,1] = mB_hat[n];
    salt2_hat[n,2] = x1_hat[n];
    salt2_hat[n,3] = cl_hat[n];
    salt2[n,1] = mB[n];
    salt2[n,2] = x1[n];
    salt2[n,3] = cl[n];
    salt2_hat[n] ~ multi_normal(salt2[n], C_hat[n]);
  }
}
