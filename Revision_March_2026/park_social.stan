functions {

    real partial_sum(
        array[] int y_slice,
        int start,
        int end,
        array[] int ii,
        array[] int jj,
        real alpha,
        vector gamma,
        vector d,
        array[] real p,
        real phi
    ) {
  
        // index
        int len = end - start + 1;
        array[len] int ii_slice = ii[start:end];
        array[len] int jj_slice = jj[start:end];
        vector[len] d_slice = d[start:end];
        
        // accumulator
        real psum = 0.0;
        
        for (n in 1:len) {
        
            real mu = exp(alpha + gamma[ii_slice[n]] + log(p[jj_slice[n]]) - d_slice[n]);

            if (y_slice[n] == 0)
                psum += neg_binomial_2_lpmf(0 | mu, phi);
            else if (y_slice[n] == 1)
                psum += neg_binomial_2_lpmf(1 | mu, phi);
            else if (y_slice[n] == 2)
                psum += log_diff_exp(neg_binomial_2_lcdf(5 | mu, phi), 
                                     neg_binomial_2_lcdf(1 | mu, phi));
            else if (y_slice[n] == 3)
                psum += log_diff_exp(neg_binomial_2_lcdf(10 | mu, phi), 
                                     neg_binomial_2_lcdf(5 | mu, phi));
            else
                psum += neg_binomial_2_lccdf(10 | mu, phi);
          
        }
  
        return psum;
        
    }

}

data {
  
    int<lower=0> N;                 // No. of item-responses
    int<lower=0> I;                 // No. of Individuals.
    int<lower=0> J;                 // No. Items.
    int<lower=1> D;
    array[N] int<lower=0> ii;             // individual identifiers
    array[N] int<lower=0> jj;             // item identifiers
    array[N] int<lower=0> y;              // responses
    array[J] real<lower=0, upper=1> pops; // group offset
    array[J] real<lower=0> sds;           // std. dev. of priors
  
}

transformed data {
    
    // let stan automatically choose grainsize
    int grainsize = 1;

}

parameters {
  
    // Population parameters
    array[J] real<lower=0, upper=1> p; 
    
    // Gregariousness parameters
    real alpha;                
    vector[I] gamma_star;      
    real<lower=0> sigma_gamma; 
    
    // Group positions
    array[J - D - 1] vector[D] pos_j_star;
    vector<lower=0>[D] const_j;
    
    // Individual positions
    cholesky_factor_corr[D] L_Rho_i; 
    array[I] vector[D] pos_i_star;
    
    // overdispersion 
    real<lower=0> phi;
  
}

transformed parameters {
  
    vector[I] gamma;
    array[J] vector[D] pos_j;
    array[I] vector[D] pos_i;
    vector[N] d;
    
    gamma = gamma_star * sigma_gamma;
    
    // fixed stimuli
    for (d1 in 1:D) {
        for (d2 in 1:D) {
            if (d1==d2) {
                pos_j[d1, d2] = const_j[d1];
            } else {
                pos_j[d1, d2] = 0;
            }
        }
    }
    
    // position for rest (mean fixed at 0)
    {
        vector[D] cum_sum;
      
        cum_sum = rep_vector(0,D);
      
        for (cc in 1:D)
            cum_sum = cum_sum + pos_j[cc];
      
        for (cc in (D+1):(J-1)) {
            pos_j[cc] = 5 * pos_j_star[cc-D];
            cum_sum = cum_sum + pos_j[cc];
        }
      
        pos_j[J] = -cum_sum;
      
    }
  
    // individual positions
    for (i in 1:I)
        pos_i[i] = L_Rho_i * pos_i_star[i];
  
    // distances
    {
        
        vector[N] d_star;
        for (nn in 1:N) 
            d_star[nn] = log(squared_distance(pos_i[ii[nn]], pos_j[jj[nn]]));
        d = d_star - mean(d_star);
        
    }
  
}

model {
  
    // population estimates priors
    for (j in 1:J) 
        p[j] ~ normal(pops[j],sds[j]);
  
    // gregariousness priors
    alpha ~ normal(0, 5);
    gamma_star ~ normal(0, 1);
    sigma_gamma ~ normal(0, 3);
  
    // group position priors
    for (j in 1:(J-D-1)) 
        pos_j_star[j] ~ normal(0, 1);
    const_j ~ gamma(2, .1);
  
    // individual position prior
    for (i in 1:I)
       pos_i_star[i] ~ normal(0, 1);
    L_Rho_i ~ lkj_corr_cholesky(2);
  
    // overdispersion
    phi ~ normal(0, 3);
  
    // log-likelihood
    target += reduce_sum(partial_sum, y, grainsize, ii, jj, alpha, gamma, d, p, phi);
    
}

generated quantities {
  
  array[N] real log_lik;
  corr_matrix[D] Rho_i;
  
  // log-likelihood
  for (i in 1:N) {

    real mu;
    mu = exp(alpha + gamma[ii[i]] + log(p[jj[i]]) - d[i]);

    if (y[i] == 0)
        log_lik[i] = neg_binomial_2_lpmf(0 | mu, phi);
    else if (y[i] == 1)
       log_lik[i] = neg_binomial_2_lpmf(1 | mu, phi);
    else if (y[i] == 2)
        log_lik[i] = log_diff_exp(neg_binomial_2_lcdf(5 | mu, phi), 
                                  neg_binomial_2_lcdf(1|mu, phi));
    else if (y[i] == 3)
       log_lik[i] = log_diff_exp(neg_binomial_2_lcdf(10 | mu, phi), 
                                 neg_binomial_2_lcdf(5 | mu, phi));
    else
        log_lik[i] = neg_binomial_2_lccdf(10 | mu, phi);
  }
  
  Rho_i = multiply_lower_tri_self_transpose(L_Rho_i);
  
}
