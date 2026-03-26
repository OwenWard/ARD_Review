functions {

    real partial_sum(
        array[] int y_slice,
        int start,
        int end,
        array[] int ii,
        array[] int jj,
        real alpha,
        vector gamma,
        array[] real p
    ) {
  
        // no of units in slice
        int len = end - start + 1;
        
        // reindex vars
        array[len] int ii_slice = ii[start:end];
        array[len] int jj_slice = jj[start:end];

        // accumulator
        real psum = 0.0;
        
        for (n in 1:len) {
        
            real mu = exp(alpha + gamma[ii_slice[n]] + log(p[jj_slice[n]]));
            
            if (y_slice[n] == 0)
                psum += poisson_lpmf(0 | mu);
            else if (y_slice[n] == 1)
                psum += poisson_lpmf(1 | mu);
            else if (y_slice[n] == 2)
                psum += log_diff_exp(poisson_lcdf(5 | mu), poisson_lcdf(1 | mu));
            else if (y_slice[n] == 3)
                psum += log_diff_exp(poisson_lcdf(10 | mu), poisson_lcdf(5 | mu));
            else
                psum += poisson_lccdf(10 | mu);
    
        }
  
        return psum;
        
    }

}
data {
  
    int<lower=0> N;                 // No. of item-responses
    int<lower=0> I;                 // No. of Individuals.
    int<lower=0> J;                 // No. Items.
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

}

transformed parameters {

    vector[I] gamma;
  
    gamma = gamma_star * sigma_gamma;
  
}

model {

    // population estimates priors
    for (j in 1:J) 
        p[j] ~ normal(pops[j],sds[j]);
  
    // gregariousness priors
    alpha ~ normal(0, 5);
    gamma_star ~ normal(0, 1);
    sigma_gamma ~ normal(0, 3);

    // log-likelihood
    target += reduce_sum(partial_sum, y, grainsize, ii, jj, alpha, gamma, p);
    
}

generated quantities {

  array[N] real log_lik;
  
  // log-likelihood
  for (i in 1:N) {

    real mu = exp(alpha + gamma[ii[i]] + log(p[jj[i]]));

    if (y[i] == 0)
        log_lik[i] = poisson_lpmf(0|mu);
    else if (y[i] == 1)
        log_lik[i] = poisson_lpmf(1|mu);
    else if (y[i] == 2)
        log_lik[i] = log_diff_exp(poisson_lcdf(5 |mu), 
                                  poisson_lcdf(1 | mu));
    else if (y[i] == 3)
        log_lik[i] = log_diff_exp(poisson_lcdf(10 | mu), 
                                  poisson_lcdf(5 | mu));
    else
        log_lik[i] = poisson_lccdf(10|mu);
  }

}
