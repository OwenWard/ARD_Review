functions {

    real partial_sum(
        array[] int y_slice,
        int start,
        int end,
        array[] int ii,
        array[] int jj,
        real alpha,
        vector gamma,
        array[] real p,
        array[] real phi
    ) {
  
        // index
        int len = end - start + 1;
        
        // reindex vars
        array[len] int ii_slice = ii[start:end];
        array[len] int jj_slice = jj[start:end];
        
        // accumulator
        real psum = 0.0;
        
        for (n in 1:len) {
        
            real mu = exp(alpha + gamma[ii_slice[n]] + log(p[jj_slice[n]]));
            real ph = phi[jj_slice[n]];
        
            if (y_slice[n] == 0)
                psum += neg_binomial_2_lpmf(0 | mu, ph);
            else if (y_slice[n] == 1)
                psum += neg_binomial_2_lpmf(1 | mu, ph);
            else if (y_slice[n] == 2)
                psum += log_diff_exp(neg_binomial_2_lcdf(5 | mu, ph), 
                                     neg_binomial_2_lcdf(1 | mu, ph));
            else if (y_slice[n] == 3)
                psum += log_diff_exp(neg_binomial_2_lcdf(10 | mu, ph), 
                                     neg_binomial_2_lcdf(5 | mu, ph));
            else
                psum += neg_binomial_2_lccdf(10 | mu, ph);
          
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
    array[J] real<lower=0> phi;

}

transformed parameters {

    vector[I] gamma;
  
    gamma = gamma_star * sigma_gamma;
  
}

model {

    // population estimates priors
    for (j in 1:J) 
        p[j] ~ normal(pops[j], sds[j]);
  
    // gregariousness priors
    alpha ~ normal(0, 5);
    gamma_star ~ normal(0, 1);
    sigma_gamma ~ normal(0, 3);
    phi ~ normal(0, 3);
  
    // log-likelihood
    target += reduce_sum(partial_sum, y, grainsize, ii, jj, alpha, gamma, p, phi);
    
}

generated quantities {

    array[N] real log_lik;
    
    // log-likelihood
    for (i in 1:N) {
  
        real mu;
        real ph;
        
        mu = exp(alpha + gamma[ii[i]] + log(p[jj[i]]));
        ph = phi[jj[i]];
    
        if (y[i]==0)
            log_lik[i] = neg_binomial_2_lpmf(0 | mu, ph);
        else if (y[i]==1)
            log_lik[i] = neg_binomial_2_lpmf(1 | mu, ph);
        else if (y[i]==2)
            log_lik[i] = log_diff_exp(neg_binomial_2_lcdf(5 | mu, ph), 
                                      neg_binomial_2_lcdf(1|mu, ph));
        else if (y[i]==3)
            log_lik[i] = log_diff_exp(neg_binomial_2_lcdf(10 | mu, ph), 
                                      neg_binomial_2_lcdf(5 | mu, ph));
        else
            log_lik[i] = neg_binomial_2_lccdf(10 | mu, ph);
            
    }
  
} 
