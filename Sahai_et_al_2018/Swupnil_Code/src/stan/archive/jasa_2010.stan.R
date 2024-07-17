
data {
	int I;                                      //individuals
	int K;                                      //sub-populations
	int E;                                      //ego groups
	int A;                                      //alter groups
	real B[A,K];                                //known alter/subpopulation proportions
	vector<lower=0>[A] theta;                   //Dirichlet hyper-parameters - fixed
	int ego_index[I];                           //ego category variable
	int y[I,K];
	}

parameters {
	real alpha[I];                              //log-degree
	simplex[A] m[E];                            //mixing probability matrix
	vector<lower=0,upper=1>[K] inv_omega;       //inverse over-dispersion
	real mu_alpha;                              // prior mean for alpha
        real<lower=0> sigma_alpha;                  // prior scale for alpha
}

transformed parameters {
	real bmat[E,K];
	real x_a[A];
	for (e in 1:E) {
		for (k in 1:K) {
			for (a in 1:A) {
				x_a[a] <- m[e,a] * B[a,k];
			}
			bmat[e,k] <- sum(x_a); //* ((exp(- 7) / sum(x_a)) * exp(1 - exp(- 7) / sum(x_a))) ^ (.5); 
		}
	}
}

model {
	real mu_i_k;
	int ego;
	alpha ~ normal(mu_alpha, sigma_alpha);
	mu_alpha ~ normal(0 , 25);               // weakly informative 
        sigma_alpha ~ normal(0 , 5);             // weakly informative 
	for (e in 1:E) {
		m[e] ~ dirichlet(theta);
	}

	for (k in 1:K) {
		real omega_k_m1;
		omega_k_m1 <- inv(inv(inv_omega[k]) - 1);
		for (i in 1:I) {
			ego <- ego_index[i];
			mu_i_k <- omega_k_m1 * exp(alpha[i]) * bmat[ego, k];
			y[i,k] ~ neg_binomial(mu_i_k, omega_k_m1);
		}
	}
}

generated quantities {
	int ego;
	real mu_new[I,K];
	for (i in 1:I) {
		for (k in 1:K) {
			ego <- ego_index[i];
			mu_new[i, k] <- exp(alpha[i]) * bmat[ego, k];
		}
	}
}