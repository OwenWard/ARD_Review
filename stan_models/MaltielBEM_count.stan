functions {
  real ll_lpdf(real y, real d, real m, real rho){
    return lchoose(d, y) + lbeta(m * (1 / rho - 1) + y, d + (1 - m) * (1 / rho - 1) - y) - lbeta(m * (1 / rho - 1), (1 - m)*(1 / rho - 1));
  }
}

data {
  int N;
  int K;
  matrix[N, K] y;
  vector<lower=0, upper=1>[K] m;
  vector[N] L;
}

parameters {
  vector <lower=0>[N] d_raw;
  real <lower=3,upper=8> mu;
  real <lower=1./4,upper=2> sigma;
  vector <lower=0,upper=1>[K] rho;
}

transformed parameters {
  vector[N] d = L + d_raw;
}

model {
  d ~ lognormal(mu, sigma);
  mu ~ uniform(3, 8);
  sigma ~ uniform(1./4, 2);
  rho ~ uniform(0, 1);
  for (i in 1:N){
    for (k in 1:K) {
      if (y[i, k] >= 0) {
        y[i, k] ~ ll(d[i], m[k], rho[k]);
      }
    }
  }
}
