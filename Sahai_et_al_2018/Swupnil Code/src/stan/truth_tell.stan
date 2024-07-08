# stan model code to estimate global lambda for general response (i.e. how many X do you know?)
truth_code = '
  data{
    int<lower=0> N;
    int<lower=0> J;
    int<lower=0> d[N];
    int<lower=0> y[N,J];
    real<lower=0, upper=1> w;
    real<lower=0, upper=1> lambda[J];
    real<lower=0> p;
  }
  parameters{
    real<lower=0, upper=5> delta[J];
  }
  model{
    for(n in 1:N){
      for(j in 1:J){
        y[n,j] ~ poisson(pow(lambda[j]*w*d[n]*delta[j], p));
      }
    }
  }
'