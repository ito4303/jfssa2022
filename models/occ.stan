data {
  int<lower=0> M;
  int<lower=0> J;
  array[M, J] int<lower=0, upper=1> Y;
  vector[M] Elev;
  vector[M] Forest;
  matrix[M, J] Wind;
}

parameters {
  array[6] real beta; 
}

transformed parameters {
  vector<lower=0, upper=1>[M] psi;
  matrix<lower=0, upper=1>[M, J] p;
  psi = inv_logit(beta[1] + beta[2] * Elev + beta[3] * Forest);
  for (m in 1:M)
    p[m] = inv_logit(beta[4] + beta[5] * Elev[m] + beta[6] * Wind[m]);
}

model {
  for (m in 1:M) {
    if (sum(Y[m]) == 0) // not detected
      target += log_sum_exp(bernoulli_lpmf(0 | psi[m]),
                            bernoulli_lpmf(1 | psi[m])
                            + bernoulli_lpmf(0 | p[m]));
    else                // detected
      target += bernoulli_lpmf(1 | psi[m])
                + bernoulli_lpmf(Y[m] | p[m]);
  }
}

generated quantities {
  array[M] int<lower=0, upper=1> z;
  int<lower=0, upper=M> Nocc;
  for (m in 1:M)
    if (sum(Y[m]) > 0) {  // detected
      z[m] = 1;
    } else {              // not detected
      real lp1 = bernoulli_lpmf(0 | psi[m]);
      real lp2 = bernoulli_lpmf(1 | psi[m])
                 + bernoulli_lpmf(Y[m] | p[m]);
      z[m] = bernoulli_rng(exp(lp2) / (exp(lp1) + exp(lp2)));
    }
  Nocc = sum(z);
}
