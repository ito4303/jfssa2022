data {
  int<lower = 0> M;
  int<lower = 0> J;
  array[M, J] int<lower = 0, upper = 1> Y;
  vector[M] Elev;
  vector[M] Forest;
  matrix[M, J] Wind;
}

parameters {
  array[6] real beta; 
}

model {
  for (m in 1:M) {
    real psi = inv_logit(beta[1] + beta[2] * Elev[m] + beta[3] * Forest[m]);
    vector[J] p = inv_logit(beta[4] + beta[5] * Elev[m] + beta[6] * Wind[m]');
    if (sum(Y[m]) == 0)
      target += log_sum_exp(bernoulli_lpmf(0 | psi),
                            bernoulli_lpmf(1 | psi) + bernoulli_lpmf(0 | p));
    else
      target += bernoulli_lpmf(1 | psi) + bernoulli_lpmf(Y[m] | p);
  }
}
