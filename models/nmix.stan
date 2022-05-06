data {
  int<lower = 0> M;
  int<lower = 0> J;
  array[M, J] int<lower = 0> Y;
  vector[M] Xsite;
  matrix[M, J] Xsurvey;
  int<lower = 0> K;  // upper limit of N for calculation
}

parameters {
  array[4] real beta;
}

transformed parameters {
  matrix[M, K + 1] lp;
  for (m in 1:M) {
    real log_lambda = beta[1] + beta[2] * Xsite[m];
    vector[J] logit_p = beta[3] + beta[4] * Xsurvey[m]';
    int y_max = max(Y[m]);
    lp[m, 1:y_max] = rep_row_vector(negative_infinity(), y_max);
    for (k in y_max:K)
      lp[m, k + 1] = poisson_log_lpmf(k | log_lambda)
                     + binomial_logit_lpmf(Y[m] | k, logit_p);
  }
}

model {
  for (m in 1:M)
    target += log_sum_exp(lp[m]);
}

generated quantities {
  array[M] int N;
  int Ntotal;
  for (m in 1:M)
    N[m] = categorical_rng(softmax(lp[m]')) - 1;
  Ntotal = sum(N);
}
