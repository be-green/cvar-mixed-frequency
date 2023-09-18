data {
  int<lower=0> K;   // number of predictors
  int<lower=0> N; // number of observations
  matrix[N, K] X; // design matrix
  vector[N] y;      // outcome vector
  real tau;
}
parameters {
  vector[N] alpha;      // intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0.0001> sigma; // AR(1) variance for alpha
}
model {
  beta ~ normal(0, 4);
  sigma ~ student_t(4, 0, 1);

  alpha[1] ~ normal(inv_Phi(tau), sigma);
  for (n in 2:N)
    alpha[n] ~ normal(alpha[n-1], sigma);

  y ~ skew_double_exponential(alpha + X * beta, 1, tau);
}
generated quantities {
    vector[N] yhat = alpha + X * beta;
}

