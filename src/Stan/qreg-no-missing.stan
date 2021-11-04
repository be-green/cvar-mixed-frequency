data {
  // we need to re-shape a giant vector into a matrix
  // because Stan doesn't have ragged data structures
  int<lower=0> K;   // number of predictors
  int<lower=0> N; // number of observations
  matrix[N, K] X; // design matrix
  vector[N] y;      // outcome vector
  real tau;
}
parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors
}
model {

  y ~ skew_double_exponential(alpha + X * beta, 1, tau);
}

