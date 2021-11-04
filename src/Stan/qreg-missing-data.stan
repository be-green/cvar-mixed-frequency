data {
  // we need to re-shape a giant vector into a matrix
  // because Stan doesn't have ragged data structures
  int<lower=0> K;   // number of predictors
  int<lower=0> N; // number of observations
  real<lower=-1, upper = 1> tau; //target quantile
  int<lower=0> N_total_known; // total number of known observations in the whole matrix
  int<lower=0> N_total_unknown; // total number of unknown observations in the whole matrix

  int<lower=0> N_known[K];   // number of known observations by column of X
  int<lower=0> N_unknown[K]; // number of missing observations by column of X
  vector[N_total_known] x_known; // vector of known obs to be sliced
  vector[N_total_unknown] x_unknown; // vector of known obs to be sliced
  // index of known observations
  int<lower=1, upper=N_total_known + N_total_unknown> ii_obs[N_total_known];
  // index of missing observations
  int<lower=1, upper=N_total_known + N_total_unknown> ii_mis[N_total_unknown];
  vector[N] y;      // outcome vector
}
transformed data {
  matrix[N, K] X;
  // we need to do this so that we can re-format the missing observations in
  // the potentially multivariate X matrix
  // complicated version of this:
  // https://mc-stan.org/docs/2_28/stan-users-guide/sliced-missing-data.html
  int known_pos;
  known_pos = 1;

  int unknown_pos;
  unknown_pos = 1;
  for (k in 1:K) {
    // assign to X the rows of the "observed index" with the known values
    // for column k
    // N_known is a vector which indexes these values
    X[segment(ii_obs, known_pos, N_known[k]), k] = segment(x_known, known_pos, N_known[k]);

    if(N_total_unknown > 0) {
      // Same deal for N_unknown
      X[segment(ii_mis, unknown_pos, N_unknown[k]), k] = segment(x_unknown, unknown_pos, N_unknown[k]);
    }

    known_pos = known_pos + N_known[k];
    unknown_pos = unknown_pos + N_unknown[k];

  }
}
parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma[K]; // variance of innovations in local-level model
}
model {

  sigma ~ cauchy(0, 1);

  for(k in 1:K) {
    // assumes scaled X variables
    // weakly informative prior in that case
    if(N_unknown[k] > 0.5) {
      X[1, k] ~ normal(0, 1);
      X[N_unknown[k], k] ~ normal(X[N_unknown[k] - 1, k], sigma[k]);
    }
  }
  y ~ skew_double_exponential(alpha + X * beta, 1, tau);
}
generated quantities {
  matrix[N, K] X_return;
  X_return = X;
}

