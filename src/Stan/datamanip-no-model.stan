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
  int start_pos_known[K];
  int start_pos_unknown[K];

  vector[N_total_known] x_known; // vector of known obs to be sliced
  // index of known observations
  int<lower=1, upper=N_total_known + N_total_unknown> ii_obs[N_total_known];
  // index of missing observations
  int<lower=1, upper=N_total_known + N_total_unknown> ii_mis[N_total_unknown];
  vector[N] y;      // outcome vector
}
parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors

  // https://mc-stan.org/docs/2_28/reference-manual/vector-and-matrix-data-types.html#cholesky-factors-of-covariance-matrices
  cholesky_factor_corr[K] L_Omega;

  // VAR(1) coefs
  row_vector[K] alpha_ar;
  matrix<lower=-1, upper=1>[K, K] beta_ar;

  vector<lower=0.0001>[K] sigma; // variance of the state-space model
  vector[N_total_unknown] x_unknown; // unknown
}
transformed parameters {
  matrix[N, K] X;
  array[N - 1] row_vector[K] var_mu;
  array[N - 1] row_vector[K] X_array;
  // we need to do this so that we can re-format the missing observations in
  // the potentially multivariate X matrix
  // complicated version of this:
  // https://mc-stan.org/docs/2_28/stan-users-guide/sliced-missing-data.html

  for (k in 1:K) {
    // assign to X the rows of the "observed index" with the known values
    // for column k
    // N_known is a vector which indexes these values
    X[segment(ii_obs, start_pos_known[k], N_known[k]), k] = segment(x_known, start_pos_known[k], N_known[k]);

    if(N_total_unknown > 0) {
      // Same deal for N_unknown
      X[segment(ii_mis, start_pos_unknown[k], N_unknown[k]), k] = segment(x_unknown, start_pos_unknown[k], N_unknown[k]);
    }
  }

  for (n in 1:(N - 1)) {
    var_mu[n] = alpha_ar + X[n, ] * beta_ar;
  }

  for(n in 1:(N - 1)) {
    X_array[n] = X[n + 1,];
  }
}
model {
  for(n in 1:N) {
    to_vector(X[n,]) ~ std_normal();
  }
  alpha_ar ~ std_normal();
  for(k in 1:K) {
    to_vector(beta_ar[k,]) ~ std_normal();
  }
  beta ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(1);
  sigma ~ std_normal();
  y ~ std_normal();
}
generated quantities {
  matrix[K, K] Sigma;
  Sigma = diag_pre_multiply(sigma, L_Omega);
  vector[N] y_pred;
  y_pred = alpha + X * beta;
}


