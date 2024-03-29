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
  real mu_beta;
  real<lower=0.0000001> sigma_beta;
  vector[K] beta_raw;
  real<lower=0.0000001> sigma[K]; // variance of the state-space model
  vector[N_total_unknown] x_unknown; // unknown
}
transformed parameters {
  vector[K] beta;       // coefficients for predictors
  beta = mu_beta + sigma_beta * beta_raw;
  matrix[N, K] X;
  // we need to do this so that we can re-format the missing observations in
  // the potentially multivariate X matrix
  // complicated version of this:
  // https://mc-stan.org/docs/2_28/stan-users-guide/sliced-missing-data.html

  for (k in 1:K) {
    // assign to X the rows of the "observed index" with the known values
    // for column k
    // N_known is a vector which indexes these values
    X[segment(ii_obs, start_pos_known[k], N_known[k]), k] =
                    segment(x_known, start_pos_known[k], N_known[k]);

    if(N_total_unknown > 0) {
      // Same deal for N_unknown
      X[segment(ii_mis, start_pos_unknown[k], N_unknown[k]), k] =
              segment(x_unknown, start_pos_unknown[k], N_unknown[k]);
    }

  }
}
model {

  beta_raw ~ std_normal();
  // this is necessary because y may be centered but the quantile
  // will not be
  alpha ~ normal(quantile(y, tau), 2);
  mu_beta ~ normal(0, 2);
  sigma_beta ~ normal(0, 3);
  for(k in 1:K) {
    // centered at data scaling, informative for testing
    sigma[k] ~ normal(1,2);

    if(N_unknown[k] > 0) {
      // assumes scaled X variables
      // weakly informative prior in that case
      X[1, k] ~ normal(0, 2);
      X[2:N, k] ~ normal(X[1:(N - 1), k], sigma[k]);
    } else {
      sigma[k] ~ normal(0.01, 0.0001);
    }
  }
  y ~ skew_double_exponential(alpha + X * beta, 1, tau);
}
generated quantities {
  vector[N] y_pred;
  y_pred = alpha + X * beta;
}

