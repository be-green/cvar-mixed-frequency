library(cmdstanr)
library(data.table)
library(magrittr)

qreg_model <- cmdstanr::cmdstan_model("src/Stan/qreg-missing-data.stan")

X = matrix(rnorm(1000), ncol = 2)
beta = rnorm(ncol(X))
y = 0.5 + X %*% beta + rnorm(nrow(X))

# so stan doesn't get mad
y = as.vector(y)

X_vector = c()
for(i in 1:ncol(X)) {
  X_vector = c(X_vector, X[,i])
}

# this is necessary because Stan is bad
# with ragged arrays...
missing_ids = as.vector(apply(X, MARGIN = 2, function(x_col) which(is.na(x_col))))
notmissing_ids = as.vector(apply(X, MARGIN = 2, function(x_col) which(!is.na(x_col))))
X_missing = X_vector[missing_ids]
X_notmissing = X_vector[notmissing_ids]

missing_by_column = colSums(is.na(X))
notmissing_by_column = colSums(!is.na(X))


standata = list(
  K = ncol(X),
  N = length(y),
  y = y,
  N_total_known = length(notmissing_ids),
  N_total_unknown = length(missing_ids),
  x_known = X_notmissing,
  x_unknown = X_missing,
  ii_obs = notmissing_ids,
  ii_mis = missing_ids,
  tau = 0.5,
  N_known = notmissing_by_column,
  N_unknown = missing_by_column
)

options(mc.cores = 4)

test = qreg_model$sample(data = standata, chains = 4, iter_sampling = 1000, iter_warmup = 1000)


qreg <- cmdstan_model("src/stan/qreg-no-missing.stan")

qreg$sample(list(y = y, X = X, tau = 0.5, N = length(y), K = ncol(X)))

coef(quantreg::rq.fit.br(y = y, x = cbind(1,X), tau = 0.5))



