library(cmdstanr)
library(data.table)
library(magrittr)

qreg_model <- cmdstanr::cmdstan_model("src/Stan/qreg-missing-data.stan")

X = matrix(rnorm(1000), ncol = 2)
beta = rnorm(ncol(X))
y = 0.5 + X %*% beta + rnorm(nrow(X))

# so stan doesn't get mad
y = as.vector(y)

# this is necessary because Stan is bad
# with ragged arrays...
missing_ids = lapply(data.frame(X), function(x_col) which(is.na(x_col)))
notmissing_ids = lapply(data.frame(X), function(x_col) which(!is.na(x_col)))

X_vector = c()
missing_ids_vector = c()
notmissing_ids_vector = c()
for(i in 1:ncol(X)) {
  X_vector = c(X_vector, X[,i])
  if(!is.null(missing_ids)) {
    missing_ids_vector = c(missing_ids_vector, missing_ids[[i]])
  }
  notmissing_ids_vector = c(notmissing_ids_vector, notmissing_ids[[i]])
}


X_missing = X_vector[which(is.na(X_vector))]
X_notmissing = X_vector[which(!is.na(X_vector))]

missing_by_column = colSums(is.na(X))
notmissing_by_column = colSums(!is.na(X))


standata = list(
  K = ncol(X),
  N = length(y),
  y = y,
  N_total_known = length(notmissing_ids_vector),
  N_total_unknown = length(missing_ids_vector),
  x_known = X_notmissing,
  x_unknown = X_missing,
  ii_obs = notmissing_ids_vector,
  ii_mis = missing_ids_vector,
  tau = 0.5,
  N_known = notmissing_by_column,
  N_unknown = missing_by_column,
  cumulative_known = cumsum(notmissing_by_column),
  cumulative_unknown = cumsum(missing_by_column)
)

options(mc.cores = 4)

test = qreg_model$sample(data = standata, chains = 4, iter_sampling = 1000, iter_warmup = 1000)

coef(quantreg::rq.fit.br(y = y, x = cbind(1,X), tau = 0.5))
test


