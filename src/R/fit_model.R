# set working directory to be the project base folder
#setwd("./Documents/GitHub/cvar-mixed-frequency")

# Install and Load relevant packages
#####

if(!require(cmdstanr)) {
  install.packages("cmdstanr",
                   repos = c("https://mc-stan.org/r-packages/",
                             getOption("repos")))
  install_cmdstan(cores = 2) # this natively sets the path
}

if(!require(data.table)) install.packages("data.table")
if(!require(magrittr)) install.packages("magrittr")
if(!require(dplyr)) install.packages("dplyr")
if(!require(tidybayes)) install.packages("tidybayes")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggplot2)) install.packages("ggthemes")
if(!require(future.apply)) install.packages("future.apply")

library(future.apply)
library(ggthemes)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(cmdstanr)
library(data.table)
library(magrittr)

#####
# Build Model
#####
qreg_model <- cmdstanr::cmdstan_model("src/Stan/qreg-missing-data-rw.stan")

#####
# Load in Model Data
#####
data <- fread("data/processed/time-series-data.csv")
data <- data[DATE > as.Date("2017-01-01")]
shift = 1
X_oos = as.matrix(data[DATE >= as.Date("2020-01-01"),lapply(.SD, function(x) scale(as.numeric(x))), .SDcols = !c("DATE","SP500")])

X = as.matrix(data[DATE < as.Date("2020-01-02"),lapply(.SD, function(x) scale(as.numeric(x))), .SDcols = !c("DATE","SP500")])
y = data[DATE < as.Date("2020-01-02")]$SP500
y = log(y[(1 + shift):length(y)]) - log(y[1:(length(y) - shift)])
X = X[1:(nrow(X) - shift),]

y_oos = data[DATE < as.Date("2020-01-01")]$SP500
y_oos = log(y_oos[(1 + shift):length(y)]) - log(y_oos[1:(length(y_oos) - shift)])

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

y_scaled <- (y - mean(y)) / sd(y)

standata = list(
  K = ncol(X),
  N = length(y),
  y = y_scaled,
  N_total_known = length(notmissing_ids_vector),
  N_total_unknown = length(missing_ids_vector),
  x_known = X_notmissing,
  ii_obs = notmissing_ids_vector,
  ii_mis = missing_ids_vector,
  tau = 0.5,
  N_known = notmissing_by_column,
  N_unknown = missing_by_column,
  start_pos_known = c(1, cumsum(notmissing_by_column) + 1)[1:ncol(X)],
  start_pos_unknown = c(1, cumsum(missing_by_column) + 1)[1:ncol(X)]
)

options(mc.cores = 4)

# variational optimizes and objective assuming that the posterior
# is a mixture of gaussians
# much faster than MCMC, but can be less accurate
# and more fragile since it's a point estimate (essentially)
test = qreg_model$sample(data = standata,iter_warmup = 400, iter_sampling = 400,
                         refresh = 10, max_treedepth = 10)

tmp = test$draws()

draws <- tidy_draws(test)

test_X <- draws[, colnames(draws) %like% "X\\["]
test_X <- test_X[1,]
nms <- colnames(test_X)

reformat_X <- matrix(ncol = ncol(X), nrow = nrow(X))
for(i in 1:ncol(test_X)) {
   eval(parse(text = paste0("reformat_", nms[i], " <- ", test_X[1,i])))
}

reconstruct_x <- draws %>%
  select(starts_with("X[")) %>%
  .[1,] %>%
  format_matrix
max(reconstruct_x - X, na.rm = T)

# @param draw_row a single draw for the parameter vector
# @param K dimension of the X matrix
parse_draw <- function(new_X, draw_row, K, N) {

  intercept <- rep(NA, K)
  coefs <- matrix(nrow = K, ncol = K)
  covmat <- matrix(nrow = K, ncol = K)
  last_X <- matrix(nrow = 1, ncol = K)
  beta <- c()

  for(i in 1:K) {
    intercept[i] <- draw_row[[paste0("alpha_ar[",i,"]")]]
    beta[i] <- draw_row[[paste0("beta[",i,"]")]]
    for(j in 1:K) {
      coefs[i, j] <- draw_row[[paste0("beta_ar[",i,",", j, "]")]]
      covmat[i, j] <- draw_row[[paste0("Sigma[",i,",", j, "]")]]
    }
    last_X[1,i] <- draw_row[[paste0("X[",N,",", i, "]")]]
  }
  covmat <- crossprod(covmat)

  alpha <- draw_row[["alpha"]]

  pred <- MASS::mvrnorm(n = 1, mu = as.vector(intercept + last_X %*% coefs), Sigma = covmat)

  for(i in 1:(nrow(new_X))) {
    missing <- which(is.na(new_X[i,]))
    new_X[i, missing] <- pred[missing]
    pred = MASS::mvrnorm(n = 1, mu = as.vector(intercept + new_X[i,] %*% coefs), Sigma = covmat)
  }
  alpha + new_X %*% beta
}

plan(multisession(workers = 8))
#
# # matrix of predicted values
# pred_y <- future_apply(draws, 1, function(x) parse_draw(X_oos, x,K =  ncol(X), N = nrow(X)), future.seed=TRUE)

in_sample <- draws[, colnames(draws) %like% "y_pred|\\."]
in_sample <- as.data.table(in_sample)
in_sample <- melt(in_sample, c(".chain", ".iteration", ".draw"), variable.name = "time_period",
     value.name = "pred_CVaR")

in_sample[, time_period := as.integer(str_extract(time_period, "[0-9]+"))]

plot_data <- in_sample %>%
  group_by(time_period) %>%
  median_qi(.width = c(0.1, 0.25, 0.75, 0.8, 0.9, 0.95))
#
# oos <- data.table(time_period = max(plot_data$time_period) + 1:498,
#                   Date = tail(data$DATE, 498),
#                   pred_y) %>%
#   melt(1:2, value.name = "pred_CVaR") %>%
#   .[,time_period := NULL] %>%
#   .[,variable := NULL] %>%
#   group_by(Date) %>%
#   median_qi(.width = c(0.1, 0.25, 0.75, 0.8, 0.9, 0.95)) %>%
#   as.data.table
#
# oos[, Sample := "Out of Sample"]

plot_data %>%
  merge(data.frame(time_period = 1:nrow(X), Date = data$DATE[1:nrow(X)],
                   y_scaled)) %>%
  as.data.table %>%
  .[,time_period := NULL] %>%
  .[,Sample := "In Sample"] %>%
  # rbind(oos) %>%
  ggplot(aes(x = Date, y = pred_CVaR, ymin = .lower,
             ymax = .upper, fill = factor(.width, levels = sort(unique(.width), decreasing = T)))) +
  geom_ribbon() +
  geom_point(size = 0.01, color = "white") +
  scale_fill_brewer() +
  labs(fill = "Interval",
         y = "Predicted 5th Percentile Return") +
  theme_clean() +
  geom_point(aes(y = y_scaled))
