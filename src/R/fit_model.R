library(cmdstanr)
library(data.table)
library(magrittr)

qreg_model <- cmdstanr::cmdstan_model("src/Stan/qreg-missing-data-var.stan")

data <- fread("data/processed/time-series-data.csv")
data <- data[DATE > as.Date("1990-01-01")]
data[,TR_CAPE := NULL]
data[,TB3SMFFM := NULL]
shift = 1

X = as.matrix(data[,lapply(.SD, function(x) scale(as.numeric(x))), .SDcols = !c("DATE","SP500")])
y = data$SP500
y = log(y[(1 + shift):length(y)]) - log(y[1:(length(y) - shift)])
X = X[1:(nrow(X) - shift),]

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
  tau = 0.05,
  N_known = notmissing_by_column,
  N_unknown = missing_by_column,
  start_pos_known = c(1, cumsum(notmissing_by_column) + 1)[1:ncol(X)],
  start_pos_unknown = c(1, cumsum(missing_by_column) + 1)[1:ncol(X)]
)

options(mc.cores = 4)

test = qreg_model$variational(data = standata)

library(tidybayes)
library(ggplot2)
draws <- tidy_draws(test)

draws <- draws[, colnames(draws) %like% "y_pred|\\."]
draws <- as.data.table(draws)
draws <- melt(draws, c(".chain", ".iteration", ".draw"), variable.name = "time_period",
     value.name = "pred_CVaR")

draws[, time_period := as.integer(str_extract(time_period, "[0-9]+"))]

plot_data <- draws %>%
  group_by(time_period) %>%
  median_qi(.width = c(0.1, 0.25, 0.75, 0.8, 0.9, 0.95))



plot_data %>%
  merge(data.frame(time_period = 1:nrow(data), Date = data$DATE)) %>%
  ggplot(aes(x = Date, y = pred_CVaR, ymin = .lower,
             ymax = .upper, fill = factor(.width, levels = sort(unique(.width), decreasing = T)))) +
  geom_ribbon() +
  scale_fill_brewer() +
  labs(fill = "Interval",
         y = "Predicted 5th Percentile Return")



