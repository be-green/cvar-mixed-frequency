library(data.table)
library(readr)
library(readxl)
library(magrittr)
library(cmdstanr)
library(dplyr)

hpi <- read_excel("data/raw/HPI_PO_monthly_hist.xls",sheet =1)
hpi$Month <- as.POSIXct(hpi$Month,format = "%Y-%m-%d")
fixed <- read.csv("data/raw/MORTGAGE30US.csv")

fixed <- fixed %>%
  group_by(ym=paste(year(DATE), month(DATE))) %>%
  slice_head(n=1)
fixed <- fixed[fixed$DATE>=min(hpi$Month),]
fixed <- fixed[,c(1:2)]
fixed$DATE2 <- as.POSIXct(paste(year(fixed$DATE),month(fixed$DATE),"1",sep="-"),format = "%Y-%m-%d", tz="UTC")
fixed <- fixed[,-1]
names(fixed)[2] <- "Month"

hpi.dat <- merge(hpi,fixed,by = "Month")
names(hpi.dat)

hpi.dat$diffNA <- c(NA,diff(hpi.dat$"USA\n\n(NSA)",1))
hpi.dat$diffM30 <- c(NA,diff(hpi.dat$MORTGAGE30US))

x <- lag(hpi.dat$diffM30[!is.na(hpi.dat$diffM30)],2)
y <- hpi.dat$diffNA[!is.na(hpi.dat$diffNA)]
y <- y[!is.na(x)]
x <- x[!is.na(x)]

# int<lower=0> K;   // number of predictors
# int<lower=0> N; // number of observations
# matrix[N, K] X; // design matrix
# vector[N] y;      // outcome vector
# real tau;
k = 1
N = length(y)
x = scale(x, center = T)
x_scale = attr(x, "scaled:scale")
x_center = attr(x, "scaled:center")
x = as.numeric(x)
x = matrix(x, nrow = N)


y_scale = attr(x, "scaled:scale")
y_center = attr(x, "scaled:center")
y = as.numeric(y)

tau = 0.05

standata = list(
  K = k,
  N = N,
  X = x,
  y = y,
  tau = tau
)

qreg_model <- cmdstanr::cmdstan_model("src/Stan/qreg-no-missing.stan")

post = qreg_model$sample(data = standata, chains = 4, parallel_chains = 4)

library(tidybayes)
library(ggplot2)
draws = tidybayes::tidy_draws(post)

pdata = draws %>%
  as.data.table %>%
  melt(.) %>%
  .[variable %like% "yhat"] %>%
  .[, time := as.numeric(stringr::str_extract(variable, "[0-9]+"))] %>%
  .[,.(m = mean(value), s = sd(value)), by = time]


alpha =
  draws %>%
  as.data.table %>%
  melt(.) %>%
  .[variable %like% "alpha"] %>%
  .[, time := as.numeric(stringr::str_extract(variable, "[0-9]+"))] %>%
  .[,.(m = mean(value), s = sd(value)), by = time]
dt = hpi.dat$Month[4:nrow(hpi.dat)]
pdata[, date := dt]
alpha[, date := dt]

g = pdata %>%
  ggplot(aes(x = date, y = m,
             ymax = m + 1.96 * s,
             ymin = m - 1.96 * s)) +
  geom_pointinterval(alpha = 0.5, color = "orange")

g +
  geom_point(data = data.table(date = dt, m = y),
             aes(y = m, ymin = NULL, ymax = NULL),
             color = "black")

pdata %>%
  ggplot(aes(x = date, y = m)) +
  geom_point(alpha = 0.5, color = "black")+
  geom_point(data = alpha, color = "blue")
