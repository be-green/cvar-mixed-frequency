# Build data from FRED and Shiller data

library(data.table)
library(readxl)
library(magrittr)
library(stringr)
library(quantmod)

if(!dir.exists("data/raw/yahoo")) {
  dir.create("data/raw/yahoo")
}

if(!file.exists("data/raw/yahoo/sp500.csv")) {
  yahoo_data <- getSymbols("^GSPC", from="1950-12-31", src='yahoo', auto.assign = F)
  yahoo_data <- as.data.table(yahoo_data)
  yahoo_data <- yahoo_data[,.(DATE = index, SP500 = `GSPC.Close`)]
  fwrite(yahoo_data, "data/raw/yahoo/sp500.csv")
  rm(yahoo_data)
}

eomonth <- function(date) {
  lubridate::ceiling_date(date, "month") - 1
}

fred_data <- lapply(list.files("data/raw/fred_source_data/", full.name = T,
                               pattern = "*.csv$"), fread)

gdp <- fred_data[[1]]$A939RX0Q048SBEA

fred_data[[1]]$A939RX0Q048SBEA <- c(NA, log(gdp[2:length(gdp)]) - log(gdp[1:(length(gdp) - 1)]))

cpi <- fred_data[[2]]$CPIAUCSL

fred_data[[2]]$CPIAUCSL <- c(NA, log(cpi[2:length(cpi)]) - log(cpi[1:(length(cpi) - 1)]))

macro_and_rates <- Reduce(function(x, y) merge(x, y, by = "DATE", all = T),fred_data)

sp500_data <- fread("data/raw/yahoo/sp500.csv")

cape_data <- read_xlsx("data/raw/shiller_data/cleaned_shiller_data.xlsx")
cape_data <- as.data.table(cape_data)
cape_data <- cape_data[,.(DATE = as.IDate(eomonth(as.Date(paste0(Year,"-", Month, "-25")))),
             CAPE = as.numeric(CAPE),
             TR_CAPE = as.numeric(TR_CAPE))][
  !is.na(CAPE)
]

all_data <- merge(sp500_data, macro_and_rates, by = "DATE", all = T) %>%
  merge(cape_data, by = "DATE", all = T)

all_data <- all_data[year(DATE) > 1951]
setorder(all_data, DATE)
all_data[, SP500 := nafill(SP500, type = "locf")]
all_data[, SP500 := nafill(SP500, type = "nocb")]

fwrite(all_data, "data/processed/time-series-data.csv")

