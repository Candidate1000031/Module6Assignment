crude_oil <- read.table("./crudeoil.txt", header=T)
gasoline <- read.table("./gasoline.txt", header=F)

crude.ts <- ts(crude_oil$US, start=c(1997, 1), end=c(2010, 39), frequency=52)
gasoline.ts <- ts(gasoline$V1, start=c(1997, 1), end=c(2010, 39), frequency=52)

# fit using ar with AIC
crude.ar <- ar(crude.ts, method='mle', aic=TRUE)
gasoline.ar <- ar(gasoline.ts, method='mle', aic=TRUE)

# Plot ACF and PACF to see if ACF decays fast and get an idea for the optimal
# order
acf(crude.ts)
pacf(crude.ts)
acf(gasoline.ts)
pacf(gasoline.ts)

print(crude.ar$order)
print(gasoline.ar$order)

crude.order <- crude.ar$order
gasoline.order <- gasoline.ar$order

# R doesn't allow me to run tsdiag on the model fitted using ar, so I'm
# performing another fit using arima, with the optimal order already known from
# above.
crude.ar_ <- arima(crude.ts, c(crude.order, 0, 0))
gasoline.ar_ <- arima(gasoline.ts, c(gasoline.order, 0, 0))

# Run tsdiag
tsdiag(crude.ar_)
tsdiag(gasoline.ar_)
