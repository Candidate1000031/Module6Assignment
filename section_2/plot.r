crude_oil <- read.table("./crudeoil.txt", header=T)
gasoline <- read.table("./gasoline.txt", header=F)

cr <- data.frame(date=ISOdate(crude_oil$Year, crude_oil$Mon, crude_oil$Day),
                 logpriceUS=log(crude_oil$US))
gs <- data.frame(date=ISOdate(crude_oil$Year, crude_oil$Mon, crude_oil$Day),
                 logprice=log(gasoline$V1))

plot(cr$date, cr$logpriceUS, type='l', xlab='Date', ylab='Log Price')
plot(cr$date, gs$logprice, type='l', xlab='Date', ylab='Log Price')
