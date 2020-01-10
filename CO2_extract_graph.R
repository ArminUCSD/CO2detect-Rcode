# Extract graph from Peters 2017 - Fig 2

library(digitize)
file.name = 'Peters2017_Fig2_edited.png'

#---------------------------------------------------------------------
# Past time series
cal = ReadAndCal(file.name)
series1 = DigitData(col = 'red')
series2 = DigitData(col = 'blue')
series1.cal = Calibrate(series1, cal, 1959, 2017, 0, 30)
series2.cal = Calibrate(series2, cal, 1959, 2017, 0, 30)

# Fix time axis
t = 1959:2017
observations = series1.cal$y
reconstructed = series2.cal$y
save(t, observations, reconstructed, file = 'Peters2017_Fig2.Rdata')

plot(t, reconstructed, type="l", col='blue')
lines(t, observations, col='red')

imb = observations - reconstructed
mean(imb)
sd(imb)

# Save as text file
write.table(data.frame(t, observations, reconstructed), row.names = F, file = 'Peters2017_Fig2_past.txt')


#---------------------------------------------------------------------
# Future scenarios
cal = ReadAndCal(file.name)
series1 = DigitData(col = 'blue')
series2 = DigitData(col = 'red')
series3 = DigitData(col = 'green')
series1.cal = Calibrate(series1, cal, 2017, 2040, 0, 30)
series2.cal = Calibrate(series2, cal, 2017, 2040, 0, 30)
series3.cal = Calibrate(series3, cal, 2017, 2040, 0, 30)
series1.cal[1] = series2.cal[1] = series3.cal[1] = (series1.cal[1] + series2.cal[1] + series3.cal[1])/3
future = list(t = c(2017, 2020, 2030, 2040), '1'=series1.cal$y, '0'=series2.cal$y,
                 '-1'=series3.cal$y)

# series1 = digitize(file.name)

save(future, file = 'Peters2017_Fig2_future.Rdata')


# Fit polynomial
load(file = 'Peters2017_Fig2_future.Rdata')
t.in = future$t
y.in = cbind(future$'1', future$'0', future$`-1`)

# Fit functional form
fit = lm(y.in ~ t.in + I(t.in^2))

t = 2017:2040
tt = cbind(rep(1, length(t)), t, t^2)
y = tt %*% fit$coefficients

# Check fit
matplot(t.in, y.in)
matlines(t, y)

tmp = data.frame(t, y); names(tmp) = c('t', '1%', '0%', '-1%')
write.table(tmp, row.names = F, file = 'Peters2017_Fig2_future.txt')

# Extrapolate to longer time frame
t = 2017:2050
tt = cbind(rep(1, length(t)), t, t^2)
y = tt %*% fit$coefficients
matplot(t, y, type='l')
matpoints(t.in, y.in)
tmp = data.frame(t, y); names(tmp) = c('t', '1%', '0%', '-1%')
write.table(tmp, row.names = F, file = 'Peters2017_Fig2_future2050.txt')
