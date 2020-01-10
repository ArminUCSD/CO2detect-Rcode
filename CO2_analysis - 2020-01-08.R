# Analysis code for: ----
# A statistical protocol for atmospheric verification of CO2 emissions
# Armin Schwartzman and Ralph Keeling

#===================================================================================
# Load data ----
library(car)

past = read.table(file = 'Peters2017_Fig2_past.txt', header = T)
names(past) = c('t', 'obs', 'rec')
past$imb = past$obs - past$rec

future = read.table(file = 'Peters2017_Fig2_future.txt', header = T, col.names = c('t', '1%', '0%', '-1%'))
names(future) = c('t', '1%', '0%', '-1%')

#===================================================================================
# Reproduce Peters' 2017 Fig. 2 ----

# Error bands around 1%
plot(NULL, xlim=c(1959,2040), ylim=c(0,30),
     xlab='', ylab='Atmospheric growth rate (GtCO2/year)')
lines(past$t, past$obs, type="l", col='black', lwd=2)
lines(past$t, past$rec, col='gray', lwd=2)
mu = future$`1%`
sigma = 3
q = qnorm(0.025)
lines(future$t, mu, col='orange', lty=1, lwd=2)
lines(future$t, mu + sigma, col='black', lty=2, lwd=2)
lines(future$t, mu - sigma, col='black', lty=2, lwd=2)
lines(future$t, mu + q*sigma, col='black', lty=3, lwd=2)
lines(future$t, mu - q*sigma, col='black', lty=3, lwd=2)
mu = future$`0%`
lines(future$t, mu, col='green', lty=1, lwd=2)
mu = future$`-1%`
lines(future$t, mu, col='cyan', lty=1, lwd=2)


#===================================================================================
# The imbalance process ----
# AR(1) process with coefficient 0.48

plot(past$t, past$imb, type='l', lwd=2, col='blue', xlab='', ylab = 'imbalance', cex=cex)
abline(0, 0)

# Basic statistics
mean(past$imb)
sd(past$imb)
sigma = 3

n = length(past$imb)
plot(past$imb[1:(n-1)], past$imb[2:n])
m.ar = lm(past$imb[2:n] ~ past$imb[1:(n-1)])
summary(m.ar)
abline(m.ar, col='blue')

# AR modeling
acf(past$imb)
pacf(past$imb)
mle = ar.mle(past$imb, aic = TRUE, order.max = NULL, demean = FALSE, intercept = FALSE)
mle$aic
mle
sqrt(mle$asy.var.coef)
sqrt(mle$var.pred)

ols = ar.ols(past$imb, aic = FALSE, order.max = 1, demean = FALSE, intercept = FALSE)
ols
acf(ols$resid[2:n])
qqPlot(ols$resid[2:n])

acf(past$imb)
k = 0:17
lines(k, mle$ar^k, col='red')
pacf(past$imb)

mean(past$imb)
sd(past$imb)/sqrt(n*(1 - ols$ar[1]^2))

# Fixed values
sigma = 3
rho = 0.5


#===================================================================================
# Simulate future imbalance processes ----

# Future needs to be extended to allow for longer detection delays
future = read.table(file = 'Peters2017_Fig2_future2050.txt', header = T, col.names = c('t', '1%', '0%', '-1%'))
names(future) = c('t', '1%', '0%', '-1%')
N = nrow(future)

# Simulation
sigma = 3
rho.est = 0.48
Nsim = 1e5

gen.AR = function(N, rho) {
  # Generate AR process with AR parameter rho and marginal variance 1
  eps = rnorm(N)
  x = eps
  for(i in 2:N) x[i] = rho*x[i-1] + sqrt(1-rho^2)*eps[i]
  return(x)
}

# Simulated imbalance (AR process with parameter rho)
# imb = array(0, dim = c(N, Nsim))
# for(i in 1:Nsim) imb[,i] = gen.AR(N, rho.est)
# save(imb, file = 'imb.Rdata')
load(file = 'imb.Rdata')

# Integrated imbalance (IAR process with AR parameter rho, variance v)
Imb = apply(imb, 2, cumsum)
Imb.var = apply(Imb, 1, var)
v = ((1:N)*(1-rho.est^2) - 2*rho.est*(1-rho.est^(1:N)))/(1-rho.est)^2  # Compare formula
# Formula from
# Hunter, A. J., & Connors, W. A. Statistics of an autoregressive correlated random walk along a return path. Electronics Letters, 53(23), 1550-1552 (2017).


#===================================================================================
# Detection delay probability functions ----

#-----------------------------------------------------------------------------------
# Detection time functions ----

# Probability that detection time is smaller than t
PT = function(x.array, thresh) {
  detect = x.array
  for(j in 1:ncol(detect))
    # for(i in 1:N) detect[i,j] = min(x.array[1:i,j]) < thresh
    detect[,j] = cumsum(x.array[,j] < thresh) > 0
  return(apply(detect, 1, mean))
}

# Create a list of detection time CDFs according to a fixed scenario and BAU
PT.list = function(scenario, z, mu, thr) {
  PT.0 = array(NA, dim = c(N, length(thr)))
  PT.1 = array(NA, dim = c(N, length(thr)))
  for(i in 1:length(thr)) {
    PT.0[,i] = PT(z, thr[i])
    PT.1[,i] = PT(z + mu, thr[i])
  }
  out = list(t = future$t, PT.0, PT.1)
  names(out) = c('t', '1%', scenario)
  return(out)
}

# Plot CDFs according to scenarios
prob.plot = function(p.frame, ymax = 1, ylab = '', main = '', probs = c(0.05, 0.5, 0.95),
                     mar = c(2.2, 4.1, 1.1, 1.1), legendloc = 'topleft') {
  par(mar = mar)
  plot(NULL, xlim=future$t[c(1,N)], ylim=c(0,ymax), xlab='', ylab=ylab, main=main)
  for(j in 1:length(probs)) abline(h=probs[j], col='gray', lty=1)
  lines(future$t, p.frame$'1%', col='orange', lty=1, lwd=2)
  if(!is.null(p.frame$'0%')) {
    lines(future$t, p.frame$'0%', col='green', lty=1, lwd=2)
    legend(legendloc, legend=c('R0','BAU'), col=c('green','orange'), lwd=2)
  }
  else if(!is.null(p.frame$'-1%')) {
    lines(future$t, p.frame$'-1%', col='cyan', lty=1, lwd=2)
    legend(legendloc, legend=c('R1','BAU'), col=c('cyan','orange'), lwd=2)
  }
}

# Time differencing
diff2 = function(x) {
  n = length(x)
  if(n > 3) dx = c(x[2]-x[1], diff(x, 2)/2, x[n]-x[n-1])
  else dx = x[2]-x[1]
}


#-----------------------------------------------------------------------------------
# Evaluate Peters et al. ----
# Cumulative detection probability

criterion = 'AGR'
sigma = 3
z = imb*sigma
thr = -2*sigma

# Green scenario (0%)
scenario = '0%'
mu = future[[scenario]] - future$`1%`
imb.PT.0.2sigma = PT.list(scenario, z, mu, thr)

prob.plot(imb.PT.0.2sigma, ylab='Detection probability', main=paste(criterion, 'sigma =', sigma, 'thr =', thr))
abline(v = future$t[which.min(abs(imb.PT.0.2sigma[[scenario]] - 0.95))], lty=2)

# Cyan scenario (-1%)
scenario = '-1%'
mu = future[[scenario]] - future$`1%`
imb.PT.1.2sigma = PT.list(scenario, z, mu, thr)

prob.plot(imb.PT.1.2sigma, ylab='Detection probability', main=paste(criterion, 'sigma =', sigma, 'thr =', thr))
abline(v = future$t[which.min(abs(imb.PT.1.2sigma[[scenario]] - 0.95))], lty=2)

# save(file = 'results.2sigma.RData', list = c('imb.PT.0.2sigma', 'imb.PT.1.2sigma'))
load(file = 'results.2sigma.RData')

#===================================================================================
# Estimate detection delays ----

#-----------------------------------------------------------------------------------
# Calibrate threshold (AGR)

criterion = 'AGR'
scenarios = c('0%', '-1%')
sigmas = c(3, 1.5, 0.75)

# Initialize
results.AGR = data.frame(array(NA, dim=c(length(scenarios)*length(sigmas), 7)))
colnames(results.AGR) = c('scenario', 'sigma', 'thr', 'delay.mode', 'delay.05', 'delay.50', 'delay.95')
imb.PT.AGR = array(data.frame(NULL), c(length(scenarios),length(sigmas)))

for(k in 1:length(scenarios)) {
  for(j in 1:length(sigmas)) {
    # BAU scenario (1%)
    z = imb*sigmas[j]
    z.min = apply(z, 2, cummin)
    q.AGR.BAU = apply(z.min, 1, function(x) quantile(x, probs = 0.05))
  
    # Alternative scenario (0% or -1%)
    mu = future[[scenarios[k]]] - future$`1%`
    z.0.min = apply(z + mu, 2, cummin)
    q.AGR = apply(z.0.min, 1, function(x) quantile(x, probs = 0.95))
    thr = q.AGR.BAU[which.min(abs(q.AGR.BAU - q.AGR))]
    
    # Compute probability distribution
    imb.PT = PT.list(scenarios[k], z, mu, thr)

    # Collect results
    imb.PT.AGR[[k,j]] = imb.PT
    ind = (k-1)*length(sigmas)+j  # row index
    results.AGR$scenario[ind] = scenarios[k]
    results.AGR$sigma[ind] = sigmas[j]
    results.AGR$thr[ind] = round(thr, 2)
    results.AGR$delay.mode[ind] = which.max(diff2(imb.PT[[scenarios[k]]])) - 1
    results.AGR$delay.05[ind] = which.min(abs(imb.PT[[scenarios[k]]] - 0.05)) - 1
    results.AGR$delay.50[ind] = which.min(abs(imb.PT[[scenarios[k]]] - 0.5)) - 1
    results.AGR$delay.95[ind] = which.min(abs(imb.PT[[scenarios[k]]] - 0.95)) - 1
  }
}
# save(file = 'results.AGR.RData', list = c('results.AGR', 'imb.PT.AGR'))
load(file = 'results.AGR.RData')


#-----------------------------------------------------------------------------------
# Calibrate threshold (ACC)

criterion = 'ACC'
scenarios = c('0%', '-1%')
sigmas = c(3, 1.5, 0.75)

# Initialize
results.ACC = data.frame(array(NA, dim=c(length(scenarios)*length(sigmas), 7)))
colnames(results.ACC) = c('scenario', 'sigma', 'thr', 'delay.mode', 'delay.05', 'delay.50', 'delay.95')
Imb.PT.ACC = array(data.frame(NULL), c(length(scenarios),length(sigmas)))

for(k in 1:length(scenarios)) {
  for(j in 1:length(sigmas)) {
    # BAU scenario (1%)
    z = Imb*sigmas[j]
    z.min = apply(z, 2, cummin)
    q.ACC.BAU = apply(z.min, 1, function(x) quantile(x, probs = 0.05))
    
    # Alternative scenario (0% or -1%)
    mu = cumsum(future[[scenarios[k]]] - future$`1%`)
    z.0.min = apply(z + mu, 2, cummin)
    q.ACC = apply(z.0.min, 1, function(x) quantile(x, probs = 0.95))
    thr = q.ACC.BAU[which.min(abs(q.ACC.BAU - q.ACC))]
    
    # Compute probability distribution
    Imb.PT = PT.list(scenarios[k], z, mu, thr)
    
    # Collect results
    Imb.PT.ACC[[k,j]] = Imb.PT
    ind = (k-1)*length(sigmas)+j  # row index
    results.ACC$scenario[ind] = scenarios[k]
    results.ACC$sigma[ind] = sigmas[j]
    results.ACC$thr[ind] = round(thr, 2)
    results.ACC$delay.mode[ind] = which.max(diff2(Imb.PT[[scenarios[k]]])) - 1
    results.ACC$delay.05[ind] = which.min(abs(Imb.PT[[scenarios[k]]] - 0.05)) - 1
    results.ACC$delay.50[ind] = which.min(abs(Imb.PT[[scenarios[k]]] - 0.5)) - 1
    results.ACC$delay.95[ind] = which.min(abs(Imb.PT[[scenarios[k]]] - 0.95)) - 1
  }
}
# save(file = 'results.ACC.RData', list = c('results.ACC', 'Imb.PT.ACC'))
load(file = 'results.ACC.RData')


#===================================================================================
# Figures

# Figure parameters
figpath = '../figures'
 res = 300
#res = 1200
wd = 360 * res/72
ht = 300 * res/72
cex = 1.2

#-----------------------------------------------------------------------------------
# The imbalance process

png(filename = "../figures/imbalance_ts.png", width = wd, height = ht, res = res)
par(cex = cex)
par(mar = c(3, 4, 1.1, 1.1))
plot(past$t, past$imb, type='l', lwd=2, col='blue', xlab='', ylab = 'imbalance', cex=cex)
abline(0, 0)
dev.off()

# Plots
png(filename = paste(figpath, "/imbalance_qqplot.png", sep=''), width = wd, height = ht, res = res)
par(cex = cex)
par(mar = c(3, 4, 1.1, 1.1))
qqPlot(past$imb, xlab='normal quantlies', ylab = 'sample quantiles')
dev.off()

png(filename = paste(figpath, "/imbalance_acf.png", sep=''), width = wd, height = ht, res = res)
par(cex = cex)
par(mar = c(5, 4, 1.1, 1.1))
acf(past$imb, main='')
# k = 0:17; lines(k, mle$ar^k, col='red')
dev.off()

png(filename = paste(figpath, "/imbalance_pacf.png", sep=''), width = wd, height = ht, res = res)
par(cex = cex)
par(mar = c(5, 4, 1.1, 1.1))
pacf(past$imb, main='')
dev.off()

#-----------------------------------------------------------------------------------
# Detection delay distributions

# Naive threshold
criterion = 'AGR'
sigma = 3; thr = -2*sigma
scenario = '0%'
p.frame = imb.PT.0.2sigma
filename = paste(figpath, "/cumDetProb_R", ifelse(scenario == '0%', 0, 1), "_", criterion, "_sigma", sigma, "_thr", thr, ".png", sep = "")
png(filename = filename, width = wd, height = ht, res = res)
par(cex = cex)
# prob.plot(p.frame, ylab='Detection probability', main=paste(criterion, 'sigma =', sigma, 'thr =', thr))
prob.plot(p.frame, ylab='Detection probability', legendloc = 'topleft')
i = which.min(abs(p.frame[[scenario]] - 0.95))
segments(x0 = future$t[i], y0=-0.1, y1=p.frame[[scenario]][i], lty=2)
dev.off()

scenario = '-1%'
p.frame = imb.PT.1.2sigma
filename = paste(figpath, "/cumDetProb_R", ifelse(scenario == '0%', 0, 1), "_", criterion, "_sigma", sigma, "_thr", thr, ".png", sep = "")
png(filename = filename, width = wd, height = ht, res = res)
par(cex = cex)
# prob.plot(p.frame, ylab='Detection probability', main=paste(criterion, 'sigma =', sigma, 'thr =', thr))
prob.plot(p.frame, ylab='Detection probability', legendloc = 'topleft')
i = which.min(abs(p.frame[[scenario]] - 0.95))
segments(x0 = future$t[i], y0=-0.1, y1=p.frame[[scenario]][i], lty=2)
dev.off()

# Calibrated threshold (AGR)
scenarios = c('0%', '-1%')
sigmas = c(3, 1.5, 0.75)
criterion = 'AGR'
for(k in 1:length(scenarios)) {
  for(j in 1:length(sigmas)) {
    scenario = scenarios[k]
    sigma = sigmas[j]
    ind = (k-1)*length(sigmas)+j  # row index
    thr = results.AGR$thr[ind]
    imb.PT = imb.PT.AGR[[k,j]]
  
    filename = paste(figpath, "/cumDetProb_", criterion, "_R", ifelse(scenario == '0%', 0, 1), "_sigma", sigma, "_thr", thr, ".png", sep = "")
    png(filename = filename, width = wd, height = ht, res = res)
    par(cex = cex)
    # prob.plot(imb.PT, ylab='Detection probability', main=paste(criterion, 'sigma =', sigma, 'thr =', thr))
    if(sigma == 3) legendloc = 'topleft' else legendloc = 'topright'
    prob.plot(imb.PT, ylab='Detection probability', legendloc = legendloc)
    i = results.AGR$delay.05[ind]+1; segments(x0 = future$t[i], y0=-0.1, y1=imb.PT[[scenario]][i], lty=2)
    i = results.AGR$delay.50[ind]+1; segments(x0 = future$t[i], y0=-0.1, y1=imb.PT[[scenario]][i], lty=2)
    i = results.AGR$delay.95[ind]+1; segments(x0 = future$t[i], y0=-0.1, y1=imb.PT[[scenario]][i], lty=2)
    dev.off()
  }
}

# Calibrated threshold (ACC)
scenarios = c('0%', '-1%')
sigmas = c(3, 1.5, 0.75)
criterion = 'ACC'
for(k in 1:length(scenarios)) {
  for(j in 1:length(sigmas)) {
    scenario = scenarios[k]
    sigma = sigmas[j]
    ind = (k-1)*length(sigmas)+j  # row index
    thr = results.ACC$thr[ind]
    Imb.PT = Imb.PT.ACC[[k,j]]
    
    filename = paste(figpath, "/cumDetProb_", criterion, "_R", ifelse(scenario == '0%', 0, 1), "_sigma", sigma, "_thr", thr, ".png", sep = "")
    png(filename = filename, width = wd, height = ht, res = res)
    par(cex = cex)
    # prob.plot(Imb.PT, ylab='Detection probability', main=paste(criterion, 'sigma =', sigma, 'thr =', thr))
    if(sigma == 3) legendloc = 'topleft' else legendloc = 'topright'
    prob.plot(Imb.PT, ylab='Detection probability', legendloc = legendloc)
    i = results.ACC$delay.05[ind]+1; segments(x0 = future$t[i], y0=-0.1, y1=Imb.PT[[scenario]][i], lty=2)
    i = results.ACC$delay.50[ind]+1; segments(x0 = future$t[i], y0=-0.1, y1=Imb.PT[[scenario]][i], lty=2)
    i = results.ACC$delay.95[ind]+1; segments(x0 = future$t[i], y0=-0.1, y1=Imb.PT[[scenario]][i], lty=2)
    dev.off()
  }
}


#-----------------------------------------------------------------------------------
# Positive predictive value

criterion = 'ACC'
pS = c(0.5, 0.2, 0.1) # prior probability of alternative scenario
# k = 1; j = 1  # index for scenario, sigma
k = 2; j = 1  # index for scenario, sigma

scenario = scenarios[k]
sigma = sigmas[j]
Imb.PT = Imb.PT.ACC[[k,j]]

filename = paste(figpath, "/PPV_", criterion, "_R", ifelse(scenario == '0%', 0, 1), "_sigma", sigma, ".png", sep = "")
png(filename = filename, width = wd, height = ht, res = res)
par(cex = cex); par(mar = c(2.2, 4.1, 1.1, 1.1))
# plot(NULL, xlim=future$t[c(1,N)], ylim=c(0,1), xlab='', ylab='Positive predictive value', main=paste(criterion, scenario, "sigma =", sigma))
plot(NULL, xlim=future$t[c(1,N)], ylim=c(0,1), xlab='', ylab='Positive predictive value')
for(i in 1:length(pS)) {
  pT.D = Imb.PT[[scenario]]*pS[i]/(Imb.PT[[scenario]]*pS[i] + Imb.PT$'1%'*(1-pS[i]))
  pT.D[pT.D == 1] = NaN
  lines(future$t, pT.D, col=ifelse(scenario == '0%', 'green', 'cyan'), lty=i, lwd=2)
}
dev.off()

#-----------------------------------------------------------------------------------
# Illustration of detection threshold

# ACC
criterion = 'ACC'
k = 1
scenario = results.ACC$scenario[k]
sigma = results.ACC$sigma[k] * 0.129
thr = results.ACC$thr[k] * 0.129
ACC.baseline = 403.3
mu1 = ACC.baseline + cumsum(future$`1%`) * 0.129
mu0 = ACC.baseline + cumsum(future$`0%`) * 0.129
delay.95 = future$t[1] + results.ACC$delay.95[k]

filename = paste(figpath, "/thr_example_", criterion, "_R", ifelse(scenario == '0%', 0, 1), "_sigma", sigma, ".png", sep = "")
png(filename = filename, width = wd, height = ht, res = res)
par(cex = cex); par(mar = c(2.2, 4.1, 1.1, 1.1))
plot(NULL, xlim=future$t[c(1,N)], ylim=c(400,500), xlab='', ylab='Atmospheric CO2 (ppm)')
lines(future$t, mu1, col='black', lty=1, lwd=2)
lines(future$t, mu1 + Imb[,38], col='orange', lty=1, lwd=2)
lines(future$t, mu0 + Imb[,3], col='green', lty=1, lwd=2)
lines(future$t, mu1 + thr, col='black', lty=2, lwd=1)
legend('topleft', legend=c('BAU','R0'), col=c('orange','green'), lwd=2)
dev.off()


#-----------------------------------------------------------------------------------
# Illustration of detection bands
two.lines = function(tt, mu, sigma=1, q=0, col.line='black', col.fill='gray', ltype=1, lwidth=2) {
  if(!is.na(col.fill))
    polygon(c(tt,rev(tt)), c(mu + q*sigma, rev(mu - q*sigma)), col=col.fill, border = NA)
  if(!is.na(col.line)) {
    lines(tt, mu + q*sigma, col=col.line, lty=ltype, lwd=lwidth)
    lines(tt, mu - q*sigma, col=col.line, lty=ltype, lwd=lwidth)
  }
}

# AGR
criterion = 'AGR'
k = 1
scenario = results.AGR$scenario[k]
sigma = results.AGR$sigma[k]
thr = results.AGR$thr[k]
mu1 = future$`1%`
mu0 = future$`0%`
delay.95 = future$t[1] + results.AGR$delay.95[k]

# q = qnorm(0.95)   # Detection time does not match separation for AGR
q = 2    # Two standard errors
my.gray = rgb(red=0.8, green=0.8, blue=0.8, alpha=0.5)

filename = paste(figpath, "/trajectories_", criterion, "_R", ifelse(scenario == '0%', 0, 1), "_sigma", sigma, "_highres.png", sep = "")
png(filename = filename, width = wd, height = ht, res = 300)
par(cex = cex); par(mar = c(2.2, 4.1, 1.1, 1.1))
plot(NULL, xlim=future$t[c(1,N)], ylim=c(0, 35), xlab='', ylab='Atmospheric growth rate (Gt CO2)')
two.lines(future$t, mu1, sigma, q=q, col.line=NA, col.fill=my.gray, ltype=3)
two.lines(future$t, mu0, sigma, q=q, col.line=NA, col.fill=my.gray, ltype=3)
lines(future$t, mu1, col='orange', lty=1, lwd=2)
lines(future$t, mu0, col='green', lty=1, lwd=2)
# lines(future$t, mu1 + thr, col='black', lty=2, lwd=1)
abline(v = delay.95, col='black', lty=2, lwd=1)
dev.off()

# ACC
criterion = 'ACC'
k = 1
scenario = results.ACC$scenario[k]
sigma = results.ACC$sigma[k] * 0.129
thr = results.ACC$thr[k] * 0.129
ACC.baseline = 403.3
mu1 = ACC.baseline + cumsum(future$`1%`) * 0.129
mu0 = ACC.baseline + cumsum(future$`0%`) * 0.129
delay.05 = future$t[1] + results.ACC$delay.05[k]
delay.50 = future$t[1] + results.ACC$delay.50[k]
delay.95 = future$t[1] + results.ACC$delay.95[k]

# q = qnorm(0.95)   # Detection time matches separation, only for ACC
q = 2    # Two standard errors
my.gray = rgb(red=0.8, green=0.8, blue=0.8, alpha=0.5)

filename = paste(figpath, "/trajectories_", criterion, "_R", ifelse(scenario == '0%', 0, 1), "_sigma", sigma, ".png", sep = "")
png(filename = filename, width = wd, height = ht, res = res)
par(cex = cex); par(mar = c(2.2, 4.1, 1.1, 1.1))
plot(NULL, xlim=future$t[c(1,N)], ylim=c(400,500), xlab='', ylab='Atmospheric CO2 (ppm)')
two.lines(future$t, mu1, sigma*sqrt(v), q=q, col.line=NA, col.fill=my.gray, ltype=3)
two.lines(future$t, mu0, sigma*sqrt(v), q=q, col.line=NA, col.fill=my.gray, ltype=3)
lines(future$t, mu1, col='orange', lty=1, lwd=2)
lines(future$t, mu0, col='green', lty=1, lwd=2)
# lines(future$t, mu1 + thr, col='black', lty=2, lwd=1)
abline(v = delay.05, col='black', lty=2, lwd=1)
abline(v = delay.50, col='black', lty=2, lwd=1)
abline(v = delay.95, col='black', lty=2, lwd=1)
legend('topleft', legend=c('BAU','R0'), col=c('orange','green'), lwd=2)
dev.off()

