###### Last try to fix the model-based interpolation:
library(tidyverse)
library(aghq)
library(TMB)
library(Matrix)
library(parallel)
library(INLA)
n_samp = 3000


source("/Users/ziangzhang/Documents/research_comps/RealData/function_used.R")
dyn.load(dynlib("Real_Smoothing"))

#### Load data, select rows
cUrl = paste0("http://scrippsco2.ucsd.edu/assets/data/atmospheric/",
              "stations/flask_co2/daily/daily_flask_co2_mlo.csv")
cFile = basename(cUrl)
if (!file.exists(cFile)) download.file(cUrl, cFile)
co2s = read.table(cFile, header = FALSE, sep = ",",
                  skip = 69, stringsAsFactors = FALSE, col.names = c("day",
                                                                     "time", "junk1", "junk2", "Nflasks", "quality",
                                                                     "co2"))

co2s$date = strptime(paste(co2s$day, co2s$time), format = "%Y-%m-%d %H:%M",
                     tz = "UTC")
# remove low-quality measurements
co2s = co2s[co2s$quality == 0, ]


co2s$day = as.Date(co2s$date)
toAdd = data.frame(day = seq(max(co2s$day) + 3, as.Date("2022/1/1"),
                             by = "7 days"), co2 = NA)
co2ext = rbind(co2s[, colnames(toAdd)], toAdd)
timeOrigin = as.Date("2000/1/1")

co2ext$timeYears = round(as.numeric(co2ext$day - timeOrigin)/365.25,
                         2)
co2ext$cos12 = cos(2 * pi * co2ext$timeYears)
co2ext$sin12 = sin(2 * pi * co2ext$timeYears)
co2ext$cos6 = cos(2 * 2 * pi * co2ext$timeYears)
co2ext$sin6 = sin(2 * 2 * pi * co2ext$timeYears)
co2ext$dayInt = as.integer(co2ext$day)

# ### Reduce the size of the dataset:
co2ext <- co2ext %>% filter(day >= "1990-01-01")
### Now use the full data, but look at a weekly observed grid.

allDays = seq(from = min(co2ext$day), to = max(co2ext$day),
              by = "7 day")

observed_dataset <- co2ext %>% filter(!is.na(co2ext$co2)) %>% select(c("co2", "cos12", "sin12", "cos6", "sin6", "dayInt", "day"))

designX <- as(as.matrix(observed_dataset[,-c(1,6,7)]), "dgTMatrix")



x_grid <- observed_dataset$dayInt
z <- as.integer(allDays)
z_grid <- z[!z %in% x_grid]
all_grid <- sort(c(x_grid, z_grid))
x_locations <- which(all_grid %in% x_grid)
myA <- construct_A(all_grid = all_grid, x_indx = x_locations)


n <- length(all_grid)
designB <- as(myA, "dgTMatrix")



d <- diff(all_grid)
H <- compute_H_rue(d,n = length(all_grid))
B <- compute_B(d,n = length(all_grid))
A <- compute_A(d, n = length(all_grid))


#####################################################################
#####################################################################
#####################################################################
#####################################################################
###### RW2:
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 1e-15), "dgTMatrix")



tmbdat <- list(
  # Design matrix
  X = designX,
  B = designB,
  # Penalty(Precision) matrix
  P = Q1,
  # Log determinant of penalty matrix (without the sigma part)
  # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
  logPdet = as.numeric(determinant(Q1,logarithm = T)$modulus),
  # Response
  y = observed_dataset$co2,
  # PC Prior params
  u1 = 0.001,
  alpha1 = 0.5,
  u2 = 1,
  alpha2 = 0.5,
  betaprec = 10^(-6)
)


tmbparams <- list(
  W = c(rep(0, n), rep(0, ncol(designX))), # W = c(U,beta); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "Real_Smoothing",
  silent = TRUE
)


# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)


# AGHQ
set.seed(123)
start_time <- Sys.time()
quad <- aghq::marginal_laplace_tmb(ff,5,c(0,0))
Sys.time() - start_time

gw <- sample_marginal(quad, n_samp)
mean_gw  <- apply(gw$samps,1, mean)
upper_gw  <- apply(gw$samps,1, quantile, p = 0.975)
lower_gw  <- apply(gw$samps,1, quantile, p = 0.025)

#### Just the random effects:
x_location <- which(all_grid %in% x_grid)
construct_X_matrix <- construct_A(all_grid, x_location)
gx <- construct_X_matrix %*% gw$samps[1:ncol(Q1), ]

mean_gx  <- apply(gx, 1, mean)
upper_gx  <- apply(gx, 1, quantile, p = 0.975)
lower_gx  <- apply(gx, 1, quantile, p = 0.025)

x_data <- data.frame(mean = mean_gx, x = x_grid, date = observed_dataset$day)
plot(mean~date, type = 'l', xlab = "Time", ylab = expression(f[np]), data = x_data)
lines(upper_gx~x_data$date, col = 'orange', lty = 'dashed')
lines(lower_gx~x_data$date, col = 'orange', lty = 'dashed')

### Parametric part:
para <- designX %*% gw$samps[-(1:ncol(Q1)), ]
f <- para + gx
mean_f  <- apply(f, 1, mean)
upper_f  <- apply(f, 1, quantile, p = 0.975)
lower_f  <- apply(f, 1, quantile, p = 0.025)

f_data <- data.frame(mean = mean_f, x = x_grid, date = observed_dataset$day)
plot(mean_f~date, type = 'l', xlab = "Time", ylab = expression(f[np] + f[p]), data = x_data)
lines(upper_f~f_data$date, col = 'orange', lty = 'dashed')
lines(lower_f~f_data$date, col = 'orange', lty = 'dashed')
points(observed_dataset$co2~observed_dataset$day, pch='.', col = 'blue')


# for (i in sample.int(n_samp,10)) {
#   lines(gx[,i] ~ x_data$date, col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
# }



### Check hyperparameter:
# Plot of theta1 posterior

theta_logprior <- function(theta, prior_alpha = tmbdat$alpha1, prior_u = tmbdat$u1) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2*log(x)))

prec_marg <- quad$marginals[[1]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam, type='l', xlab = expression(sigma[s]), ylab = "Posterior density" ))
lines(priorfuncsigma(seq(0.0001,0.001, by = 0.0001))~seq(0.0001,0.001, by = 0.0001), lty = 'dashed')




# Plot of theta2 posterior

theta_logprior <- function(theta, prior_alpha = tmbdat$alpha2, prior_u = tmbdat$u2) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2*log(x)))

prec_marg <- quad$marginals[[2]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l', xlab = expression(sigma[epsilon]), ylab = "Posterior density" ))
lines(priorfuncsigma(seq(0,1, by = 0.1))~seq(0,1, by = 0.1), lty = 'dashed')



### Takeout 
z_location <- which(all_grid %in% z)
construct_Z_matrix <- construct_A(all_grid, z_location)
gz <- construct_Z_matrix %*% gw$samps[1:ncol(Q1), ]

U_1st_Deriv <- apply(gz, 2, compute_deriv, order = 1) * 52
U_2nd_Deriv <- apply(gz, 2, compute_deriv, order = 2) * (52^2)


mean_gz  <- apply(gz,1, mean)
upper_gz  <- apply(gz,1, quantile, p = 0.975)
lower_gz  <- apply(gz,1, quantile, p = 0.025)

z_days <- allDays

z_data <- data.frame(mean = mean_gz, z = z, date = z_days)

mean_gz_1st  <- apply(U_1st_Deriv,1, mean)
upper_gz_1st  <- apply(U_1st_Deriv,1, quantile, p = 0.975)
lower_gz_1st  <- apply(U_1st_Deriv,1, quantile, p = 0.025)

z_1st <- data.frame(mean = mean_gz_1st, z = z[-1], date = z_days[-1])



plot(mean_gz ~ z_days, type = 'l', data = z_data)
for (i in sample.int(n_samp,10)) {
  lines(gz[,i] ~ z_days, col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
}


plot(mean_gz_1st ~ date, type = 'l', data = z_1st, ylim = c(-0.01,0.1))
lines(upper_gz_1st ~ date, lty = 'dashed', data = z_1st, col = 'orange')
lines(lower_gz_1st ~ date, lty = 'dashed', data = z_1st, col = 'orange')
for (i in sample.int(n_samp,10)) {
  lines(U_1st_Deriv[,i] ~ z_days[-1], col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
}


dev.copy(pdf,'2010_RW2_1st.pdf')
dev.off()


mean_gz_2nd  <- apply(U_2nd_Deriv,1, mean)
upper_gz_2nd  <- apply(U_2nd_Deriv,1, quantile, p = 0.975)
lower_gz_2nd  <- apply(U_2nd_Deriv,1, quantile, p = 0.025)

z_2nd <- data.frame(mean = mean_gz_2nd, z = z[-c(1,2)], date = z_days[-c(1,2)])

plot(mean_gz_2nd ~ date, type = 'l', data = z_2nd, ylim = c(-0.01,0.01))
lines(upper_gz_2nd~date, lty = 'dashed', data = z_2nd, col = 'orange')
lines(lower_gz_2nd~date, lty = 'dashed', data = z_2nd, col = 'orange')

for (i in sample.int(n_samp,3)) {
  lines(U_2nd_Deriv[,i] ~ z_days[-c(1,2)], col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
}

dev.copy(pdf,'2010_RW2_2nd.pdf')
dev.off()






#####################################################################
#####################################################################
#####################################################################
#####################################################################
###### ARIMA: Takes too long on the original data, reduce the dataset
cFile = basename(cUrl)
if (!file.exists(cFile)) download.file(cUrl, cFile)
co2s = read.table(cFile, header = FALSE, sep = ",",
                  skip = 69, stringsAsFactors = FALSE, col.names = c("day",
                                                                     "time", "junk1", "junk2", "Nflasks", "quality",
                                                                     "co2"))

co2s$date = strptime(paste(co2s$day, co2s$time), format = "%Y-%m-%d %H:%M",
                     tz = "UTC")
# remove low-quality measurements
co2s = co2s[co2s$quality == 0, ]


co2s$day = as.Date(co2s$date)
toAdd = data.frame(day = seq(max(co2s$day) + 3, as.Date("2022/1/1"),
                             by = "7 days"), co2 = NA)
co2ext = rbind(co2s[, colnames(toAdd)], toAdd)
timeOrigin = as.Date("2000/1/1")

co2ext$timeYears = round(as.numeric(co2ext$day - timeOrigin)/365.25,
                         2)
co2ext$cos12 = cos(2 * pi * co2ext$timeYears)
co2ext$sin12 = sin(2 * pi * co2ext$timeYears)
co2ext$cos6 = cos(2 * 2 * pi * co2ext$timeYears)
co2ext$sin6 = sin(2 * 2 * pi * co2ext$timeYears)
co2ext$dayInt = as.integer(co2ext$day)

# ### Reduce the size of the dataset:
co2ext <- co2ext %>% filter(day >= "2010-01-01")
### Now use the full data, but look at a weekly observed grid.
allDays = seq(from = min(co2ext$day), to = max(co2ext$day),
              by = "7 day")
observed_dataset <- co2ext %>% filter(!is.na(co2ext$co2)) %>% select(c("co2", "cos12", "sin12", "cos6", "sin6", "dayInt", "day"))
designX <- as(as.matrix(observed_dataset[,-c(1,6,7)]), "dgTMatrix")
x_grid <- observed_dataset$dayInt
z <- as.integer(allDays)
z_grid <- z[!z %in% x_grid]
all_grid <- sort(c(x_grid, z_grid))
x_locations <- which(all_grid %in% x_grid)
myA <- construct_A(all_grid = all_grid, x_indx = x_locations)
n <- length(all_grid)
designB <- as(myA, "dgTMatrix")

d <- diff(all_grid)
H <- compute_H_rue(d,n = length(all_grid))
B <- compute_B(d,n = length(all_grid))

D <- H[-c(1,n),]
R <- B[-c(1,n), -c(1,n)]

Q2 <- as(t(D) %*% solve(R) %*% D, 'dgTMatrix')
Q2 <- as(Q2 + Diagonal(n, x = 1e-15), "dgTMatrix")



tmbdat <- list(
  # Design matrix
  X = designX,
  B = designB,
  # Penalty(Precision) matrix
  P = Q2,
  # Log determinant of penalty matrix (without the sigma part)
  # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
  logPdet = as.numeric(determinant(Q1,logarithm = T)$modulus),
  # Response
  y = observed_dataset$co2,
  # PC Prior params
  u1 = 0.001,
  alpha1 = 0.5,
  u2 = 1,
  alpha2 = 0.5,
  betaprec = 10^(-6)
)


tmbparams <- list(
  W = c(rep(0, n), rep(0, ncol(designX))), # W = c(U,beta); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "Real_Smoothing",
  silent = TRUE
)


# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)


# AGHQ
set.seed(123)
start_time <- Sys.time()
quad2 <- aghq::marginal_laplace_tmb(ff,5,c(0,0))
Sys.time() - start_time
### Takes 11.87 minutes to run the code
set.seed(123)
gw <- sample_marginal(quad2, n_samp)

mean_gw  <- apply(gw$samps,1, mean)
upper_gw  <- apply(gw$samps,1, quantile, p = 0.975)
lower_gw  <- apply(gw$samps,1, quantile, p = 0.025)





### Check hyperparameter:
# Plot of theta1 posterior
prec_marg <- quad2$marginals[[1]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))




# Plot of theta2 posterior
prec_marg <- quad2$marginals[[2]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))




### Takeout 
z_location <- which(all_grid %in% z)
construct_Z_matrix <- construct_A(all_grid, z_location)
gz <- construct_Z_matrix %*% gw$samps[1:ncol(Q2), ]

U_1st_Deriv <- apply(gz, 2, compute_deriv, order = 1) * 52
U_2nd_Deriv <- apply(gz, 2, compute_deriv, order = 2) * (52^2)


mean_gz  <- apply(gz,1, mean)
upper_gz  <- apply(gz,1, quantile, p = 0.975)
lower_gz  <- apply(gz,1, quantile, p = 0.025)

z_days <- allDays

z_data <- data.frame(mean = mean_gz, z = z, date = z_days)

mean_gz_1st  <- apply(U_1st_Deriv,1, mean)
upper_gz_1st  <- apply(U_1st_Deriv,1, quantile, p = 0.975)
lower_gz_1st  <- apply(U_1st_Deriv,1, quantile, p = 0.025)

z_1st <- data.frame(mean = mean_gz_1st, z = z[-1], date = z_days[-1])



plot(mean_gz ~ z_days, type = 'l', data = z_data, xlab = "Time", ylab = expression(f[np]))
lines(upper_gz ~ z_days, lty = 'dashed', col = 'orange')
lines(lower_gz ~ z_days, lty = 'dashed', col = 'orange')
for (i in sample.int(n_samp,10)) {
  lines(gz[,i] ~ z_days, col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
}
mean(upper_gz - lower_gz)
### MCI: 0.6455676 for 99 percent
### MCI: 0.486293 for 95 percent
### MCI: 0.4071836 for 90 percent


plot(mean_gz_1st ~ date, type = 'l', data = z_1st, ylim = c(0,5), xlab = "Time", ylab = expression(paste(f[np], ":1st derivative ", sep = "\ ")))
lines(upper_gz_1st ~ date, lty = 'dashed', data = z_1st, col = 'orange')
lines(lower_gz_1st ~ date, lty = 'dashed', data = z_1st, col = 'orange')
for (i in sample.int(n_samp,10)) {
  lines(U_1st_Deriv[,i] ~ z_days[-1], col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
}

dev.copy(pdf,'2010_ARIMA_1st.pdf')
dev.off()


mean(upper_gz_1st - lower_gz_1st)
### MCI: 0.04734562 for 99 percent
### MCI: 0.03461681 for 95 percent
### MCI: 0.02857216 for 90 percent




mean_gz_2nd  <- apply(U_2nd_Deriv,1, mean)
upper_gz_2nd  <- apply(U_2nd_Deriv,1, quantile, p = 0.975)
lower_gz_2nd  <- apply(U_2nd_Deriv,1, quantile, p = 0.025)

z_2nd <- data.frame(mean = mean_gz_2nd, z = z[-c(1,2)], date = z_days[-c(1,2)])

plot(mean_gz_2nd ~ date, type = 'l', data = z_2nd, ylim = c(-25,25), xlab = "Time", ylab = expression(paste(f[np], ":2nd derivative ", sep = "\ ")))
lines(upper_gz_2nd~date, lty = 'dashed', data = z_2nd, col = 'orange')
lines(lower_gz_2nd~date, lty = 'dashed', data = z_2nd, col = 'orange')

for (i in sample.int(n_samp,3)) {
  lines(U_2nd_Deriv[,i] ~ z_days[-c(1,2)], col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
}

dev.copy(pdf,'2010_ARIMA_2nd.pdf')
dev.off()

mean(upper_gz_2nd - lower_gz_2nd)



### MCI: 0.01654983 for 99 percent
### MCI: 0.01174789 for 95 percent
### MCI: 0.009553871 for 90 percent







###### RW2 again, just to compare:
d <- diff(all_grid)
H <- compute_H_rue(d,n = length(all_grid))
B <- compute_B(d,n = length(all_grid))
A <- compute_A(d, n = length(all_grid))

Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 1e-15), "dgTMatrix")



tmbdat <- list(
  # Design matrix
  X = designX,
  B = designB,
  # Penalty(Precision) matrix
  P = Q1,
  # Log determinant of penalty matrix (without the sigma part)
  # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
  logPdet = as.numeric(determinant(Q1,logarithm = T)$modulus),
  # Response
  y = observed_dataset$co2,
  # PC Prior params
  u1 = 0.001,
  alpha1 = 0.5,
  u2 = 1,
  alpha2 = 0.5,
  betaprec = 10^(-6)
)


tmbparams <- list(
  W = c(rep(0, n), rep(0, ncol(designX))), # W = c(U,beta); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "Real_Smoothing",
  silent = TRUE
)


# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)


# AGHQ
set.seed(123)
start_time <- Sys.time()
quad <- aghq::marginal_laplace_tmb(ff,5,c(0,0))
Sys.time() - start_time
## takes 2.55 secs




set.seed(123)
gw <- sample_marginal(quad, n_samp)
mean_gw  <- apply(gw$samps,1, mean)
upper_gw  <- apply(gw$samps,1, quantile, p = 0.975)
lower_gw  <- apply(gw$samps,1, quantile, p = 0.025)





### Check hyperparameter:
# Plot of theta1 posterior
prec_marg <- quad$marginals[[1]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))




# Plot of theta2 posterior
prec_marg <- quad$marginals[[2]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))




### Takeout 
z_location <- which(all_grid %in% z)
construct_Z_matrix <- construct_A(all_grid, z_location)
gz <- construct_Z_matrix %*% gw$samps[1:ncol(Q1), ]

U_1st_Deriv <- apply(gz, 2, compute_deriv, order = 1) * 52
U_2nd_Deriv <- apply(gz, 2, compute_deriv, order = 2) * (52^2)


mean_gz  <- apply(gz,1, mean)
upper_gz  <- apply(gz,1, quantile, p = 0.975)
lower_gz  <- apply(gz,1, quantile, p = 0.025)

z_days <- allDays

z_data <- data.frame(mean = mean_gz, z = z, date = z_days)

mean_gz_1st  <- apply(U_1st_Deriv,1, mean)
upper_gz_1st  <- apply(U_1st_Deriv,1, quantile, p = 0.975)
lower_gz_1st  <- apply(U_1st_Deriv,1, quantile, p = 0.025)

z_1st <- data.frame(mean = mean_gz_1st, z = z[-1], date = z_days[-1])



plot(mean_gz ~ z_days, type = 'l', data = z_data, xlab = "Time", ylab = expression(f[np]))
lines(upper_gz ~ z_days, lty = 'dashed', col = 'orange')
lines(lower_gz ~ z_days, lty = 'dashed', col = 'orange')
for (i in sample.int(n_samp,10)) {
  lines(gz[,i] ~ z_days, col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
}
mean(upper_gz - lower_gz)


### MCI: 0.6481862 for 99 percent
### MCI: 0.4896828 for 95 percent
### MCI: 0.4084033 for 90 percent


plot(mean_gz_1st ~ date, type = 'l', data = z_1st, ylim = c(0,5), xlab = "Time", ylab = expression(paste(f[np], ":1st derivative ", sep = "\ ")))
lines(upper_gz_1st ~ date, lty = 'dashed', data = z_1st, col = 'orange')
lines(lower_gz_1st ~ date, lty = 'dashed', data = z_1st, col = 'orange')
for (i in sample.int(n_samp,10)) {
  lines(U_1st_Deriv[,i] ~ z_days[-1], col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
}

dev.copy(pdf,'2010_RW2_1st.pdf')
dev.off()

mean(upper_gz_1st - lower_gz_1st)


### MCI: 0.04770736 for 99 percent
### MCI: 0.03491584 for 95 percent
### MCI: 0.02883124 for 90 percent


mean_gz_2nd  <- apply(U_2nd_Deriv,1, mean)
upper_gz_2nd  <- apply(U_2nd_Deriv,1, quantile, p = 0.975)
lower_gz_2nd  <- apply(U_2nd_Deriv,1, quantile, p = 0.025)

z_2nd <- data.frame(mean = mean_gz_2nd, z = z[-c(1,2)], date = z_days[-c(1,2)])

plot(mean_gz_2nd ~ date, type = 'l', data = z_2nd, ylim = c(-25,25), xlab = "Time", ylab = expression(paste(f[np], ":2nd derivative ", sep = "\ ")))
lines(upper_gz_2nd~date, lty = 'dashed', data = z_2nd, col = 'orange')
lines(lower_gz_2nd~date, lty = 'dashed', data = z_2nd, col = 'orange')

for (i in sample.int(n_samp,3)) {
  lines(U_2nd_Deriv[,i] ~ z_days[-c(1,2)], col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
}

dev.copy(pdf,'2010_RW2_2nd.pdf')
dev.off()

mean(upper_gz_2nd - lower_gz_2nd)
### MCI: 0.01997148 for 99 percent
### MCI: 0.01419562 for 95 percent
### MCI: 0.01156123 for 90 percent















##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
#### Compare with INLA:


mm = get("inla.models", INLA:::inla.get.inlaEnv())
if(class(mm) == 'function') mm = mm()
mm$latent$crw2$min.diff = NULL
mm$latent$rw2$min.diff = NULL

co2res_full = inla(co2 ~ -1 + sin12 + cos12 + sin6 + cos6 +
                     f(dayInt, model = 'rw2', values = sort(c(z_grid,x_grid)), diagonal = 1e-10,
                       prior='pc.prec', param = c(0.001, 0.5), constr = FALSE),
                   data = observed_dataset, family='gaussian',
                   control.family = list(hyper=list(prec=list(
                     prior='pc.prec', param=c(1, 0.5)))),
                   # add this line if your computer has trouble
                   control.inla = list(strategy='gaussian'),
                   control.fixed = list(prec = 10^(-6)))

z_Grid_loc <- which(z %in% co2res_full$summary.random$dayInt$ID)
inla_pos_mean <- co2res_full$summary.random$dayInt$mean[z_Grid_loc]
inla_pos_ID <- co2res_full$summary.random$dayInt$ID[z_Grid_loc]
z_grid_inla <- data.frame(ID = inla_pos_ID, Mean = inla_pos_mean) %>% arrange(ID)
z_grid_inla$time <- allDays
z_grid_inla <- z_grid_inla[z_location,]
plot(Mean~time, data = z_grid_inla, type = 'l')

qCols = c('0.5quant','0.025quant','0.975quant')
1/sqrt(co2res$summary.hyperpar[,qCols])












### Study the maximum of curves of each method:
set.seed(123)
gw_rw2 <- sample_marginal(quad, n_samp)
gw_arima <- sample_marginal(quad2, n_samp)

gz_rw2 <- construct_Z_matrix %*% gw_rw2$samps[1:ncol(Q1), ]
gz_arima <- construct_Z_matrix %*% gw_arima$samps[1:ncol(Q2), ]

### maximum of the random effects:
gz_rw2_f_max <- apply(gz_rw2, 2, max)
summary(gz_rw2_f_max)
gz_arima_f_max <- apply(gz_arima, 2, max)
summary(gz_arima_f_max)

### maximum of the random effects 1st deriv:
U_1st_Deriv_rw2 <- apply(gz_rw2, 2, compute_deriv, order = 1)
U_2nd_Deriv_rw2 <- apply(gz_rw2, 2, compute_deriv, order = 2)

U_1st_Deriv_arima <- apply(gz_arima, 2, compute_deriv, order = 1)
U_2nd_Deriv_arima <- apply(gz_arima, 2, compute_deriv, order = 2)

gz_rw2_f1st_max <- apply(U_1st_Deriv_rw2, 2, max)
summary(gz_rw2_f1st_max)
hist(gz_rw2_f1st_max)
gz_arima_f1st_max <- apply(U_1st_Deriv_arima, 2, max)
summary(gz_arima_f1st_max)
hist(gz_arima_f1st_max)

gz_rw2_f2nd_max <- apply(U_2nd_Deriv_rw2, 2, max)
summary(gz_rw2_f2nd_max)
hist(gz_rw2_f2nd_max)
gz_arima_f2nd_max <- apply(U_2nd_Deriv_arima, 2, max)
summary(gz_arima_f2nd_max)
hist(gz_arima_f2nd_max)


f2nd_max_data <- tibble(Max = c(gz_arima_f2nd_max,gz_rw2_f2nd_max), Types = c(rep("ARIMA", length(gz_arima_f2nd_max)), rep("RW2", length(gz_rw2_f2nd_max))))

f2nd_max_data %>% ggplot(aes(x = Max, y = Types)) + geom_boxplot()



f1st_max_data <- tibble(Max = c(gz_arima_f1st_max,gz_rw2_f1st_max), Types = c(rep("ARIMA", length(gz_arima_f1st_max)), rep("RW2", length(gz_rw2_f1st_max))))

f1st_max_data %>% ggplot(aes(x = Max, y = Types)) + geom_boxplot()

























