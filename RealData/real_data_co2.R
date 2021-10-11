### Real Data Analysis:
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
plot(co2s$date, co2s$co2, log = "y", cex = 0.3, col = "#00000040",
     xlab = "time", ylab = "ppm")
plot(co2s[co2s$date > ISOdate(2015, 3, 1, tz = "UTC"),
          c("date", "co2")], log = "y", type = "o", xlab = "time",
     ylab = "ppm", cex = 0.5)


co2s$day = as.Date(co2s$date)
toAdd = data.frame(day = seq(max(co2s$day) + 3, as.Date("2025/1/1"),
                             by = "10 days"), co2 = NA)
co2ext = rbind(co2s[, colnames(toAdd)], toAdd)
timeOrigin = as.Date("2000/1/1")

co2ext$timeYears = round(as.numeric(co2ext$day - timeOrigin)/365.25,
                         2)
co2ext$cos12 = cos(2 * pi * co2ext$timeYears)
co2ext$sin12 = sin(2 * pi * co2ext$timeYears)
co2ext$cos6 = cos(2 * 2 * pi * co2ext$timeYears)
co2ext$sin6 = sin(2 * 2 * pi * co2ext$timeYears)
allDays = seq(from = min(co2ext$day), to = max(co2ext$day),
              by = "1 day")
co2ext$dayInt = as.integer(co2ext$day)




### Get the x_grid and z_grid:
library(tidyverse)
library(aghq)
library(TMB)
library(Matrix)
library(tidyverse)




observed_dataset <- co2ext %>% filter(!is.na(co2ext$co2)) %>% select(c("co2", "cos12", "sin12", "cos6", "sin6", "dayInt"))
### Start with debugging a simpler example:
# observed_dataset <- observed_dataset[1:1500, ]
x_grid <- observed_dataset$dayInt
z_grid <- co2ext$dayInt[is.na(co2ext$co2)]
n <- length(x_grid)

designX <- as(cbind(rep(1,n),as.matrix(observed_dataset[,-c(1,6)])), "dgTMatrix")
designB <- as(as.matrix(Diagonal(n)), "dgTMatrix")


compile("Real_Smoothing.cpp")
dyn.load(dynlib("Real_Smoothing"))











############ Implement:
n_samp = 1000
compute_H_rue <- function(d,n){
  H <- matrix(data = 0, nrow = n, ncol = n)
  for (i in 2:(nrow(H)-1)) {
    H[i,i] <- -(1/d[i-1]+1/d[i])
    H[i,i-1] <- 1/d[i-1]
    H[i,i+1] <- 1/d[i]
  }
  H
}
compute_B <- function(d,n){
  B <-matrix(0, nrow = n, ncol = n)
  B[1,1] <- d[1]/3
  B[1,2] <- d[1]/6
  B[n,n-1] <- d[n-1]/6
  B[n,n] <- d[n-1]/3
  for (i in 2:(nrow(B)-1)) {
    B[i,i-1] <- d[i-1]/6
    B[i,i] <- (d[i-1]+d[i])/3
    B[i,i+1] <- d[i]/6
  }
  B
}
compute_A <- function(d,n){
  A <-matrix(0, nrow = n, ncol = n)
  A[1,1] <- d[1]/2
  A[n,n] <- d[n-1]/2
  for (i in 2:(nrow(A)-1)) {
    A[i,i] <- (d[i-1]+d[i])/2
  }
  A
}


switch_matrix_once <- function(M, from, to){
  L <- diag(1, nrow = nrow(M), ncol = nrow(M))
  original_row <- L[to,]
  L[to,] <- L[from,]
  L[from,] <- original_row
  L%*%M%*%t(L)
}


switch_matrix <- function(z, x, M){
  all_grid <- sort(c(z,x))
  n <- length(x)
  for (i in 1:length(z)) {
    z_location <- which(all_grid %in% z[i])
    M <- switch_matrix_once(M,from = z_location, to = n + i)
    all_grid[z_location] <- all_grid[n+i]
    all_grid[n + i] <- z[i]
  }
  for (j in 1:length(x)) {
    x_location <- which(all_grid %in% x[j])
    M <- switch_matrix_once(M,from = x_location, to = j)
    all_grid[x_location] <- all_grid[j]
    all_grid[j] <- x[j]
  }
  M
}


Interpolation_vec_v1 <- function(t, x_grid, gx, GP, n_samples = 500){
  all_grid <- sort(c(t,x_grid))
  t_location <- which(all_grid %in% t)
  n <- length(all_grid)
  all_d <- diff(all_grid)
  H <- compute_H_rue(all_d, n)
  B <- compute_B(all_d, n)
  A <- compute_A(all_d, n)
  if(GP == "RW2"){
    Q <- t(H) %*% solve(A) %*% H
  }
  else{
    D <- H[-c(1,length(all_grid)),]
    R <- B[-c(1,length(all_grid)), -c(1,length(all_grid))]
    Q <- t(D) %*% solve(R) %*% D
  }
  
  C <- solve(Q + Diagonal(n, x = 0.00001))
  C_trans <- switch_matrix(t,x_grid,C)
  Q_trans <- Matrix::forceSymmetric(solve(C_trans))
  QAA <- Q_trans[(length(x_grid) + 1):n, (length(x_grid) + 1):n]
  QAB <- Q_trans[(length(x_grid) + 1):n, 1:(length(x_grid))]
  QBB <- Q[1:(length(x_grid)), 1:(length(x_grid))]
  compute_conditional_mean <- function(gx){
    conditional_mean <- -solve(QAA)%*% QAB %*% (gx)
    conditional_mean
  }
  conditional_mean_list <- apply(gx, 2, compute_conditional_mean)
  
  simulated_gz <- function(cmu){
    Lt <- chol(QAA)
    z <- rnorm(n = length(t))
    v <- solve(a = Lt, b = z)
    v + cmu
  }
  
  samples <- lapply(conditional_mean_list, simulated_gz)
  samples
}






### First Case: RW2 method
d <- diff(observed_dataset$dayInt)
H <- compute_H_rue(d,n = length(observed_dataset$dayInt))
B <- compute_B(d,n = length(observed_dataset$dayInt))
H <- compute_H_rue(d,n = length(observed_dataset$dayInt))
B <- compute_B(d,n = length(observed_dataset$dayInt))
A <- compute_A(d, n = length(observed_dataset$dayInt))


###### RW2:
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 0.00001), "dgTMatrix")


tmbdat <- list(
  # Design matrix
  X = designX,
  B = designB,
  # Penalty(Precision) matrix
  P = Q1,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q1,logarithm = TRUE)$modulus),
  # Response
  y = observed_dataset$co2,
  # PC Prior params
  u1 = 0.001,
  alpha1 = 0.5,
  u2 = 1,
  alpha2 = 0.5,
  betaprec = 0.001
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
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))
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



### Beta part:
beta_mean <- mean_gw[(ncol(Q1)+1):(ncol(Q1)+5)]


### Smoothing part:
U <- gw$samps[1:n,]
U_mean <- mean_gw[1:ncol(Q1)]
U_upper <- upper_gw[1:ncol(Q1)]
U_lower <- lower_gw[1:ncol(Q1)]

plot(x = co2ext$day[1:ncol(Q1)], y = U_mean, type = 'l', lty = 'solid', ylim = c(-60,60), ylab = "Random Effects", xlab = "time")
lines(U_upper~co2ext$day[1:ncol(Q1)], lty = 'dashed')
lines(U_lower~co2ext$day[1:ncol(Q1)], lty = 'dashed')


### Overall Curve
f_par <-  designX %*% c(gw$samps[(ncol(Q1)+1):(ncol(Q1)+5),1])
for (i in 2:ncol(gw$samps)) {
  f_par <- cbind(f_par, designX %*% c(gw$samps[(ncol(Q1)+1):(ncol(Q1)+5),i]))
}
mean_f_par <- apply(f_par,1, mean)

plot(x = co2ext$day[1:ncol(Q1)], y = mean_f_par, type = 'l', lty = 'solid', ylab = "Parametric Effects", xlab = "time")

f_overall <- f_par + U
mean_f_overall <- apply(f_overall,1, mean)
upper_f_overall <- apply(f_overall,1, quantile, p = 0.975)
lower_f_overall <- apply(f_overall,1, quantile, p = 0.025)

plot(x = co2ext$day[1:ncol(Q1)], y = mean_f_overall, type = 'l', lty = 'solid', ylab = "Overall Effects", xlab = "time")
lines(upper_f_overall~co2ext$day[1:ncol(Q1)], lty = 'dashed')
lines(lower_f_overall~co2ext$day[1:ncol(Q1)], lty = 'dashed')

### RE Derivative and Functional Derivative:
compute_deriv <- function(vec, order){diff(vec, differences = order)}

## For random effects:
U_1st_Deriv <- apply(U, 2, compute_deriv, order = 1)
U_2nd_Deriv <- apply(U, 2, compute_deriv, order = 2)

U_1st_mean <- apply(U_1st_Deriv, 1, mean)
U_1st_upper <- apply(U_1st_Deriv, 1, quantile, p = 0.975)
U_1st_lower <- apply(U_1st_Deriv, 1, quantile, p = 0.025)

U_2nd_mean <- apply(U_2nd_Deriv, 1, mean)
U_2nd_upper <- apply(U_2nd_Deriv, 1, quantile, p = 0.975)
U_2nd_lower <- apply(U_2nd_Deriv, 1, quantile, p = 0.025)

plot(x = co2ext$day[2:ncol(Q1)], y = U_1st_mean, type = 'l', lty = 'solid', ylim = c(-4,4), ylab = "Random Effects 1st Deriv", xlab = "time")
lines(U_1st_upper~co2ext$day[2:ncol(Q1)], lty = 'dashed')
lines(U_1st_lower~co2ext$day[2:ncol(Q1)], lty = 'dashed')
# for (i in sample.int(n_samp,30)) {
#   lines(unlist(U_1st_Deriv[,i]) ~ co2ext$day[2:ncol(Q1)], col = rgb(240, 0, 0, max = 255, alpha = 40, names = "grey"))
# }

plot(x = co2ext$day[3:ncol(Q1)][-c(1:2000)], y = U_2nd_mean[-c(1:2000)], type = 'l', lty = 'solid', ylim = c(-4,4), ylab = "Random Effects 2nd Deriv", xlab = "time", xlim = )
lines(U_2nd_upper~co2ext$day[3:ncol(Q1)], lty = 'dashed')
lines(U_2nd_lower~co2ext$day[3:ncol(Q1)], lty = 'dashed')

## Above derivatives are not ideal, should be based on the Z_grid instead.

##################### Inference on high-resolution grid Z:
z <- as.integer(allDays)
z_grid <- z[!z %in% x_grid]
gz_list <- Interpolation_vec_v1(t = z_grid, x_grid, U, "RW2")
gz <- Matrix(0,nrow = length(z_grid), ncol = n_samp)
for (i in 1:length(gz_list)) {
  gz[,i] <- gz_list[[i]]
}











