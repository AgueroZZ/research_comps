library(INLA)
library(mgcv)
library(TMB)
library(aghq)



#################Generating Gaussian Data for Smoothing ##############
######################################################################
#### Assuming true sigma is known to be 1
n <- 500
x <- seq(0.1,50, by = 0.1)
a <- min(x)
b <- max(x)
y <- 5*sin(0.5*x) + rnorm(n, sd = 1)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")



#### Approach 1: Using RW2 model with Diaognal Approximation
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
H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 0.0001), "dgTMatrix")



compile("01_RW2Comparison.cpp")
dyn.load(dynlib("01_RW2Comparison"))

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q1,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q1,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u = 2,
  alpha = 0.5
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta = 0 # -2log(sigma)
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "01_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,0)

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))

# Use evaluation vector to reconstruct the fitted function:
construct_fitted <- function(weights){
  knots <- x
  splineFunc <- splines::bs(location_interest, degree = 1, knots = knots)
  fitted_func <- splineFunc[,-(length(knots) + 1)] %*% matrix(weights, ncol = 1)
  as.numeric(fitted_func)
}

compute_deriv <- function(values){
  values <- diff(values)
  values
}



# Inference for W
samps1 <- sample_marginal(quad,1e03)
resolution = 0.001
location_interest <- seq(a,b, by = resolution)
fitted_func <- construct_fitted(weights = samps1$samps[,1])
plot(fitted_func~location_interest, col = "blue", type = "l")

# Construct samples for g(.):
sample_func <- apply(samps1$samps, MARGIN = 2, FUN = construct_fitted)
sample_deriv <- apply(sample_func, 2, compute_deriv)
sample_deriv_second <- apply(sample_deriv, 2, compute_deriv)



# Posterior mean
W1 <- apply(samps1$samps,1,mean)
mean_func <- apply(sample_func,1,mean)
upper_func <- apply(sample_func,1,quantile, probs = 0.975)
lower_func <- apply(sample_func,1,quantile, probs = 0.025)


### Plot of function:
plot(mean_func~location_interest, col = "red", type = "l", ylim=c(-6,6), xlab = "x", ylab = "g(.)")
lines(upper_func ~ location_interest, col = "orange", lty = "dashed")
lines(lower_func ~ location_interest, col = "orange", lty = "dashed")

for (i in sample.int(1000,50)) {
  lines(sample_func[,i] ~ location_interest, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}


### Plot of function derivative:
mean_deriv <- apply(sample_deriv,1,mean)
upper_deriv <- apply(sample_deriv,1,quantile, probs = 0.975)
lower_deriv <- apply(sample_deriv,1,quantile, probs = 0.025)


plot(mean_deriv~location_interest[-1], col = "red", type = "l", ylab = "1st derivative", xlab = "x")
lines(upper_deriv ~ location_interest[-1], col = "orange", lty = "dashed")
lines(lower_deriv ~ location_interest[-1], col = "orange", lty = "dashed")

for (i in sample.int(1000,5)) {
  lines(sample_deriv[,i] ~ location_interest[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}


### Plot of function second derivative:
mean_deriv_second <- apply(sample_deriv_second,1,mean)
upper_deriv_second <- apply(sample_deriv_second,1,quantile, probs = 0.975)
lower_deriv_second <- apply(sample_deriv_second,1,quantile, probs = 0.025)


plot(mean_deriv_second~location_interest[-(1:2)], col = "blue", type = "l", ylim=c(-0.0008,0.0008), ylab = "2nd derivative", xlab = "x")
# lines(upper_deriv_second ~ location_interest[-(1:2)], col = "orange", lty = "dashed")
# lines(lower_deriv_second ~ location_interest[-(1:2)], col = "orange", lty = "dashed")

for (i in sample.int(1000,1)) {
  lines(sample_deriv_second[,i] ~ location_interest[-(1:2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

### Look at the second order differences at original locations (since spacings are equal)
second_diff <- apply(samps1$samps,2,diff, differences = 2)
mean_second_diff <- apply(second_diff,1,mean)
upper_second_diff <- apply(second_diff,1,quantile, probs = 0.975)
lower_second_diff <- apply(second_diff,1,quantile, probs = 0.025)

### Plot of second order differences:
plot(mean_second_diff ~ x[-c(1,2)], col = "red", type = "l", ylim=c(-0.1,0.1), xlab = "x", ylab = "second order difference")
lines(upper_second_diff ~ x[-c(1,2)], col = "orange", lty = "dashed")
lines(lower_second_diff ~ x[-c(1,2)], col = "orange", lty = "dashed")

for (i in sample.int(1000,10)) {
  lines(second_diff[,i] ~ x[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

acf(second_diff[,1], lag.max = 10)
acf(mean_second_diff, lag.max = 10)



#### Approach 2: Using RW2 model without Diagonal Approximation
Q2 <- t(H) %*% solve(B) %*% H
Q2 <- as(as.matrix(Q2 + Diagonal(n, x = 0.0001)), "dgTMatrix")

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q2,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q2,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u = 2,
  alpha = 0.5
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta = 0 # -2log(sigma)
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "01_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,0)

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))


# Inference for W
samps2 <- sample_marginal(quad,1e03)
resolution = 0.001
location_interest <- seq(a,b, by = resolution)
fitted_func <- construct_fitted(weights = samps2$samps[,1])
plot(fitted_func~location_interest, col = "blue", type = "l")

# Construct samples for g(.):
sample_func <- apply(samps2$samps, MARGIN = 2, FUN = construct_fitted)
sample_deriv <- apply(sample_func, 2, compute_deriv)
sample_deriv_second <- apply(sample_deriv, 2, compute_deriv)




# Posterior mean
W2 <- apply(samps2$samps,1,mean)
mean_func <- apply(sample_func,1,mean)
upper_func <- apply(sample_func,1,quantile, probs = 0.975)
lower_func <- apply(sample_func,1,quantile, probs = 0.025)


### Plot of function:
plot(mean_func~location_interest, col = "red", type = "l", ylim=c(-6,6), xlab = "x", ylab = "g(.)")
lines(upper_func ~ location_interest, col = "orange", lty = "dashed")
lines(lower_func ~ location_interest, col = "orange", lty = "dashed")

for (i in sample.int(1000,50)) {
  lines(sample_func[,i] ~ location_interest, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}


### Plot of function derivative:
mean_deriv <- apply(sample_deriv,1,mean)
upper_deriv <- apply(sample_deriv,1,quantile, probs = 0.975)
lower_deriv <- apply(sample_deriv,1,quantile, probs = 0.025)


plot(mean_deriv~location_interest[-1], col = "red", type = "l", ylab = "1st derivative", xlab = "x")
lines(upper_deriv ~ location_interest[-1], col = "orange", lty = "dashed")
lines(lower_deriv ~ location_interest[-1], col = "orange", lty = "dashed")

for (i in sample.int(1000,5)) {
  lines(sample_deriv[,i] ~ location_interest[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}


### Plot of function second derivative:
mean_deriv_second <- apply(sample_deriv_second,1,mean)
upper_deriv_second <- apply(sample_deriv_second,1,quantile, probs = 0.975)
lower_deriv_second <- apply(sample_deriv_second,1,quantile, probs = 0.025)


plot(mean_deriv_second~location_interest[-(1:2)], col = "blue", type = "l", ylim=c(-0.0008,0.0008), ylab = "2nd derivative", xlab = "x")
# lines(upper_deriv_second ~ location_interest[-(1:2)], col = "orange", lty = "dashed")
# lines(lower_deriv_second ~ location_interest[-(1:2)], col = "orange", lty = "dashed")

for (i in sample.int(1000,1)) {
  lines(sample_deriv_second[,i] ~ location_interest[-(1:2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

### Look at the second order differences at original locations (since spacings are equal)
second_diff <- apply(samps1$samps,2,diff, differences = 2)
mean_second_diff <- apply(second_diff,1,mean)
upper_second_diff <- apply(second_diff,1,quantile, probs = 0.975)
lower_second_diff <- apply(second_diff,1,quantile, probs = 0.025)

### Plot of second order differences:
plot(mean_second_diff ~ x[-c(1,2)], col = "red", type = "l", ylim=c(-0.1,0.1), xlab = "x", ylab = "second order difference")
lines(upper_second_diff ~ x[-c(1,2)], col = "orange", lty = "dashed")
lines(lower_second_diff ~ x[-c(1,2)], col = "orange", lty = "dashed")

for (i in sample.int(1000,10)) {
  lines(second_diff[,i] ~ x[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}



#### Approach 3: Using ARIMA method
D <- H[-c(1,n),]
R <- B[-c(1,n), -c(1,n)]

Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(n, x = 0.0001)), "dgTMatrix")

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q3,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q3,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u = 2,
  alpha = 0.5
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta = 0 # -2log(sigma)
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "01_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,0)

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))


# Inference for W
samps3 <- sample_marginal(quad,1e03)
W3 <- apply(samps3$samps,1,mean)


construct_fitted_cubic <- function(weights){
  spline_fitted <- spline(x, y = weights, xout = location_interest, method = "natural",
                          xmin = min(x), xmax = max(x), ties = mean)
  spline_fitted$y
}




######## Cubic way:
sample_func_cubic <- apply(samps3$samps, MARGIN = 2, FUN = construct_fitted_cubic)
mean_func_cubic <- apply(sample_func_cubic,1,mean)
sample_deriv_cubic <- apply(sample_func_cubic, 2, compute_deriv)
sample_deriv_second_cubic <- apply(sample_deriv_cubic, 2, compute_deriv)
upper_func_cubic <- apply(sample_func_cubic,1,quantile, probs = 0.975)
lower_func_cubic <- apply(sample_func_cubic,1,quantile, probs = 0.025)


### Plot of function:
plot(mean_func_cubic~location_interest, col = "red", type = "l", ylim=c(-6,6), xlab = "x", ylab = "g(.)")
lines(upper_func_cubic ~ location_interest, col = "orange", lty = "dashed")
lines(lower_func_cubic ~ location_interest, col = "orange", lty = "dashed")

for (i in sample.int(1000,50)) {
  lines(sample_func[,i] ~ location_interest, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

### Plot of function derivative:
mean_deriv_cubic <- apply(sample_deriv_cubic,1,mean)
upper_deriv_cubic <- apply(sample_deriv_cubic,1,quantile, probs = 0.975)
lower_deriv_cubic <- apply(sample_deriv_cubic,1,quantile, probs = 0.025)


plot(mean_deriv_cubic~location_interest[-1], col = "red", type = "l", ylab = "1st derivative", xlab = "x")
lines(upper_deriv_cubic ~ location_interest[-1], col = "orange", lty = "dashed")
lines(lower_deriv_cubic ~ location_interest[-1], col = "orange", lty = "dashed")

for (i in sample.int(1000,5)) {
  lines(sample_deriv_cubic[,i] ~ location_interest[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

### Plot of function 2nd derivative:
mean_deriv_second_cubic <- apply(sample_deriv_second_cubic,1,mean)
upper_deriv_second_cubic <- apply(sample_deriv_second_cubic,1,quantile, probs = 0.975)
lower_deriv_second_cubic <- apply(sample_deriv_second_cubic,1,quantile, probs = 0.025)


plot(mean_deriv_second_cubic~location_interest[-c(1,2)], col = "red", type = "l", ylab = "2nd derivative", xlab = "x", ylim = c(-9e-06,9e-06))
# lines(upper_deriv_second_cubic ~ location_interest[-c(1,2)], col = "orange", lty = "dashed")
# lines(lower_deriv_second_cubic ~ location_interest[-c(1,2)], col = "orange", lty = "dashed")

for (i in sample.int(1000,1)) {
  lines(sample_deriv_second_cubic[,i] ~ location_interest[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}


### Look at the second order differences at original locations (since spacings are equal)
second_diff_cubic <- apply(samps3$samps,2,diff, differences = 2)
mean_second_diff_cubic <- apply(second_diff_cubic,1,mean)
upper_second_diff_cubic <- apply(second_diff_cubic,1,quantile, probs = 0.975)
lower_second_diff_cubic <- apply(second_diff_cubic,1,quantile, probs = 0.025)

### Plot of second order differences:
plot(mean_second_diff_cubic ~ x[-c(1,2)], col = "red", type = "l", ylim=c(-0.1,0.1), xlab = "x", ylab = "second order difference")
lines(upper_second_diff_cubic ~ x[-c(1,2)], col = "orange", lty = "dashed")
lines(lower_second_diff_cubic ~ x[-c(1,2)], col = "orange", lty = "dashed")

for (i in sample.int(1000,10)) {
  lines(second_diff_cubic[,i] ~ x[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}




######## Linear way:
sample_func_linear <- apply(samps3$samps, MARGIN = 2, FUN = construct_fitted)
mean_func_linear <- apply(sample_func_linear,1,mean)
upper_func_linear <- apply(sample_func_linear,1,quantile, prob = 0.975)
lower_func_linear <- apply(sample_func_linear,1,quantile, prob = 0.025)


### Plot of function:
plot(mean_func_linear~location_interest, col = "red", type = "l", ylim=c(-6,6), xlab = "x", ylab = "g(.)")
lines(upper_func_linear ~ location_interest, col = "orange", lty = "dashed")
lines(lower_func_linear ~ location_interest, col = "orange", lty = "dashed")

for (i in sample.int(1000,50)) {
  lines(sample_func_linear[,i] ~ location_interest, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}


sample_deriv_linear <- apply(sample_func_linear, 2, compute_deriv)
sample_deriv_second_linear <- apply(sample_deriv_linear, 2, compute_deriv)
upper_deriv_linear <- apply(sample_deriv_linear,1,quantile, probs = 0.975)
lower_deriv_linear <- apply(sample_deriv_linear,1,quantile, probs = 0.025)
mean_deriv_linear <- apply(sample_deriv_linear, 1, mean)

### Plot of function derivative:

plot(mean_deriv_linear~location_interest[-1], col = "red", type = "l", ylab = "1st derivative", xlab = "x")
lines(upper_deriv_linear ~ location_interest[-1], col = "orange", lty = "dashed")
lines(lower_deriv_linear ~ location_interest[-1], col = "orange", lty = "dashed")

for (i in sample.int(1000,5)) {
  lines(sample_deriv_linear[,i] ~ location_interest[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

### Plot of function 2nd derivative:
mean_deriv_second_linear <- apply(sample_deriv_second_linear,1,mean)
upper_deriv_second_linear <- apply(sample_deriv_second_linear,1,quantile, probs = 0.975)
lower_deriv_second_linear <- apply(sample_deriv_second_linear,1,quantile, probs = 0.025)


plot(mean_deriv_second_linear~location_interest[-c(1,2)], col = "red", type = "l", ylab = "2nd derivative", xlab = "x", ylim = c(-0.001,0.001))
# lines(upper_deriv_second_cubic ~ location_interest[-c(1,2)], col = "orange", lty = "dashed")
# lines(lower_deriv_second_cubic ~ location_interest[-c(1,2)], col = "orange", lty = "dashed")

for (i in sample.int(1000,1)) {
  lines(sample_deriv_second_linear[,i] ~ location_interest[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}


### Look at the second order differences at original locations (since spacings are equal)
second_diff_linear <- apply(samps3$samps,2,diff, differences = 2)
mean_second_diff_linear <- apply(second_diff_linear,1,mean)
upper_second_diff_linear <- apply(second_diff_linear,1,quantile, probs = 0.975)
lower_second_diff_linear <- apply(second_diff_linear,1,quantile, probs = 0.025)

### Plot of second order differences:
plot(mean_second_diff_linear ~ x[-c(1,2)], col = "red", type = "l", ylim=c(-0.1,0.1), xlab = "x", ylab = "second order difference")
lines(upper_second_diff_linear ~ x[-c(1,2)], col = "orange", lty = "dashed")
lines(lower_second_diff_linear ~ x[-c(1,2)], col = "orange", lty = "dashed")

for (i in sample.int(1000,10)) {
  lines(second_diff_linear[,i] ~ x[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
































#####################################################
#####################################################
####### If true sigma is unknown with a PC prior ####


n <- 500
x <- seq(0.1,50, by = 0.1)
y <- 5*sin(0.5*x) + rnorm(n, sd = 1)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")



#### Approach 1: Using RW2 model with Diaognal Approximation
H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 0.0001), "dgTMatrix")



compile("02_RW2Comparison.cpp")
dyn.load(dynlib("02_RW2Comparison"))

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q1,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q1,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u1 = 2,
  alpha1 = 0.5,
  u2 = 2,
  alpha2 = 0.5
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))



# Plot of theta posterior
logpostsigma1 <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma1,plot(transparam,pdf_transparam,type='l'))
logpostsigma2 <- compute_pdf_and_cdf(quad$marginals[[2]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma2,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps1 <- sample_marginal(quad,1e03)
# Posterior mean
W1 <- apply(samps1$samps,1,mean)
L1 <- apply(samps1$samps,1,quantile, probs = 0.025)
U1 <- apply(samps1$samps,1,quantile, probs = 0.975)



#### Approach 2: Using RW2 model without Diaognal Approximation
Q2 <- t(H) %*% solve(B) %*% H
Q2 <- as(as.matrix(Q2 + Diagonal(n, x = 0.0001)), "dgTMatrix")
tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q2,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q2,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u1 = 2,
  alpha1 = 0.5,
  u2 = 2,
  alpha2 = 0.5
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "03_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))
# Plot of theta posterior
logpostsigma1 <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma1,plot(transparam,pdf_transparam,type='l'))
logpostsigma2 <- compute_pdf_and_cdf(quad$marginals[[2]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma2,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps2 <- sample_marginal(quad,1e03)
# Posterior mean
W2 <- apply(samps2$samps,1,mean)
L2 <- apply(samps2$samps,1,quantile, probs = 0.025)
U2 <- apply(samps2$samps,1,quantile, probs = 0.975)













#### Approach 3: Using ARIMA method
D <- H[-c(1,n),]
R <- B[-c(1,n), -c(1,n)]

Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(n, x = 0.0001)), "dgTMatrix")

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q3,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q3,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u1 = 2,
  alpha1 = 0.5,
  u2 = 2,
  alpha2 = 0.5
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))
# Plot of theta posterior
logpostsigma1 <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma1,plot(transparam,pdf_transparam,type='l'))
logpostsigma2 <- compute_pdf_and_cdf(quad$marginals[[2]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma2,plot(transparam,pdf_transparam,type='l'))
# Inference for W
samps3 <- sample_marginal(quad,1e03)
# Posterior mean
W3 <- apply(samps3$samps,1,mean)
L3 <- apply(samps3$samps,1,quantile, probs = 0.025)
U3 <- apply(samps3$samps,1,quantile, probs = 0.975)



#### Comparison
plot(y~x)
lines(W1~x, col = "Red")
lines(U1~x, col = "Red", lty=2)
lines(L1~x, col = "Red", lty=2)


lines(W2~x, col = "blue")
lines(U2~x, col = "blue", lty=2)
lines(L2~x, col = "blue", lty=2)

lines(W3~x, col = "green")
lines(U3~x, col = "green", lty=2)
lines(L3~x, col = "green", lty=2)
lines(5*sin(0.5*x)~x, col = "black")






#####################################################
####### try a different PC prior ####

n <- 500
x <- seq(0.1,50, by = 0.1)
y <- 5*sin(0.5*x) + rnorm(n, sd = 1)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")



#### Approach 1: Using RW2 model with Diaognal Approximation
H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 0.0001), "dgTMatrix")



compile("02_RW2Comparison.cpp")
dyn.load(dynlib("02_RW2Comparison"))

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q1,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q1,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u1 = 5,
  alpha1 = 0.1,
  u2 = 5,
  alpha2 = 0.1
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))



# Plot of theta posterior
logpostsigma1 <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma1,plot(transparam,pdf_transparam,type='l'))
logpostsigma2 <- compute_pdf_and_cdf(quad$marginals[[2]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma2,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps1 <- sample_marginal(quad,1e03)
# Posterior mean
W1 <- apply(samps1$samps,1,mean)
L1 <- apply(samps1$samps,1,quantile, probs = 0.025)
U1 <- apply(samps1$samps,1,quantile, probs = 0.975)



#### Approach 2: Using RW2 model without Diaognal Approximation
Q2 <- t(H) %*% solve(B) %*% H
Q2 <- as(as.matrix(Q2 + Diagonal(n, x = 0.0001)), "dgTMatrix")
tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q2,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q2,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u1 = 5,
  alpha1 = 0.1,
  u2 = 5,
  alpha2 = 0.1
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))
# Plot of theta posterior
logpostsigma1 <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma1,plot(transparam,pdf_transparam,type='l'))
logpostsigma2 <- compute_pdf_and_cdf(quad$marginals[[2]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma2,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps2 <- sample_marginal(quad,1e03)
# Posterior mean
W2 <- apply(samps2$samps,1,mean)
L2 <- apply(samps2$samps,1,quantile, probs = 0.025)
U2 <- apply(samps2$samps,1,quantile, probs = 0.975)













#### Approach 3: Using ARIMA method
D <- H[-c(1,n),]
R <- B[-c(1,n), -c(1,n)]

Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(n, x = 0.0001)), "dgTMatrix")

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q3,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q3,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u1 = 5,
  alpha1 = 0.1,
  u2 = 5,
  alpha2 = 0.1
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))
# Plot of theta posterior
logpostsigma1 <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma1,plot(transparam,pdf_transparam,type='l'))
logpostsigma2 <- compute_pdf_and_cdf(quad$marginals[[2]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma2,plot(transparam,pdf_transparam,type='l'))
# Inference for W
samps3 <- sample_marginal(quad,1e03)
# Posterior mean
W3 <- apply(samps3$samps,1,mean)
L3 <- apply(samps3$samps,1,quantile, probs = 0.025)
U3 <- apply(samps3$samps,1,quantile, probs = 0.975)



#### Comparison
plot(y~x)
lines(W1~x, col = "Red")
lines(U1~x, col = "Red", lty=2)
lines(L1~x, col = "Red", lty=2)


lines(W2~x, col = "blue")
lines(U2~x, col = "blue", lty=2)
lines(L2~x, col = "blue", lty=2)

lines(W3~x, col = "green")
lines(U3~x, col = "green", lty=2)
lines(L3~x, col = "green", lty=2)
lines(5*sin(0.5*x)~x, col = "black")


































########################################
########################################
############ Very Sparse knots #########
n <- 500
x <- seq(2,50, by = 2)
x <- rep(x, each = 20)
y <- 5*sin(0.5*x) + rnorm(n, sd = 3)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")




#### Approach 1: Using RW2 model with Diaognal Approximation
H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 0.0001), "dgTMatrix")



compile("02_RW2Comparison.cpp")
dyn.load(dynlib("02_RW2Comparison"))

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q1,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q1,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u1 = 5,
  alpha1 = 0.1,
  u2 = 5,
  alpha2 = 0.1
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))



# Plot of theta posterior
logpostsigma1 <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma1,plot(transparam,pdf_transparam,type='l'))
logpostsigma2 <- compute_pdf_and_cdf(quad$marginals[[2]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma2,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps1 <- sample_marginal(quad,1e03)
# Posterior mean
W1 <- apply(samps1$samps,1,mean)
L1 <- apply(samps1$samps,1,quantile, probs = 0.025)
U1 <- apply(samps1$samps,1,quantile, probs = 0.975)



#### Approach 2: Using RW2 model without Diaognal Approximation
Q2 <- t(H) %*% solve(B) %*% H
Q2 <- as(as.matrix(Q2 + Diagonal(n, x = 0.0001)), "dgTMatrix")
tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q2,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q2,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u1 = 5,
  alpha1 = 0.1,
  u2 = 5,
  alpha2 = 0.1
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))
# Plot of theta posterior
logpostsigma1 <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma1,plot(transparam,pdf_transparam,type='l'))
logpostsigma2 <- compute_pdf_and_cdf(quad$marginals[[2]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma2,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps2 <- sample_marginal(quad,1e03)
# Posterior mean
W2 <- apply(samps2$samps,1,mean)
L2 <- apply(samps2$samps,1,quantile, probs = 0.025)
U2 <- apply(samps2$samps,1,quantile, probs = 0.975)




#### Approach 3: Using ARIMA method
D <- H[-c(1,n),]
R <- B[-c(1,n), -c(1,n)]

Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(n, x = 0.0001)), "dgTMatrix")

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q3,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q3,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u1 = 5,
  alpha1 = 0.1,
  u2 = 5,
  alpha2 = 0.1
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta1 = 0, # -2log(sigma)
  theta2 = 0
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))
# Plot of theta posterior
logpostsigma1 <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma1,plot(transparam,pdf_transparam,type='l'))
logpostsigma2 <- compute_pdf_and_cdf(quad$marginals[[2]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma2,plot(transparam,pdf_transparam,type='l'))
# Inference for W
samps3 <- sample_marginal(quad,1e03)
# Posterior mean
W3 <- apply(samps3$samps,1,mean)
L3 <- apply(samps3$samps,1,quantile, probs = 0.025)
U3 <- apply(samps3$samps,1,quantile, probs = 0.975)



#### Comparison
plot(y~x)

lines(predict(loess(W1~x, span = 2)), col = "red")
lines(predict(loess(U1~x, span = 2)), col = "red", lty=2)
lines(predict(loess(L1~x, span = 2)), col = "red", lty=2)


lines(predict(loess(W2~x, span = 2)), col = "blue")
lines(predict(loess(U2~x, span = 2)), col = "blue", lty=2)
lines(predict(loess(L2~x, span = 2)), col = "blue", lty=2)


lines(predict(loess(W3~x, span = 2)), col = "green")
lines(predict(loess(U3~x, span = 2)), col = "green", lty=2)
lines(predict(loess(L3~x, span = 2)), col = "green", lty=2)


lines(5*sin(0.5*x)~x, col = "black")





#### Comparison
plot(y~x)
lines(W1~x, col = "Red")
lines(U1~x, col = "Red", lty=2)
lines(L1~x, col = "Red", lty=2)


lines(W2~x, col = "blue")
lines(U2~x, col = "blue", lty=2)
lines(L2~x, col = "blue", lty=2)

lines(W3~x, col = "green")
lines(U3~x, col = "green", lty=2)
lines(L3~x, col = "green", lty=2)
lines(5*sin(0.5*x)~x, col = "black")











######################################################
######################################################
######## For Non-Gaussian Data #######################
######################################################

n <- 500
x <- seq(0.1,50, by = 0.1)
ylat <- log(5*sin(0.5*x) + 6)
y <- rpois(n,lambda = exp(ylat))
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")

H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 0.0001), "dgTMatrix")

compile("03_RW2Comparison.cpp")
dyn.load(dynlib("03_RW2Comparison"))



### Approach 1:
tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q1,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q1,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u = 2,
  alpha = 0.5
)
tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta = 0 # -2log(sigma)
)
ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "03_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,c(0))


# Inference for W
samps1 <- sample_marginal(quad,1e03)
# Posterior mean
W1 <- apply(samps1$samps,1,mean)
L1 <- apply(samps1$samps,1,quantile, probs = 0.025)
U1 <- apply(samps1$samps,1,quantile, probs = 0.975)


#### Comparison
plot(ylat~x, type = "l")
lines(W1~x, col = "Red")
lines(U1~x, col = "Red", lty=2)
lines(L1~x, col = "Red", lty=2)






### Approach 2:
Q2 <- t(H) %*% solve(B) %*% H
Q2 <- as(as.matrix(Q2 + Diagonal(n, x = 0.0001)), "dgTMatrix")

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q2,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q2,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u = 2,
  alpha = 0.5
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta = 0 # -2log(sigma)
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "03_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,0)

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))


# Inference for W
samps2 <- sample_marginal(quad,1e03)
# Posterior mean
W2 <- apply(samps2$samps,1,mean)
L2 <- apply(samps2$samps,1,quantile, probs = 0.025)
U2 <- apply(samps2$samps,1,quantile, probs = 0.975)


#### Comparison
plot(ylat~x, type = "l")
lines(W1~x, col = "Red")
lines(U1~x, col = "Red", lty=2)
lines(L1~x, col = "Red", lty=2)
lines(W2~x, col = "blue")
lines(U2~x, col = "blue", lty=2)
lines(L2~x, col = "blue", lty=2)





#### Approach 3:
D <- H[-c(1,n),]
R <- B[-c(1,n), -c(1,n)]

Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(n, x = 0.0001)), "dgTMatrix")

tmbdat <- list(
  # Design matrix
  X = X,
  # Penalty(Precision) matrix
  P = Q3,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(Q3,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # PC Prior params
  u = 2,
  alpha = 0.5
)

tmbparams <- list(
  W = rep(0, n), # W = c(U); U = B-Spline coefficients
  theta = 0 # -2log(sigma)
)

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "03_RW2Comparison",
  silent = TRUE
)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,0)

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))


# Inference for W
samps3 <- sample_marginal(quad,1e03)
# Posterior mean
W3 <- apply(samps3$samps,1,mean)
L3 <- apply(samps3$samps,1,quantile, probs = 0.025)
U3 <- apply(samps3$samps,1,quantile, probs = 0.975)



#### Comparison
plot(ylat~x, type = "l")
lines(W1~x, col = "Red")
lines(U1~x, col = "Red", lty=2)
lines(L1~x, col = "Red", lty=2)
lines(W2~x, col = "blue")
lines(U2~x, col = "blue", lty=2)
lines(L2~x, col = "blue", lty=2)
lines(W3~x, col = "green")
lines(U3~x, col = "green", lty=2)
lines(L3~x, col = "green", lty=2)






