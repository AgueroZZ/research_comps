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



# Inference for W
samps1 <- sample_marginal(quad,1e03)
resolution = 0.001
location_interest <- seq(a,b, by = resolution)
fitted_func <- construct_fitted(weights = samps1$samps[,1])
plot(fitted_func~location_interest, col = "blue", type = "l")

# Construct samples for g(.):
sample_func <- apply(samps1$samps, MARGIN = 2, FUN = construct_fitted)




# Posterior mean
W1 <- apply(samps1$samps,1,mean)
mean_func <- apply(sample_func,1,mean)

### Plot:
plot(mean_func~location_interest, col = "blue", type = "l")
for (i in sample.int(1000,50)) {
  lines(sample_func[,i] ~ location_interest, col = "grey")
}




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




# Posterior mean
W2 <- apply(samps2$samps,1,mean)
mean_func <- apply(sample_func,1,mean)

### Plot:
plot(mean_func~location_interest, col = "blue", type = "l")
for (i in sample.int(1000,50)) {
  lines(sample_func[,i] ~ location_interest, col = "grey")
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
### Plot:
plot(mean_func_cubic~location_interest, col = "blue", type = "l")
for (i in sample.int(1000,50)) {
  lines(sample_func_cubic[,i] ~ location_interest, col = "grey")
}


######## Linear way:
sample_func <- apply(samps3$samps, MARGIN = 2, FUN = construct_fitted)
mean_func <- apply(sample_func,1,mean)
### Plot:
plot(mean_func~location_interest, col = "blue", type = "l")
for (i in sample.int(1000,50)) {
  lines(sample_func[,i] ~ location_interest, col = "grey")
}







#### Comparison
plot(y~x)
lines(W1~x, col = "Red")
lines(W2~x, col = "blue")
lines(W3~x, col = "green")
lines(5*sin(0.5*x)~x, col = "black")





























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






