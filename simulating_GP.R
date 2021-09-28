##### Simulation from conditional distribution of GP:


### Backward Prediction
### Inputs: sort(t_grid,x_grid) must be a equally spaced grid
## t_grid: a vector of locations to be simulated, must not be the same as x_grid
## x_grid: a vector of locations that are observed and inferred
## aghq: the aghq object obtained using the dataset {(y,x)} 
## GP: a string either be "ARIMA" or "RW2"
## n_samples has to be larger than the # col in gx
# (Note that this version requires t_grid to be placed in front of the x_grid)
simulation_GP_backward <- function(t_grid, x_grid, gx, GP, n_samples = 500){
  n <- length(x_grid)
  all_grid <- sort(c(t_grid, x_grid))
  all_d <- diff(all_grid)
  H <- compute_H_rue(all_d,n = length(all_grid))
  B <- compute_B(all_d,n = length(all_grid))
  A <- compute_A(all_d, n = length(all_grid))
  Q <- t(H) %*% solve(A) %*% H
  if(GP == "ARIMA"){
    D <- H[-c(1,length(all_grid)),]
    R <- B[-c(1,length(all_grid)), -c(1,length(all_grid))]
    Q <- t(D) %*% solve(R) %*% D
  }
  QAA <- Q[1:length(t_grid), 1:length(t_grid)]
  QAB <- Q[1:length(t_grid), -c(1:length(t_grid))]
  QBB <- Q[-c(1:length(t_grid)), -c(1:length(t_grid))]
  
  compute_conditional_mean <- function(gx){
    conditional_mean <- -solve(QAA)%*% QAB %*% (gx)
    conditional_mean
  }
  
  conditional_mean_list <- apply(gx, 2, compute_conditional_mean)
  
  simulated_gz <- function(cmu){
    Lt <- chol(QAA)
    z <- rnorm(n = length(t_grid))
    v <- solve(a = Lt, b = z)
    v + cmu
  }
  
  samples <- apply(conditional_mean_list, 2, simulated_gz)
  samples
}

### Forward Prediction
simulation_GP_forward <- function(t_grid, x_grid, gx, GP, n_samples = 500){
  n <- length(x_grid)
  all_grid <- sort(c(x_grid, t_grid))
  all_d <- diff(all_grid)
  H <- compute_H_rue(all_d,n = length(all_grid))
  B <- compute_B(all_d,n = length(all_grid))
  A <- compute_A(all_d, n = length(all_grid))
  Q <- t(H) %*% solve(A) %*% H
  if(GP == "ARIMA"){
    D <- H[-c(1,length(all_grid)),]
    R <- B[-c(1,length(all_grid)), -c(1,length(all_grid))]
    Q <- t(D) %*% solve(R) %*% D
  }
  QAA <- Q[(n+1):length(all_grid), (n+1):length(all_grid)]
  QAB <- Q[(n+1):length(all_grid), 1:n]
  QBB <- Q[1:n, 1:n]
  
  compute_conditional_mean <- function(gx){
    conditional_mean <- -solve(QAA)%*% QAB %*% (gx)
    conditional_mean
  }
  
  conditional_mean_list <- apply(gx, 2, compute_conditional_mean)
  
  simulated_gz <- function(cmu){
    Lt <- chol(QAA)
    z <- rnorm(n = length(t_grid))
    v <- solve(a = Lt, b = z)
    v + cmu
  }
  
  samples <- apply(conditional_mean_list, 2, simulated_gz)
  samples
}



#### Try Backward Extrapolation:
n <- 40
x <- seq(11,50, by = 1)
compute_g <- function(x){
  5*sin(0.5*x)
}
y <- compute_g(x) + rnorm(n, sd = 3)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")

H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))





################# RW2:
H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 0.0001), "dgTMatrix")

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
set.seed(123)
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))

### z_grid
z_grid <- seq(1,10, by = 1)
gx <- sample_marginal(quad, 1000)
gx <- gx$samps
gz <- simulation_GP_backward(z_grid, x, gx, "RW2")

mean_x <- apply(gx, 1, mean)
mean_z <- apply(gz, 1, mean)

plot(mean_z~z_grid, type = 'l', col = 'red', ylim = c(-80,80), xlim = c(0,50), xlab = "x", ylab = "y")
for (i in 1:1000) {
  lines(gz[,i]~z_grid, col = rgb(100, 0, 20, max = 255, alpha = 20, names = "grey"))
}
for (i in 1:1000) {
  lines(gx[,i]~x, col = rgb(0, 0, 20, max = 255, alpha = 20, names = "grey"))
}
lines(mean_x~x, col = "blue")

upper_pred <- apply(gz, 1, quantile, p = 0.975)
lower_pred <- apply(gz, 1, quantile, p = 0.025)

lines(upper_pred ~ z_grid, col = "orange")
lines(lower_pred ~ z_grid, col = "orange")
title(main = "RW2 Prediction with 400 samples")


mean(upper_pred-lower_pred)
mean((mean_z - compute_g(z_grid))^2)

### RW2 Performance:
## When n = 400:
## Prediction MSE: 127.6841
## Prediction Interval Width: 45.93148
## When n = 200:
## Prediction MSE: 362.9792
## Prediction Interval Width: 49.90343
## When n = 40:
## Prediction MSE: 238.9684
## Prediction Interval Width: 57.60811






################# ARIMA
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
set.seed(123)
quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))


### z_grid
z_grid <- seq(1,10, by = 1)
gx <- sample_marginal(quad, 500)
gx <- gx$samps
gz <- simulation_GP_backward(z_grid, x, gx, "ARIMA")

mean_x <- apply(gx, 1, mean)
mean_z <- apply(gz, 1, mean)

plot(mean_z~z_grid, type = 'l', col = 'red', ylim = c(-80,80), xlim = c(0,50))
for (i in 1:500) {
  lines(gz[,i]~z_grid, col = rgb(100, 0, 20, max = 255, alpha = 20, names = "grey"))
}
for (i in 1:500) {
  lines(gx[,i]~x, col = rgb(0, 0, 20, max = 255, alpha = 20, names = "grey"))
}
lines(mean_x~x, col = "blue")


upper_pred <- apply(gz, 1, quantile, p = 0.975)
lower_pred <- apply(gz, 1, quantile, p = 0.025)

lines(upper_pred ~ z_grid, col = "orange")
lines(lower_pred ~ z_grid, col = "orange")
title(main = "ARIMA Prediction with 400 samples")



mean(upper_pred-lower_pred)
mean((mean_z - compute_g(z_grid))^2)

### ARIMA Performance:
## When n = 400:
## Prediction MSE: 128.407
## Prediction Interval Width: 45.7697
## When n = 200:
## Prediction MSE: 471.1251
## Prediction Interval Width: 49.81479
## When n = 40:
## Prediction MSE: 236.7882
## Prediction Interval Width: 56.38347


###############################
###############################
########## Conclusion: ########
## When sample size is around 400, these methods give similar backward extrapolation results.
##############################
## The prediction results do not differ by much, regardless of the sample size.








#### Try Forward Extrapolation:
n <- 200
x <- seq(0.1,40, by = 0.2)
compute_g <- function(x){
  5*sin(0.5*x)
}
y <- compute_g(x) + rnorm(n, sd = 3)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")

H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))





################# RW2:
H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(n, x = 0.0001), "dgTMatrix")

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

### z_grid
z_grid <- seq(40.1,50, by = 0.1)
gx <- sample_marginal(quad, 500)
gx <- gx$samps
gz <- simulation_GP_forward(z_grid, x, gx, "RW2")

mean_x <- apply(gx, 1, mean)
mean_z <- apply(gz, 1, mean)

plot(mean_z~z_grid, type = 'l', col = 'red', ylim = c(-80,80), xlim = c(0,50), xlab = "x", ylab = "y")
for (i in 1:500) {
  lines(gz[,i]~z_grid, col = rgb(100, 0, 20, max = 255, alpha = 20, names = "grey"))
}
for (i in 1:500) {
  lines(gx[,i]~x, col = rgb(0, 0, 20, max = 255, alpha = 20, names = "grey"))
}
lines(mean_x~x, col = "blue")

upper_pred <- apply(gz, 1, quantile, p = 0.975)
lower_pred <- apply(gz, 1, quantile, p = 0.025)

lines(upper_pred ~ z_grid, col = "orange")
lines(lower_pred ~ z_grid, col = "orange")
title(main = "RW2 Prediction with 400 samples")

### RW2 Performance:
## When n = 400:
## Prediction MSE: 210.3667
## Prediction Interval Width: 42.21388






################# ARIMA
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


### z_grid
z_grid <- seq(40.1,50, by = 0.1)
gx <- sample_marginal(quad, 500)
gx <- gx$samps
gz <- simulation_GP_forward(z_grid, x, gx, "ARIMA")

mean_x <- apply(gx, 1, mean)
mean_z <- apply(gz, 1, mean)

plot(mean_z~z_grid, type = 'l', col = 'red', ylim = c(-80,80), xlim = c(0,50), xlab = "x", ylab = "y")
for (i in 1:500) {
  lines(gz[,i]~z_grid, col = rgb(100, 0, 20, max = 255, alpha = 20, names = "grey"))
}
for (i in 1:500) {
  lines(gx[,i]~x, col = rgb(0, 0, 20, max = 255, alpha = 20, names = "grey"))
}
lines(mean_x~x, col = "blue")


upper_pred <- apply(gz, 1, quantile, p = 0.975)
lower_pred <- apply(gz, 1, quantile, p = 0.025)

lines(upper_pred ~ z_grid, col = "orange")
lines(lower_pred ~ z_grid, col = "orange")
title(main = "ARIMA Prediction with 400 samples")

### ARIMA Performance:
## When n = 400:
## Prediction MSE: 243.2342
## Prediction Interval Width: 45.42289















