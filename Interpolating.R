##### Interpolating GP simulation:

n_samp <- 8000

#################### Interpolate a single point

Interpolation_point_v1 <- function(t, x_grid, gx, GP, n_samples = 500){
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

  C <- solve(Q + Diagonal(n, x = 0.0001))
  C_trans <- C[-t_location, -t_location]
  C_trans <- bdiag(C_trans, C[t_location, t_location])
  C_trans[1:(n-1), n] <- C[, t_location][-t_location]
  C_trans[n, 1:(n-1)] <- C[t_location, ][-t_location]
  Q_trans <- solve(C_trans)
  QAA <- Q_trans[n,n]
  QAB <- Q_trans[n,1:(n-1)]
  QBB <- Q[1:(n-1), 1:(n-1)]
  
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
  
  samples <- apply(matrix(conditional_mean_list, nrow = 1), 2, simulated_gz)
  samples
}

n <- 50
x <- seq(0,50, by = 1)
x <- x[-31]
compute_g <- function(x){
  5*sin(0.5*x)
}
y <- compute_g(x) + rnorm(n, sd = 3)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")

H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))


#################### RW2

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
gx <- sample_marginal(quad, 2000)
gx <- gx$samps
z <- 30

gz <- Interpolation_point_v1(t = z, x_grid = x, gx = gx, GP = "RW2", n_samples = 2000)


mean_x <- apply(gx, 1, mean)
upper_x <- apply(gx, 1, quantile, p = 0.975)
lower_x <- apply(gx, 1, quantile, p = 0.025)

mean_z <- mean(gz)


plot(mean_x~x, type = 'l', col = 'red', xlab = "x", ylab = "y", ylim = c(-10,10))
lines(compute_g(seq(0,50,by = 0.1))~seq(0,50,by = 0.1), type = 'l', col = "black")
lines(c(upper_x[1:30], quantile(gz, 0.975), upper_x[31:50]) ~ c(x[1:30], z, x[31:50]), col = 'orange')
lines(c(lower_x[1:30], quantile(gz, 0.025), lower_x[31:50]) ~ c(x[1:30], z, x[31:50]), col = 'orange')
points(mean_z~z, col = 'blue')

### Overall mean CI width:
mean(upper_x - lower_x) # 6.163772
### PI width:
quantile(gz, 0.975) - quantile(gz, 0.025) # 6.794209

### Overall MSE using Posterior Mean:
mean((mean_x - compute_g(x))^2) # 0.9147514

### MSE using each posterior sample:
compute_sample_MSE <- function(gx){mean((gx - compute_g(x))^2)}
hist(apply(gx, 2, compute_sample_MSE))
mean(apply(gx, 2, compute_sample_MSE)) # 3.386134

### Point prediction MSE using Posterior Mean:
(mean_z - compute_g(z))^2 # 0.05308926

### Point prediction MSE using each Posterior sample:
hist((gz - compute_g(z))^2)
mean((gz - compute_g(z))^2) # 3.001648






################################### ARIMA

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
gx <- sample_marginal(quad, 2000)
gx <- gx$samps
z <- 30

gz <- Interpolation_point_v1(t = z, x_grid = x, gx = gx, GP = "ARIMA", n_samples = 2000)


mean_x <- apply(gx, 1, mean)
upper_x <- apply(gx, 1, quantile, p = 0.975)
lower_x <- apply(gx, 1, quantile, p = 0.025)

mean_z <- mean(gz)


plot(mean_x~x, type = 'l', col = 'red', xlab = "x", ylab = "y", ylim = c(-10,10))
lines(compute_g(seq(0,50,by = 0.1))~seq(0,50,by = 0.1), type = 'l', col = "black")
lines(c(upper_x[1:30], quantile(gz, 0.975), upper_x[31:50]) ~ c(x[1:30], z, x[31:50]), col = 'orange')
lines(c(lower_x[1:30], quantile(gz, 0.025), lower_x[31:50]) ~ c(x[1:30], z, x[31:50]), col = 'orange')
points(mean_z~z, col = 'blue')

### Overall mean CI width:
mean(upper_x - lower_x) # 6.001914
### PI width:
quantile(gz, 0.975) - quantile(gz, 0.025) # 6.592026

### Overall MSE using Posterior Mean:
mean((mean_x - compute_g(x))^2) # 0.9020369

### MSE using each posterior sample:
compute_sample_MSE <- function(gx){mean((gx - compute_g(x))^2)}
hist(apply(gx, 2, compute_sample_MSE))
mean(apply(gx, 2, compute_sample_MSE)) # 3.272603

### Point prediction MSE using Posterior Mean:
(mean_z - compute_g(z))^2 # 0.1014458

### Point prediction MSE using each Posterior sample:
hist((gz - compute_g(z))^2)
mean((gz - compute_g(z))^2) # 3.058983
































###### Interpolate a whole vector:

### Switch Position function:

switch_matrix_once <- function(M, from, to){
  L <- diag(1, nrow = nrow(M), ncol = nrow(M))
  original_row <- L[to,]
  L[to,] <- L[from,]
  L[from,] <- original_row
  L%*%M%*%t(L)
}

### Test
Q
switch_matrix_once(Q, 2,3)


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


### Test
Q
switch_matrix(c(2,4,6,8),c(1,3,5,7), Q)

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
  Q_trans <- forceSymmetric(solve(C_trans))
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











### Simulation with sparse data: (n = 50)
n <- 50
x <- seq(1,99, by = 2)
compute_g <- function(x){
  5*sin(0.1*x)
}
y <- compute_g(x) + rnorm(n, sd = 3)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")


H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))


###### RW2:
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
z_grid <- seq(2,100, by = 2)
gx <- sample_marginal(quad, n_samp)
gx <- gx$samps
gz_list <- Interpolation_vec_v1(t = z_grid, x, gx, "RW2")
gz <- Matrix(0,nrow = 50, ncol = n_samp)
for (i in 1:length(gz_list)) {
  gz[,i] <- gz_list[[i]]
}

mean_x <- apply(gx, 1, mean)
upper_x <- apply(gx, 1, quantile, p = 0.975)
lower_x <- apply(gx, 1, quantile, p = 0.025)

mean_z <- apply(gz, 1, mean)
upper_z <- apply(gz, 1, quantile, p = 0.975)
lower_z <- apply(gz, 1, quantile, p = 0.025)

### Observed grid:
plot(mean_x ~ x, type = 'l', xlab = "observed grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_x ~ x, lty = 2, col = 'orange')
lines(lower_x ~ x, lty = 2, col = 'orange')
lines(compute_g(x) ~ x, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gx[,i] ~ x, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "RW2 for sparse covariate: Inference")
cat(paste("MSE = ", round(mean((mean_x - compute_g(x))^2), 3))) ## 0.644
cat(paste("CI-width = ", round(mean((upper_x - lower_x)), 3))) ## 5.058



### Interpolated grid:
plot(mean_z ~ z_grid, type = 'l', xlab = "interpolated grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_z ~ z_grid, lty = 2, col = 'orange')
lines(lower_z ~ z_grid, lty = 2, col = 'orange')
lines(compute_g(z_grid) ~ z_grid, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gz[,i] ~ z_grid, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

title(main = "RW2 for sparse covariate: Interpolation")
cat(paste("MSE = ", round(mean((mean_z - compute_g(z_grid))^2), 3))) ## 0.609
cat(paste("CI-width = ", round(mean((upper_z - lower_z)), 3))) ## 5.199





###### ARIMA:
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
z_grid <- seq(2,100, by = 2)
gx <- sample_marginal(quad, n_samp)
gx <- gx$samps
gz_list <- Interpolation_vec_v1(t = z_grid, x, gx, "ARIMA")
gz <- Matrix(0,nrow = 50, ncol = n_samp)
for (i in 1:length(gz_list)) {
  gz[,i] <- gz_list[[i]]
}

mean_x <- apply(gx, 1, mean)
upper_x <- apply(gx, 1, quantile, p = 0.975)
lower_x <- apply(gx, 1, quantile, p = 0.025)

mean_z <- apply(gz, 1, mean)
upper_z <- apply(gz, 1, quantile, p = 0.975)
lower_z <- apply(gz, 1, quantile, p = 0.025)

### Observed grid:
plot(mean_x ~ x, type = 'l', xlab = "observed grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_x ~ x, lty = 2, col = 'orange')
lines(lower_x ~ x, lty = 2, col = 'orange')
lines(compute_g(x) ~ x, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gx[,i] ~ x, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "ARIMA for sparse covariate: Inference")
cat(paste("MSE = ", round(mean((mean_x - compute_g(x))^2), 3))) ## 0.638
cat(paste("CI-width = ", round(mean((upper_x - lower_x)), 3))) ## 5.033



### Interpolated grid:
plot(mean_z ~ z_grid, type = 'l', xlab = "interpolated grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_z ~ z_grid, lty = 2, col = 'orange')
lines(lower_z ~ z_grid, lty = 2, col = 'orange')
lines(compute_g(z_grid) ~ z_grid, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gz[,i] ~ z_grid, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

title(main = "ARIMA for sparse covariate: Interpolation")
cat(paste("MSE = ", round(mean((mean_z - compute_g(z_grid))^2), 3))) ## 0.609
cat(paste("CI-width = ", round(mean((upper_z - lower_z)), 3))) ## 5.199










### Simulation with very sparse data: (n = 20)
n <- 20
x <- seq(5,100, by = 5)
compute_g <- function(x){
  5*sin(0.1*x)
}
y <- compute_g(x) + rnorm(n, sd = 3)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")


H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))


###### RW2:
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
z_grid <- seq(2,100, by = 5)
gx <- sample_marginal(quad, n_samp)
gx <- gx$samps
gz_list <- Interpolation_vec_v1(t = z_grid, x, gx, "RW2")
gz <- Matrix(0,nrow = n, ncol = n_samp)
for (i in 1:length(gz_list)) {
  gz[,i] <- gz_list[[i]]
}

mean_x <- apply(gx, 1, mean)
upper_x <- apply(gx, 1, quantile, p = 0.975)
lower_x <- apply(gx, 1, quantile, p = 0.025)

mean_z <- apply(gz, 1, mean)
upper_z <- apply(gz, 1, quantile, p = 0.975)
lower_z <- apply(gz, 1, quantile, p = 0.025)

### Observed grid:
plot(mean_x ~ x, type = 'l', xlab = "observed grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_x ~ x, lty = 2, col = 'orange')
lines(lower_x ~ x, lty = 2, col = 'orange')
lines(compute_g(x) ~ x, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gx[,i] ~ x, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "RW2 for very sparse covariate: Inference")
cat(paste("MSE = ", round(mean((mean_x - compute_g(x))^2), 3))) ## 4.994
cat(paste("CI-width = ", round(mean((upper_x - lower_x)), 3))) ## 5.146



### Interpolated grid:
plot(mean_z ~ z_grid, type = 'l', xlab = "interpolated grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_z ~ z_grid, lty = 2, col = 'orange')
lines(lower_z ~ z_grid, lty = 2, col = 'orange')
lines(compute_g(z_grid) ~ z_grid, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gz[,i] ~ z_grid, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

title(main = "RW2 for very sparse covariate: Interpolation")
cat(paste("MSE = ", round(mean((mean_z - compute_g(z_grid))^2), 3))) ## 4.816
cat(paste("CI-width = ", round(mean((upper_z - lower_z)), 3))) ## 8.696





###### ARIMA:
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
z_grid <- seq(2,100, by = 5)
gx <- sample_marginal(quad, n_samp)
gx <- gx$samps
gz_list <- Interpolation_vec_v1(t = z_grid, x, gx, "ARIMA")
gz <- Matrix(0,nrow = n, ncol = n_samp)
for (i in 1:length(gz_list)) {
  gz[,i] <- gz_list[[i]]
}

mean_x <- apply(gx, 1, mean)
upper_x <- apply(gx, 1, quantile, p = 0.975)
lower_x <- apply(gx, 1, quantile, p = 0.025)

mean_z <- apply(gz, 1, mean)
upper_z <- apply(gz, 1, quantile, p = 0.975)
lower_z <- apply(gz, 1, quantile, p = 0.025)

### Observed grid:
plot(mean_x ~ x, type = 'l', xlab = "observed grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_x ~ x, lty = 2, col = 'orange')
lines(lower_x ~ x, lty = 2, col = 'orange')
lines(compute_g(x) ~ x, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gx[,i] ~ x, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "ARIMA for very sparse covariate: Inference")
cat(paste("MSE = ", round(mean((mean_x - compute_g(x))^2), 3))) ## 4.726
cat(paste("CI-width = ", round(mean((upper_x - lower_x)), 3))) ## 5.209



### Interpolated grid:
plot(mean_z ~ z_grid, type = 'l', xlab = "interpolated grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_z ~ z_grid, lty = 2, col = 'orange')
lines(lower_z ~ z_grid, lty = 2, col = 'orange')
lines(compute_g(z_grid) ~ z_grid, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gz[,i] ~ z_grid, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

title(main = "ARIMA for very sparse covariate: Interpolation")
cat(paste("MSE = ", round(mean((mean_z - compute_g(z_grid))^2), 3))) ## 4.57
cat(paste("CI-width = ", round(mean((upper_z - lower_z)), 3))) ## 7.379








### Simulation with dense data: (n = 100)
n <- 100
x <- seq(0.5,99.5, by = 1)
compute_g <- function(x){
  5*sin(0.1*x)
}
y <- compute_g(x) + rnorm(n, sd = 3)
d <- diff(x)
X <- as(as.matrix(Diagonal(n)), "dgTMatrix")


H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
H <- compute_H_rue(d,n = length(x))
B <- compute_B(d,n = length(x))
A <- compute_A(d, n = length(x))


###### RW2:
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
z_grid <- seq(0,99, by = 1)
gx <- sample_marginal(quad, n_samp)
gx <- gx$samps
gz_list <- Interpolation_vec_v1(t = z_grid, x, gx, "RW2")
gz <- Matrix(0,nrow = n, ncol = n_samp)
for (i in 1:length(gz_list)) {
  gz[,i] <- gz_list[[i]]
}

mean_x <- apply(gx, 1, mean)
upper_x <- apply(gx, 1, quantile, p = 0.975)
lower_x <- apply(gx, 1, quantile, p = 0.025)

mean_z <- apply(gz, 1, mean)
upper_z <- apply(gz, 1, quantile, p = 0.975)
lower_z <- apply(gz, 1, quantile, p = 0.025)

### Observed grid:
plot(mean_x ~ x, type = 'l', xlab = "observed grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_x ~ x, lty = 2, col = 'orange')
lines(lower_x ~ x, lty = 2, col = 'orange')
lines(compute_g(x) ~ x, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gx[,i] ~ x, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "RW2 for dense covariate: Inference")
cat(paste("MSE = ", round(mean((mean_x - compute_g(x))^2), 3))) ## 0.543
cat(paste("CI-width = ", round(mean((upper_x - lower_x)), 3))) ## 3.802



### Interpolated grid:
plot(mean_z ~ z_grid, type = 'l', xlab = "interpolated grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_z ~ z_grid, lty = 2, col = 'orange')
lines(lower_z ~ z_grid, lty = 2, col = 'orange')
lines(compute_g(z_grid) ~ z_grid, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gz[,i] ~ z_grid, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

title(main = "RW2 for dense covariate: Interpolation")
cat(paste("MSE = ", round(mean((mean_z - compute_g(z_grid))^2), 3))) ## 0.54
cat(paste("CI-width = ", round(mean((upper_z - lower_z)), 3))) ## 3.849





###### ARIMA:
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
z_grid <- seq(0,99, by = 1)
gx <- sample_marginal(quad, n_samp)
gx <- gx$samps
gz_list <- Interpolation_vec_v1(t = z_grid, x, gx, "ARIMA")
gz <- Matrix(0,nrow = n, ncol = n_samp)
for (i in 1:length(gz_list)) {
  gz[,i] <- gz_list[[i]]
}

mean_x <- apply(gx, 1, mean)
upper_x <- apply(gx, 1, quantile, p = 0.975)
lower_x <- apply(gx, 1, quantile, p = 0.025)

mean_z <- apply(gz, 1, mean)
upper_z <- apply(gz, 1, quantile, p = 0.975)
lower_z <- apply(gz, 1, quantile, p = 0.025)

### Observed grid:
plot(mean_x ~ x, type = 'l', xlab = "observed grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_x ~ x, lty = 2, col = 'orange')
lines(lower_x ~ x, lty = 2, col = 'orange')
lines(compute_g(x) ~ x, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gx[,i] ~ x, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "ARIMA for dense covariate: Inference")
cat(paste("MSE = ", round(mean((mean_x - compute_g(x))^2), 3))) ## 0.426
cat(paste("CI-width = ", round(mean((upper_x - lower_x)), 3))) ## 3.752



### Interpolated grid:
plot(mean_z ~ z_grid, type = 'l', xlab = "interpolated grid", ylab = "y", col = "red", ylim = c(-8,8), lty = 1)
lines(upper_z ~ z_grid, lty = 2, col = 'orange')
lines(lower_z ~ z_grid, lty = 2, col = 'orange')
lines(compute_g(z_grid) ~ z_grid, lty = 3, col = 'black')
for (i in sample.int(n_samp,50)) {
  lines(gz[,i] ~ z_grid, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}

title(main = "ARIMA for dense covariate: Interpolation")
cat(paste("MSE = ", round(mean((mean_z - compute_g(z_grid))^2), 3))) ## 0.448
cat(paste("CI-width = ", round(mean((upper_z - lower_z)), 3))) ## 3.782




























