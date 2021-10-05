library(aghq)
library(TMB)
library(Matrix)
library(tidyverse)
library(foreach)
library(parallel)


cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)


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
compile(file = "02_RW2Comparison.cpp")





### B: number of replications
### dis: distance between knots x_i and x_{i+1} (integer that will be multipied by 0.5)
### k_n: number of measurements at one location
### M: number of locations in the high resolution grids
### nsamp: number of samples to be drawn from the posterior distribution
### return: dataframe with four columns, two column for each method
replicate_for_summary_once <- function(dis = 20, k_n = 1, n_samp = 2000){
  dyn.load(dynlib("02_RW2Comparison"))
  z <- seq(0.5,100,0.5)
  x <- seq(1,100, dis*0.5)
  z_grid <- z[!z %in% x]
  compute_g <- function(x){
    5*sin(0.1*x)
  }
  k <- length(x)
  n <- k*k_n
  y <- compute_g(rep(x, each = k_n)) + rnorm(n, sd = sqrt(3))
  d <- diff(x)
  X <- as(matrix(0,nrow = n, ncol = k), "dgTMatrix")
  for (i in 1:k) {
    X[,i] <- c(rep(0, (i-1)*k_n), rep(1,k_n), rep(0, ((n-k_n)-(i-1)*k_n)))
  }
  
  
  H <- compute_H_rue(d,n = length(x))
  B <- compute_B(d,n = length(x))
  H <- compute_H_rue(d,n = length(x))
  B <- compute_B(d,n = length(x))
  A <- compute_A(d, n = length(x))
  
  
  ###### RW2:
  Q1 <- t(H) %*% solve(A) %*% H
  Q1 <- as(Q1 + Diagonal(length(x), x = 0.00001), "dgTMatrix")
  
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
    W = rep(0, length(x)), # W = c(U); U = B-Spline coefficients
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
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))
  gx <- sample_marginal(quad, n_samp)
  gx <- gx$samps
  gz_list <- Interpolation_vec_v1(t = z_grid, x, gx, "RW2")
  mean_x <- apply(gx, 1, mean)
  sample_path <- Matrix(0,nrow = length(z), ncol = n_samp)
  for (i in 1:length(gz_list)) {
    sample_path[,i] <- rbind(gz_list[[i]], matrix(gx[,i],ncol = 1))
  }
  sample_path_rw2 <- as.tibble(as.matrix(sample_path))
  sample_path_rw2$locations <- c(z_grid,x)
  sample_path_rw2 <- arrange(sample_path_rw2, by = locations)
  mean_y_rw2 <- apply(sample_path_rw2, 1, mean)
  upper_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.95)
  lower_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.05)
  rIAE_rw2 <- sqrt(mean(abs(compute_g(z) - mean_y_rw2)))
  MCI_rw2 <- mean(upper_y_rw2 - lower_y_rw2)
  CR_rw2 <- mean(upper_y_rw2 >= compute_g(z) & lower_y_rw2 <= compute_g(z))
  ######
  
  D <- H[-c(1,length(x)),]
  R <- B[-c(1,length(x)), -c(1,length(x))]
  Q3 <- t(D) %*% solve(R) %*% D
  Q3 <- as(as.matrix(Q3 + Diagonal(length(x), x = 0.00001)), "dgTMatrix")
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
  quad <- aghq::marginal_laplace_tmb(ff,7,c(0,0))
  gx <- sample_marginal(quad, n_samp)
  gx <- gx$samps
  gz_list <- Interpolation_vec_v1(t = z_grid, x, gx, "ARIMA")
  gz <- Matrix(0,nrow = length(z_grid), ncol = n_samp)
  for (i in 1:length(gz_list)) {
    gz[,i] <- gz_list[[i]]
  }
  sample_path <- Matrix(0,nrow = length(z), ncol = n_samp)
  for (i in 1:length(gz_list)) {
    sample_path[,i] <- rbind(gz_list[[i]], matrix(gx[,i],ncol = 1))
  }
  mean_x <- apply(gx, 1, mean)
  sample_path_arima <- as.tibble(as.matrix(sample_path))
  sample_path_arima$locations <- c(z_grid,x)
  sample_path_arima <- arrange(sample_path_arima, by = locations)
  mean_y_arima <- apply(sample_path_arima, 1, mean)
  upper_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.95)
  lower_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.05)
  rIAE_arima <- sqrt(mean(abs(compute_g(z) - mean_y_arima)))
  MCI_arima <- mean(upper_y_arima - lower_y_arima)
  CR_arima <- mean(upper_y_arima >= compute_g(z) & lower_y_arima <= compute_g(z))
  
  result <- tibble(rIAE_RW2 = rIAE_rw2, rIAE_ARIMA = rIAE_arima, MCI_RW2 = MCI_rw2, MCI_ARIMA = MCI_arima, CR_RW2 = CR_rw2, CR_ARIMA = CR_arima)
  
  #### First deriv:
  sample_path_rw2_1stDeriv <- apply(sample_path_rw2, 2, diff)
  sample_path_rw2_1stDeriv[,(n_samp + 1)] <- unlist(sample_path_rw2[,(n_samp + 1)])[-1]
  sample_path_rw2_1stDeriv_mean <- apply(sample_path_rw2_1stDeriv, 1, mean)
  sample_path_rw2_1stDeriv_upper <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.95)
  sample_path_rw2_1stDeriv_lower <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.05)
  CR_rw2_1st <- mean(sample_path_rw2_1stDeriv_upper >= diff(compute_g(z)) & sample_path_rw2_1stDeriv_lower <= diff(compute_g(z)) )
  
  
  sample_path_arima_1stDeriv <- apply(sample_path_arima, 2, diff)
  sample_path_arima_1stDeriv[,(n_samp + 1)] <- unlist(sample_path_arima[,(n_samp + 1)])[-1]
  sample_path_arima_1stDeriv_mean <- apply(sample_path_arima_1stDeriv, 1, mean)
  sample_path_arima_1stDeriv_upper <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.95)
  sample_path_arima_1stDeriv_lower <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.05)
  CR_arima_1st <- mean(sample_path_arima_1stDeriv_upper >= diff(compute_g(z)) & sample_path_arima_1stDeriv_lower <= diff(compute_g(z)) )
  
  
  result_1st <- tibble(rIAE_RW2_1st = sqrt(mean(abs(diff(compute_g(z)) - sample_path_rw2_1stDeriv_mean))), rIAE_ARIMA_1st = sqrt(mean(abs(diff(compute_g(z)) - sample_path_arima_1stDeriv_mean))), MCI_RW2_1st = mean(sample_path_rw2_1stDeriv_upper - sample_path_rw2_1stDeriv_lower), MCI_ARIMA_1st = mean(sample_path_arima_1stDeriv_upper - sample_path_arima_1stDeriv_lower), CR_RW2_1st = CR_rw2_1st, CR_ARIMA_1st = CR_arima_1st)
  result <- cbind(result, result_1st)
  
  ### Second deriv
  sample_path_rw2_2ndDeriv <- apply(sample_path_rw2, 2, diff, differences = 2)
  sample_path_rw2_2ndDeriv[,(n_samp + 1)] <- unlist(sample_path_rw2[,(n_samp + 1)])[-c(1,2)]
  sample_path_rw2_2ndDeriv_mean <- apply(sample_path_rw2_2ndDeriv, 1, mean)
  sample_path_rw2_2ndDeriv_upper <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.95)
  sample_path_rw2_2ndDeriv_lower <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.05)
  CR_rw2_2nd <- mean(sample_path_rw2_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_rw2_2ndDeriv_lower <= diff(compute_g(z), differences = 2))
  
  sample_path_arima_2ndDeriv <- apply(sample_path_arima, 2, diff, differences = 2)
  sample_path_arima_2ndDeriv[,(n_samp + 1)] <- unlist(sample_path_arima[,(n_samp + 1)])[-c(1,2)]
  sample_path_arima_2ndDeriv_mean <- apply(sample_path_arima_2ndDeriv, 1, mean)
  sample_path_arima_2ndDeriv_upper <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.95)
  sample_path_arima_2ndDeriv_lower <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.05)
  CR_arima_2nd <- mean(sample_path_arima_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_arima_2ndDeriv_lower <= diff(compute_g(z), differences = 2))
  
  result_2nd <- tibble(rIAE_RW2_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2) - sample_path_rw2_2ndDeriv_mean))), rIAE_ARIMA_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2) - sample_path_arima_2ndDeriv_mean))), MCI_RW2_2nd = mean(sample_path_rw2_2ndDeriv_upper - sample_path_rw2_2ndDeriv_lower), MCI_ARIMA_2nd = mean(sample_path_arima_2ndDeriv_upper - sample_path_arima_2ndDeriv_lower), CR_RW2_2nd = CR_rw2_2nd, CR_ARIMA_2nd = CR_arima_2nd)
  result <- cbind(result, result_2nd)
  result
  
}




set.seed(123)
time_begin <- Sys.time()
result <- replicate_for_summary_once()
# time_begin <- Sys.time()
# result <- rbind(result, replicate_for_summary_once())
# time_end <- Sys.time()
# time_end - time_begin
## takes around 20 seconds per iteration

# for (i in 1:99) {
#   result <- rbind(result, replicate_for_summary_once())
# }
# time_end <- Sys.time()
# time_end - time_begin
# ### five iterations take 3.602 mins
# save(file = "resultCase1.rda", result)


time_begin <- Sys.time()
result  <- foreach(i = 1:100, .combine='rbind', .packages = c("aghq", "TMB", "Matrix", "tidyverse")) %dopar% {
  replicate_for_summary_once()
}
time_end <- Sys.time()
time_end - time_begin
save(file = "resultCase1.rda", result)

### five iterations take 2.182191 mins

# 
# result_n50 <- replicate_for_summary_once(dis = 4, k_n = 1, n_samp = 2000)
# for (i in 1:99) {
#   result_n50 <- rbind(result_n50, replicate_for_summary_once(dis = 4, k_n = 1, n_samp = 2000))
# }
# save(file = "resultCase2.rda", result_n50)

time_begin <- Sys.time()
result_n50  <- foreach(i = 1:100, .combine='rbind', .packages = c("aghq", "TMB", "Matrix", "tidyverse")) %dopar% {
  replicate_for_summary_once(dis = 4, k_n = 1, n_samp = 2000)
}
time_end <- Sys.time()
time_end - time_begin
save(file = "resultCase2.rda", result_n50)



# result_n50 <- foreach(i = 1:99, .combine = 'rbind', packages = c("Matrix", 'aghq')) %dopar% {
#   replicate_for_summary_once(dis = 4, k_n = 1, n_samp = 2000)
# }

time_begin <- Sys.time()
result_n50_repeat5  <- foreach(i = 1:100, .combine='rbind', .packages = c("aghq", "TMB", "Matrix", "tidyverse")) %dopar% {
  replicate_for_summary_once(dis = 20, k_n = 5, n_samp = 2000)
}
time_end <- Sys.time()
time_end - time_begin
save(file = "resultCase3.rda", result_n50_repeat5)



parallel::stopCluster(cl)





