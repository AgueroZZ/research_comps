library(aghq)
library(TMB)
library(Matrix)
library(tidyverse)
library(foreach)
library(parallel)


cl <- parallel::makeCluster(8)
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
construct_A <- function(all_grid, x_indx){
  A <- matrix(0, nrow = length(x_indx), ncol = length(all_grid))
  for (i in 1:nrow(A)) {
    A[i,x_indx[i]] <- 1
  }
  A
}



compile(file = "02_RW2Comparison.cpp")




### B: number of replications
### dis: distance between knots x_i and x_{i+1} (integer that will be multipied by 0.5)
### k_n: number of measurements at one location
### M: number of locations in the high resolution grids
### nsamp: number of samples to be drawn from the posterior distribution
### return: dataframe with four columns, two column for each method
dis <- 20
z <- round(seq(0.5,100, 0.5),2)
x <- seq(1, 100, dis*0.5)
x_indx <- which(z %in% x)
construct_X <- construct_A(z, x_indx)
d <- diff(z)
H <- compute_H_rue(d,n = length(z))
B <- compute_B(d,n = length(z))
A <- compute_A(d, n = length(z))

###### RW2:
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(length(z), x = 0.00001), "dgTMatrix")
###### ARIMA:
D <- H[-c(1,length(z)),]
R <- B[-c(1,length(z)), -c(1,length(z))]
Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(length(z), x = 0.00001)), "dgTMatrix")


replicate_for_summary_once <- function(z, k_n = 1, x_indx, n_samp = 2000, Q1, Q3){
  dyn.load(dynlib("02_RW2Comparison"))
  x <- z[x_indx]
  compute_g <- function(x){
    5*sin(0.1*x)
  }
  k <- length(x)
  n <- k * k_n
  y <- compute_g(rep(x, each = k_n)) + rnorm(n, sd = sqrt(3))
  d <- diff(z)
  h <- mean(d)
  m <- length(z[-x_indx])
  X <- as(construct_X, "dgTMatrix")
  X <- X[rep(c(1:k), each = k_n),]
  
  # H <- compute_H_rue(d,n = length(z))
  # B <- compute_B(d,n = length(z))
  # A <- compute_A(d, n = length(z))
  # 
  
  # ###### RW2:
  # Q1 <- t(H) %*% solve(A) %*% H
  # Q1 <- as(Q1 + Diagonal(length(z), x = 0.00001), "dgTMatrix")
  # 
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
    W = rep(0, length(z)), # W = c(U); U = B-Spline coefficients
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
  gz <- sample_marginal(quad, n_samp)
  sample_path_rw2 <- gz$samps
  sample_path_rw2 <- as.tibble(as.matrix(sample_path_rw2))
  mean_y_rw2 <- apply(sample_path_rw2, 1, mean)
  upper_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.95)
  lower_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.05)
  rIAE_rw2 <- sqrt(mean(abs(compute_g(z) - mean_y_rw2)))
  MCI_rw2 <- mean(upper_y_rw2 - lower_y_rw2)
  CR_rw2 <- mean(upper_y_rw2 >= compute_g(z) & lower_y_rw2 <= compute_g(z))
  ######
  # D <- H[-c(1,length(z)),]
  # R <- B[-c(1,length(z)), -c(1,length(z))]
  # Q3 <- t(D) %*% solve(R) %*% D
  # Q3 <- as(as.matrix(Q3 + Diagonal(length(z), x = 0.00001)), "dgTMatrix")
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
  gz <- sample_marginal(quad, n_samp)
  sample_path_arima <- gz$samps
  mean_z <- apply(sample_path_arima, 1, mean)
  sample_path_arima <- as.tibble(as.matrix(sample_path_arima))

  mean_y_arima <- apply(sample_path_arima, 1, mean)
  upper_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.95)
  lower_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.05)
  rIAE_arima <- sqrt(mean(abs(compute_g(z) - mean_y_arima)))
  MCI_arima <- mean(upper_y_arima - lower_y_arima)
  CR_arima <- mean(upper_y_arima >= compute_g(z) & lower_y_arima <= compute_g(z))
  
  result <- tibble(rIAE_RW2 = rIAE_rw2, rIAE_ARIMA = rIAE_arima, MCI_RW2 = MCI_rw2, MCI_ARIMA = MCI_arima, CR_RW2 = CR_rw2, CR_ARIMA = CR_arima)
  
  #### First deriv:
  sample_path_rw2_1stDeriv <- apply(sample_path_rw2, 2, diff)
  sample_path_rw2_1stDeriv_mean <- apply(sample_path_rw2_1stDeriv, 1, mean)
  sample_path_rw2_1stDeriv_upper <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.95)
  sample_path_rw2_1stDeriv_lower <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.05)
  CR_rw2_1st <- mean(sample_path_rw2_1stDeriv_upper >= diff(compute_g(z)) & sample_path_rw2_1stDeriv_lower <= diff(compute_g(z)) )
  
  
  sample_path_arima_1stDeriv <- apply(sample_path_arima, 2, diff)
  sample_path_arima_1stDeriv_mean <- apply(sample_path_arima_1stDeriv, 1, mean)
  sample_path_arima_1stDeriv_upper <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.95)
  sample_path_arima_1stDeriv_lower <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.05)
  CR_arima_1st <- mean(sample_path_arima_1stDeriv_upper >= diff(compute_g(z)) & sample_path_arima_1stDeriv_lower <= diff(compute_g(z)) )
  
  
  result_1st <- tibble(rIAE_RW2_1st = sqrt(mean(abs(diff(compute_g(z))/h - sample_path_rw2_1stDeriv_mean/h))), rIAE_ARIMA_1st = sqrt(mean(abs(diff(compute_g(z))/h - sample_path_arima_1stDeriv_mean/h))), MCI_RW2_1st = mean(sample_path_rw2_1stDeriv_upper - sample_path_rw2_1stDeriv_lower)/h, MCI_ARIMA_1st = mean(sample_path_arima_1stDeriv_upper - sample_path_arima_1stDeriv_lower)/h, CR_RW2_1st = CR_rw2_1st, CR_ARIMA_1st = CR_arima_1st)
  result <- cbind(result, result_1st)
  
  ### Second deriv
  sample_path_rw2_2ndDeriv <- apply(sample_path_rw2, 2, diff, differences = 2)
  sample_path_rw2_2ndDeriv_mean <- apply(sample_path_rw2_2ndDeriv, 1, mean)
  sample_path_rw2_2ndDeriv_upper <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.95)
  sample_path_rw2_2ndDeriv_lower <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.05)
  CR_rw2_2nd <- mean(sample_path_rw2_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_rw2_2ndDeriv_lower <= diff(compute_g(z), differences = 2))
  
  sample_path_arima_2ndDeriv <- apply(sample_path_arima, 2, diff, differences = 2)
  sample_path_arima_2ndDeriv_mean <- apply(sample_path_arima_2ndDeriv, 1, mean)
  sample_path_arima_2ndDeriv_upper <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.95)
  sample_path_arima_2ndDeriv_lower <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.05)
  CR_arima_2nd <- mean(sample_path_arima_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_arima_2ndDeriv_lower <= diff(compute_g(z), differences = 2))
  
  result_2nd <- tibble(rIAE_RW2_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2)/(h^2) - sample_path_rw2_2ndDeriv_mean/(h^2)))), rIAE_ARIMA_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2)/(h^2) - sample_path_arima_2ndDeriv_mean/(h^2)))), MCI_RW2_2nd = mean(sample_path_rw2_2ndDeriv_upper - sample_path_rw2_2ndDeriv_lower)/(h^2), MCI_ARIMA_2nd = mean(sample_path_arima_2ndDeriv_upper - sample_path_arima_2ndDeriv_lower)/(h^2), CR_RW2_2nd = CR_rw2_2nd, CR_ARIMA_2nd = CR_arima_2nd)
  result <- cbind(result, result_2nd)
  result
  
}




set.seed(123)
time_begin <- Sys.time()
result <- replicate_for_summary_once(z = z, k_n = 10, x_indx = x_indx, Q1 = Q1, Q3 = Q3)
Sys.time() - time_begin
### take 3 secs for a single run




############ Updated Simulations:

#### very Sparse case:
dis <- 20
z <- round(seq(0.5,100, 0.5),2)
x <- seq(1, 100, dis*0.5)
x_indx <- which(z %in% x)
construct_X <- construct_A(z, x_indx)
d <- diff(z)
H <- compute_H_rue(d,n = length(z))
B <- compute_B(d,n = length(z))
A <- compute_A(d, n = length(z))

###### RW2:
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(length(z), x = 0.00001), "dgTMatrix")
###### ARIMA:
D <- H[-c(1,length(z)),]
R <- B[-c(1,length(z)), -c(1,length(z))]
Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(length(z), x = 0.00001)), "dgTMatrix")

set.seed(123)
time_begin <- Sys.time()
result_n100_repeat5  <- foreach(i = 1:1000, .combine='rbind', .packages = c("aghq", "TMB", "Matrix", "tidyverse")) %dopar% {
  replicate_for_summary_once(z = z, x_indx = x_indx, k_n = 10, n_samp = 2000, Q1 = Q1, Q3 = Q3)
}
time_end <- Sys.time()
time_end - time_begin
save(file = "resultCase1_updated.rda", result_n100_repeat5)

result_n100_repeat5 %>% apply(MARGIN = 2, mean)




#### sparse case: (will include)
dis <- 4
z <- round(seq(0.5,100, 0.5),2)
x <- seq(1, 100, dis*0.5)
x_indx <- which(z %in% x)
construct_X <- construct_A(z, x_indx)
d <- diff(z)
H <- compute_H_rue(d,n = length(z))
B <- compute_B(d,n = length(z))
A <- compute_A(d, n = length(z))

###### RW2:
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(length(z), x = 0.00001), "dgTMatrix")
###### ARIMA:
D <- H[-c(1,length(z)),]
R <- B[-c(1,length(z)), -c(1,length(z))]
Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(length(z), x = 0.00001)), "dgTMatrix")

set.seed(123)
time_begin <- Sys.time()
result_n100_repeated_2times  <- foreach(i = 1:1000, .combine='rbind', .packages = c("aghq", "TMB", "Matrix", "tidyverse")) %dopar% {
  replicate_for_summary_once(z = z, x_indx = x_indx, k_n = 2, n_samp = 2000, Q1 = Q1, Q3 = Q3)
}
time_end <- Sys.time()
time_end - time_begin
save(file = "resultCase2_updated.rda", result_n100_repeated_2times)

result_n100_repeated_2times %>% apply(MARGIN = 2, mean)






#### very dense case: (will include)
dis <- 2
z <- round(seq(0.5,100, 0.5),2)
x <- seq(1, 100, dis*0.5)
x_indx <- which(z %in% x)
construct_X <- construct_A(z, x_indx)
d <- diff(z)
H <- compute_H_rue(d,n = length(z))
B <- compute_B(d,n = length(z))
A <- compute_A(d, n = length(z))

###### RW2:
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(length(z), x = 0.00001), "dgTMatrix")
###### ARIMA:
D <- H[-c(1,length(z)),]
R <- B[-c(1,length(z)), -c(1,length(z))]
Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(length(z), x = 0.00001)), "dgTMatrix")

set.seed(123)
time_begin <- Sys.time()
result_n100_unique  <- foreach(i = 1:1000, .combine='rbind', .packages = c("aghq", "TMB", "Matrix", "tidyverse")) %dopar% {
  replicate_for_summary_once(z = z, x_indx = x_indx, k_n = 1, n_samp = 2000, Q1 = Q1, Q3 = Q3)
}
time_end <- Sys.time()
time_end - time_begin
save(file = "result_VeryDense.rda", result_n100_unique)

result_n100_unique %>% apply(MARGIN = 2, mean)







parallel::stopCluster(cl)






#### plotting:

### Case 1:
load(file = "resultCase1_updated.rda")
rIAE_data <- result_n100_repeat5 %>% select(c("rIAE_RW2", "rIAE_ARIMA")) %>% gather(type, rIAE, rIAE_RW2:rIAE_ARIMA, factor_key=TRUE)
rIAE_data %>% ggplot() + geom_boxplot(aes(y = rIAE, fill = type)) + theme_classic() + ggtitle("rIAE for the function")

MCI_data <- result_n100_repeat5 %>% select(c("MCI_RW2", "MCI_ARIMA")) %>% gather(type, MCI, MCI_RW2:MCI_ARIMA, factor_key=TRUE)
MCI_data %>% ggplot() + geom_boxplot(aes(y = MCI, fill = type)) + theme_classic() + ggtitle("MCI for the function")

CR_data <- result_n100_repeat5 %>% select(c("CR_RW2", "CR_ARIMA")) %>% gather(type, CR, CR_RW2:CR_ARIMA, factor_key=TRUE)
CR_data %>% ggplot() + geom_boxplot(aes(y = CR, fill = type)) + theme_classic() + ggtitle("CR for the function")


rIAE_data <- result_n100_repeat5 %>% select(c("rIAE_RW2_1st", "rIAE_ARIMA_1st")) %>% gather(type, rIAE, rIAE_RW2_1st:rIAE_ARIMA_1st, factor_key=TRUE)
rIAE_data %>% ggplot() + geom_boxplot(aes(y = rIAE, fill = type)) + theme_classic() + ggtitle("rIAE for 1st derivative")


MCI_data <- result_n100_repeat5 %>% select(c("MCI_RW2_1st", "MCI_ARIMA_1st")) %>% gather(type, MCI, MCI_RW2_1st:MCI_ARIMA_1st, factor_key=TRUE)
MCI_data %>% ggplot() + geom_boxplot(aes(y = MCI, fill = type)) + theme_classic() + ggtitle("MCI for 1st derivative")

CR_data <- result_n100_repeat5 %>% select(c("CR_RW2_1st", "CR_ARIMA_1st")) %>% gather(type, CR, CR_RW2_1st:CR_ARIMA_1st, factor_key=TRUE)
CR_data %>% ggplot() + geom_boxplot(aes(y = CR, fill = type)) + theme_classic() + ggtitle("CR for 1st derivative")


rIAE_data <- result_n100_repeat5 %>% select(c("rIAE_RW2_2nd", "rIAE_ARIMA_2nd")) %>% gather(type, rIAE, rIAE_RW2_2nd:rIAE_ARIMA_2nd, factor_key=TRUE)
rIAE_data %>% ggplot() + geom_boxplot(aes(y = rIAE, fill = type)) + theme_classic() + ggtitle("rIAE for 2nd derivative")

MCI_data <- result_n100_repeat5 %>% select(c("MCI_RW2_2nd", "MCI_ARIMA_2nd")) %>% gather(type, MCI, MCI_RW2_2nd:MCI_ARIMA_2nd, factor_key=TRUE)
MCI_data %>% ggplot() + geom_boxplot(aes(y = MCI, fill = type)) + theme_classic() + ggtitle("MCI for 2nd derivative")

CR_data <- result_n100_repeat5 %>% select(c("CR_RW2_2nd", "CR_ARIMA_2nd")) %>% gather(type, CR, CR_RW2_2nd:CR_ARIMA_2nd, factor_key=TRUE)
CR_data %>% ggplot() + geom_boxplot(aes(y = CR, fill = type)) + theme_classic() + ggtitle("CR for 2nd derivative")



### Case 2:
load(file = "resultCase1_updated.rda")
rIAE_data <- result_n100_repeated_2times %>% select(c("rIAE_RW2", "rIAE_ARIMA")) %>% gather(type, rIAE, rIAE_RW2:rIAE_ARIMA, factor_key=TRUE)
rIAE_data %>% ggplot() + geom_boxplot(aes(y = rIAE, fill = type)) + theme_classic() + ggtitle("rIAE for the function")

MCI_data <- result_n100_repeated_2times %>% select(c("MCI_RW2", "MCI_ARIMA")) %>% gather(type, MCI, MCI_RW2:MCI_ARIMA, factor_key=TRUE)
MCI_data %>% ggplot() + geom_boxplot(aes(y = MCI, fill = type)) + theme_classic() + ggtitle("MCI for the function")

CR_data <- result_n100_repeated_2times %>% select(c("CR_RW2", "CR_ARIMA")) %>% gather(type, CR, CR_RW2:CR_ARIMA, factor_key=TRUE)
CR_data %>% ggplot() + geom_boxplot(aes(y = CR, fill = type)) + theme_classic() + ggtitle("CR for the function")


rIAE_data <- result_n100_repeated_2times %>% select(c("rIAE_RW2_1st", "rIAE_ARIMA_1st")) %>% gather(type, rIAE, rIAE_RW2_1st:rIAE_ARIMA_1st, factor_key=TRUE)
rIAE_data %>% ggplot() + geom_boxplot(aes(y = rIAE, fill = type)) + theme_classic() + ggtitle("rIAE for 1st derivative")


MCI_data <- result_n100_repeated_2times %>% select(c("MCI_RW2_1st", "MCI_ARIMA_1st")) %>% gather(type, MCI, MCI_RW2_1st:MCI_ARIMA_1st, factor_key=TRUE)
MCI_data %>% ggplot() + geom_boxplot(aes(y = MCI, fill = type)) + theme_classic() + ggtitle("MCI for 1st derivative")

CR_data <- result_n100_repeated_2times %>% select(c("CR_RW2_1st", "CR_ARIMA_1st")) %>% gather(type, CR, CR_RW2_1st:CR_ARIMA_1st, factor_key=TRUE)
CR_data %>% ggplot() + geom_boxplot(aes(y = CR, fill = type)) + theme_classic() + ggtitle("CR for 1st derivative")


rIAE_data <- result_n100_repeated_2times %>% select(c("rIAE_RW2_2nd", "rIAE_ARIMA_2nd")) %>% gather(type, rIAE, rIAE_RW2_2nd:rIAE_ARIMA_2nd, factor_key=TRUE)
rIAE_data %>% ggplot() + geom_boxplot(aes(y = rIAE, fill = type)) + theme_classic() + ggtitle("rIAE for 2nd derivative")

MCI_data <- result_n100_repeated_2times %>% select(c("MCI_RW2_2nd", "MCI_ARIMA_2nd")) %>% gather(type, MCI, MCI_RW2_2nd:MCI_ARIMA_2nd, factor_key=TRUE)
MCI_data %>% ggplot() + geom_boxplot(aes(y = MCI, fill = type)) + theme_classic() + ggtitle("MCI for 2nd derivative")

CR_data <- result_n100_repeated_2times %>% select(c("CR_RW2_2nd", "CR_ARIMA_2nd")) %>% gather(type, CR, CR_RW2_2nd:CR_ARIMA_2nd, factor_key=TRUE)
CR_data %>% ggplot() + geom_boxplot(aes(y = CR, fill = type)) + theme_classic() + ggtitle("CR for 2nd derivative")







####
### Case 3:
load(file = "result_VeryDense.rda")
rIAE_data <- result_n100_unique %>% select(c("rIAE_RW2", "rIAE_ARIMA")) %>% gather(type, rIAE, rIAE_RW2:rIAE_ARIMA, factor_key=TRUE)
rIAE_data %>% ggplot() + geom_boxplot(aes(y = rIAE, fill = type)) + theme_classic() + ggtitle("rIAE for the function")

MCI_data <- result_n100_unique %>% select(c("MCI_RW2", "MCI_ARIMA")) %>% gather(type, MCI, MCI_RW2:MCI_ARIMA, factor_key=TRUE)
MCI_data %>% ggplot() + geom_boxplot(aes(y = MCI, fill = type)) + theme_classic() + ggtitle("MCI for the function")

CR_data <- result_n100_unique %>% select(c("CR_RW2", "CR_ARIMA")) %>% gather(type, CR, CR_RW2:CR_ARIMA, factor_key=TRUE)
CR_data %>% ggplot() + geom_boxplot(aes(y = CR, fill = type)) + theme_classic() + ggtitle("CR for the function")


rIAE_data <- result_n100_unique %>% select(c("rIAE_RW2_1st", "rIAE_ARIMA_1st")) %>% gather(type, rIAE, rIAE_RW2_1st:rIAE_ARIMA_1st, factor_key=TRUE)
rIAE_data %>% ggplot() + geom_boxplot(aes(y = rIAE, fill = type)) + theme_classic() + ggtitle("rIAE for 1st derivative")


MCI_data <- result_n100_unique %>% select(c("MCI_RW2_1st", "MCI_ARIMA_1st")) %>% gather(type, MCI, MCI_RW2_1st:MCI_ARIMA_1st, factor_key=TRUE)
MCI_data %>% ggplot() + geom_boxplot(aes(y = MCI, fill = type)) + theme_classic() + ggtitle("MCI for 1st derivative")

CR_data <- result_n100_unique %>% select(c("CR_RW2_1st", "CR_ARIMA_1st")) %>% gather(type, CR, CR_RW2_1st:CR_ARIMA_1st, factor_key=TRUE)
CR_data %>% ggplot() + geom_boxplot(aes(y = CR, fill = type)) + theme_classic() + ggtitle("CR for 1st derivative")


rIAE_data <- result_n100_unique %>% select(c("rIAE_RW2_2nd", "rIAE_ARIMA_2nd")) %>% gather(type, rIAE, rIAE_RW2_2nd:rIAE_ARIMA_2nd, factor_key=TRUE)
rIAE_data %>% ggplot() + geom_boxplot(aes(y = rIAE, fill = type)) + theme_classic() + ggtitle("rIAE for 2nd derivative")

MCI_data <- result_n100_unique %>% select(c("MCI_RW2_2nd", "MCI_ARIMA_2nd")) %>% gather(type, MCI, MCI_RW2_2nd:MCI_ARIMA_2nd, factor_key=TRUE)
MCI_data %>% ggplot() + geom_boxplot(aes(y = MCI, fill = type)) + theme_classic() + ggtitle("MCI for 2nd derivative")

CR_data <- result_n100_unique %>% select(c("CR_RW2_2nd", "CR_ARIMA_2nd")) %>% gather(type, CR, CR_RW2_2nd:CR_ARIMA_2nd, factor_key=TRUE)
CR_data %>% ggplot() + geom_boxplot(aes(y = CR, fill = type)) + theme_classic() + ggtitle("CR for 2nd derivative")






############################################################################






















































####### Plot for a single iteration: for each scenario above

#############################################################
#############################################################
### Case 1: sparse spacing
dis <- 20
z <- round(seq(0.5,100, 0.5),2)
x <- seq(1, 100, dis*0.5)
x_indx <- which(z %in% x)
construct_X <- construct_A(z, x_indx)
d <- diff(z)
H <- compute_H_rue(d,n = length(z))
B <- compute_B(d,n = length(z))
A <- compute_A(d, n = length(z))
k_n <- 10

###### RW2:
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(length(z), x = 1e-5), "dgTMatrix")
###### ARIMA:
D <- H[-c(1,length(z)),]
R <- B[-c(1,length(z)), -c(1,length(z))]
Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(length(z), x = 1e-5)), "dgTMatrix")


x <- z[x_indx]
compute_g <- function(x){
  5*sin(0.1*x)
}
k <- length(x)
n <- k * k_n
y <- compute_g(rep(x, each = k_n)) + rnorm(n, sd = sqrt(3))
d <- diff(z)
h <- mean(d)
m <- length(z[-x_indx])
X <- as(construct_X, "dgTMatrix")
X <- X[rep(c(1:k), each = k_n),]

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
  W = rep(0, length(z)), # W = c(U); U = B-Spline coefficients
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
gz <- sample_marginal(quad, n_samp)
sample_path_rw2 <- gz$samps
sample_path_rw2 <- as.tibble(as.matrix(sample_path_rw2))
mean_y_rw2 <- apply(sample_path_rw2, 1, mean)
upper_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.95)
lower_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.05)
rIAE_rw2 <- sqrt(mean(abs(compute_g(z) - mean_y_rw2)))
MCI_rw2 <- mean(upper_y_rw2 - lower_y_rw2)
CR_rw2 <- mean(upper_y_rw2 >= compute_g(z) & lower_y_rw2 <= compute_g(z))
sample_path_rw2_x <- construct_X %*% gz$samps
mean_x_rw2 <- sample_path_rw2_x %>% apply(MARGIN = 1, FUN = mean)

plot(mean_y_rw2 ~ z, type = 'l', xlab = "region of interest", ylab = "y", col = "red", ylim = c(-10,10), lty = 1)
lines(upper_y_rw2 ~ z, lty = 2, col = 'orange')
lines(lower_y_rw2 ~ z, lty = 2, col = 'orange')
lines(compute_g(z) ~ z, lty = 3, col = 'black')
points(mean_x_rw2 ~ x, col = 'red')
for (i in sample.int(n_samp,30)) {
  lines(unlist(sample_path_rw2[,i]) ~ z, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "RW2 for sparse covariate: g(.)")


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
gz <- sample_marginal(quad, n_samp)
sample_path_arima <- gz$samps



mean_y_arima <- apply(sample_path_arima, 1, mean)
upper_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.95)
lower_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.05)
rIAE_arima <- sqrt(mean(abs(compute_g(z) - mean_y_arima)))
MCI_arima <- mean(upper_y_arima - lower_y_arima)
CR_arima <- mean(upper_y_arima >= compute_g(z) & lower_y_arima <= compute_g(z))
result <- tibble(rIAE_RW2 = rIAE_rw2, rIAE_ARIMA = rIAE_arima, MCI_RW2 = MCI_rw2, MCI_ARIMA = MCI_arima, CR_RW2 = CR_rw2, CR_ARIMA = CR_arima)

sample_path_arima_x <- construct_X %*% gz$samps
mean_x_arima <- sample_path_arima_x %>% apply(MARGIN = 1, FUN = mean)

plot(mean_y_arima ~ z, type = 'l', xlab = "region of interest", ylab = "y", col = "red", ylim = c(-10,10), lty = 1)
lines(upper_y_arima ~ z, lty = 2, col = 'orange')
lines(lower_y_arima ~ z, lty = 2, col = 'orange')
lines(compute_g(z) ~ z, lty = 3, col = 'black')
points(mean_x_arima ~ x, col = 'red')
for (i in sample.int(n_samp,30)) {
  lines(unlist(sample_path_arima[,i]) ~ z, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "ARIMA for sparse covariate: g(.)")





#### First deriv:
sample_path_rw2_1stDeriv <- apply(sample_path_rw2, 2, diff)
sample_path_rw2_1stDeriv_mean <- apply(sample_path_rw2_1stDeriv, 1, mean)
sample_path_rw2_1stDeriv_upper <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.95)
sample_path_rw2_1stDeriv_lower <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.05)
CR_rw2_1st <- mean(sample_path_rw2_1stDeriv_upper >= diff(compute_g(z)) & sample_path_rw2_1stDeriv_lower <= diff(compute_g(z)) )


sample_path_arima_1stDeriv <- apply(sample_path_arima, 2, diff)
sample_path_arima_1stDeriv_mean <- apply(sample_path_arima_1stDeriv, 1, mean)
sample_path_arima_1stDeriv_upper <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.95)
sample_path_arima_1stDeriv_lower <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.05)
CR_arima_1st <- mean(sample_path_arima_1stDeriv_upper >= diff(compute_g(z)) & sample_path_arima_1stDeriv_lower <= diff(compute_g(z)) )


result_1st <- tibble(rIAE_RW2_1st = sqrt(mean(abs(diff(compute_g(z))/h - sample_path_rw2_1stDeriv_mean/h))), rIAE_ARIMA_1st = sqrt(mean(abs(diff(compute_g(z))/h - sample_path_arima_1stDeriv_mean/h))), MCI_RW2_1st = mean(sample_path_rw2_1stDeriv_upper - sample_path_rw2_1stDeriv_lower)/h, MCI_ARIMA_1st = mean(sample_path_arima_1stDeriv_upper - sample_path_arima_1stDeriv_lower)/h, CR_RW2_1st = CR_rw2_1st, CR_ARIMA_1st = CR_arima_1st)
result <- cbind(result, result_1st)

### Second deriv
sample_path_rw2_2ndDeriv <- apply(sample_path_rw2, 2, diff, differences = 2)
sample_path_rw2_2ndDeriv_mean <- apply(sample_path_rw2_2ndDeriv, 1, mean)
sample_path_rw2_2ndDeriv_upper <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.95)
sample_path_rw2_2ndDeriv_lower <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.05)
CR_rw2_2nd <- mean(sample_path_rw2_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_rw2_2ndDeriv_lower <= diff(compute_g(z), differences = 2))

sample_path_arima_2ndDeriv <- apply(sample_path_arima, 2, diff, differences = 2)
sample_path_arima_2ndDeriv_mean <- apply(sample_path_arima_2ndDeriv, 1, mean)
sample_path_arima_2ndDeriv_upper <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.95)
sample_path_arima_2ndDeriv_lower <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.05)
CR_arima_2nd <- mean(sample_path_arima_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_arima_2ndDeriv_lower <= diff(compute_g(z), differences = 2))

result_2nd <- tibble(rIAE_RW2_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2)/(h^2) - sample_path_rw2_2ndDeriv_mean/(h^2)))), rIAE_ARIMA_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2)/(h^2) - sample_path_arima_2ndDeriv_mean/(h^2)))), MCI_RW2_2nd = mean(sample_path_rw2_2ndDeriv_upper - sample_path_rw2_2ndDeriv_lower)/(h^2), MCI_ARIMA_2nd = mean(sample_path_arima_2ndDeriv_upper - sample_path_arima_2ndDeriv_lower)/(h^2), CR_RW2_2nd = CR_rw2_2nd, CR_ARIMA_2nd = CR_arima_2nd)




plot(diff(compute_g(z))/h ~ z[-1], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "1st deriv", ylim = c(-3,3), main = "RW2 for sparse covariate: 1st deriv")
lines(sample_path_rw2_1stDeriv_mean/h ~ z[-1], col = 'red')
lines(sample_path_rw2_1stDeriv_upper/h ~ z[-1], col = 'orange', lty = 2)
lines(sample_path_rw2_1stDeriv_lower/h ~ z[-1], col = 'orange', lty = 2)

for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_rw2_1stDeriv[,i])/h ~ z[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}






plot(diff(compute_g(z))/h ~ z[-1], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "1st deriv", ylim = c(-3,3), main = "ARIMA for sparse covariate: 1st deriv")
lines(sample_path_arima_1stDeriv_mean/h ~ z[-1], col = 'red')
lines(sample_path_arima_1stDeriv_upper/h ~ z[-1], col = 'orange', lty = 2)
lines(sample_path_arima_1stDeriv_lower/h ~ z[-1], col = 'orange', lty = 2)

for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_arima_1stDeriv[,i])/h ~ z[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}





plot((diff(compute_g(z), differences = 2)/(h^2)) ~ z[-c(1,2)], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "2nd deriv", ylim = c(-1,1), main = "RW2 for sparse covariate: 2nd deriv")
lines(sample_path_rw2_2ndDeriv_mean/(h^2) ~ z[-c(1,2)], col = 'red')
lines(sample_path_rw2_2ndDeriv_upper/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)
lines(sample_path_rw2_2ndDeriv_lower/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)


for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_rw2_2ndDeriv[,i])/(h^2) ~ z[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}



plot(diff(compute_g(z), differences = 2)/(h^2) ~ z[-c(1,2)], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "2nd deriv", ylim = c(-1,1), main = "ARIMA for sparse covariate: 2nd deriv")
lines(sample_path_arima_2ndDeriv_mean/(h^2) ~ z[-c(1,2)], col = 'red')
lines(sample_path_arima_2ndDeriv_upper/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)
lines(sample_path_arima_2ndDeriv_lower/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)

for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_arima_2ndDeriv[,i])/(h^2) ~ z[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}







#############################################################
#############################################################
### Case 2: medium spacing
dis <- 4
z <- round(seq(0.5,100, 0.5),2)
x <- seq(1, 100, dis*0.5)
x_indx <- which(z %in% x)
construct_X <- construct_A(z, x_indx)
d <- diff(z)
H <- compute_H_rue(d,n = length(z))
B <- compute_B(d,n = length(z))
A <- compute_A(d, n = length(z))
k_n <- 2

###### RW2:
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(length(z), x = 1e-5), "dgTMatrix")
###### ARIMA:
D <- H[-c(1,length(z)),]
R <- B[-c(1,length(z)), -c(1,length(z))]
Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(length(z), x = 1e-5)), "dgTMatrix")


x <- z[x_indx]
compute_g <- function(x){
  5*sin(0.1*x)
}
k <- length(x)
n <- k * k_n
y <- compute_g(rep(x, each = k_n)) + rnorm(n, sd = sqrt(3))
d <- diff(z)
h <- mean(d)
m <- length(z[-x_indx])
X <- as(construct_X, "dgTMatrix")
X <- X[rep(c(1:k), each = k_n),]

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
  W = rep(0, length(z)), # W = c(U); U = B-Spline coefficients
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
gz <- sample_marginal(quad, n_samp)
sample_path_rw2 <- gz$samps
sample_path_rw2 <- as.tibble(as.matrix(sample_path_rw2))
mean_y_rw2 <- apply(sample_path_rw2, 1, mean)
upper_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.95)
lower_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.05)
rIAE_rw2 <- sqrt(mean(abs(compute_g(z) - mean_y_rw2)))
MCI_rw2 <- mean(upper_y_rw2 - lower_y_rw2)
CR_rw2 <- mean(upper_y_rw2 >= compute_g(z) & lower_y_rw2 <= compute_g(z))
sample_path_rw2_x <- construct_X %*% gz$samps
mean_x_rw2 <- sample_path_rw2_x %>% apply(MARGIN = 1, FUN = mean)

plot(mean_y_rw2 ~ z, type = 'l', xlab = "region of interest", ylab = "y", col = "red", ylim = c(-10,10), lty = 1)
lines(upper_y_rw2 ~ z, lty = 2, col = 'orange')
lines(lower_y_rw2 ~ z, lty = 2, col = 'orange')
lines(compute_g(z) ~ z, lty = 3, col = 'black')
points(mean_x_rw2 ~ x, col = 'red')
for (i in sample.int(n_samp,30)) {
  lines(unlist(sample_path_rw2[,i]) ~ z, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "RW2 for regular covariate: g(.)")


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
gz <- sample_marginal(quad, n_samp)
sample_path_arima <- gz$samps



mean_y_arima <- apply(sample_path_arima, 1, mean)
upper_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.95)
lower_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.05)
rIAE_arima <- sqrt(mean(abs(compute_g(z) - mean_y_arima)))
MCI_arima <- mean(upper_y_arima - lower_y_arima)
CR_arima <- mean(upper_y_arima >= compute_g(z) & lower_y_arima <= compute_g(z))
result <- tibble(rIAE_RW2 = rIAE_rw2, rIAE_ARIMA = rIAE_arima, MCI_RW2 = MCI_rw2, MCI_ARIMA = MCI_arima, CR_RW2 = CR_rw2, CR_ARIMA = CR_arima)

sample_path_arima_x <- construct_X %*% gz$samps
mean_x_arima <- sample_path_arima_x %>% apply(MARGIN = 1, FUN = mean)

plot(mean_y_arima ~ z, type = 'l', xlab = "region of interest", ylab = "y", col = "red", ylim = c(-10,10), lty = 1)
lines(upper_y_arima ~ z, lty = 2, col = 'orange')
lines(lower_y_arima ~ z, lty = 2, col = 'orange')
lines(compute_g(z) ~ z, lty = 3, col = 'black')
points(mean_x_arima ~ x, col = 'red')
for (i in sample.int(n_samp,30)) {
  lines(unlist(sample_path_arima[,i]) ~ z, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "ARIMA for regular covariate: g(.)")





#### First deriv:
sample_path_rw2_1stDeriv <- apply(sample_path_rw2, 2, diff)
sample_path_rw2_1stDeriv_mean <- apply(sample_path_rw2_1stDeriv, 1, mean)
sample_path_rw2_1stDeriv_upper <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.95)
sample_path_rw2_1stDeriv_lower <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.05)
CR_rw2_1st <- mean(sample_path_rw2_1stDeriv_upper >= diff(compute_g(z)) & sample_path_rw2_1stDeriv_lower <= diff(compute_g(z)) )


sample_path_arima_1stDeriv <- apply(sample_path_arima, 2, diff)
sample_path_arima_1stDeriv_mean <- apply(sample_path_arima_1stDeriv, 1, mean)
sample_path_arima_1stDeriv_upper <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.95)
sample_path_arima_1stDeriv_lower <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.05)
CR_arima_1st <- mean(sample_path_arima_1stDeriv_upper >= diff(compute_g(z)) & sample_path_arima_1stDeriv_lower <= diff(compute_g(z)) )


result_1st <- tibble(rIAE_RW2_1st = sqrt(mean(abs(diff(compute_g(z))/h - sample_path_rw2_1stDeriv_mean/h))), rIAE_ARIMA_1st = sqrt(mean(abs(diff(compute_g(z))/h - sample_path_arima_1stDeriv_mean/h))), MCI_RW2_1st = mean(sample_path_rw2_1stDeriv_upper - sample_path_rw2_1stDeriv_lower)/h, MCI_ARIMA_1st = mean(sample_path_arima_1stDeriv_upper - sample_path_arima_1stDeriv_lower)/h, CR_RW2_1st = CR_rw2_1st, CR_ARIMA_1st = CR_arima_1st)
result <- cbind(result, result_1st)

### Second deriv
sample_path_rw2_2ndDeriv <- apply(sample_path_rw2, 2, diff, differences = 2)
sample_path_rw2_2ndDeriv_mean <- apply(sample_path_rw2_2ndDeriv, 1, mean)
sample_path_rw2_2ndDeriv_upper <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.95)
sample_path_rw2_2ndDeriv_lower <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.05)
CR_rw2_2nd <- mean(sample_path_rw2_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_rw2_2ndDeriv_lower <= diff(compute_g(z), differences = 2))

sample_path_arima_2ndDeriv <- apply(sample_path_arima, 2, diff, differences = 2)
sample_path_arima_2ndDeriv_mean <- apply(sample_path_arima_2ndDeriv, 1, mean)
sample_path_arima_2ndDeriv_upper <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.95)
sample_path_arima_2ndDeriv_lower <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.05)
CR_arima_2nd <- mean(sample_path_arima_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_arima_2ndDeriv_lower <= diff(compute_g(z), differences = 2))

result_2nd <- tibble(rIAE_RW2_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2)/(h^2) - sample_path_rw2_2ndDeriv_mean/(h^2)))), rIAE_ARIMA_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2)/(h^2) - sample_path_arima_2ndDeriv_mean/(h^2)))), MCI_RW2_2nd = mean(sample_path_rw2_2ndDeriv_upper - sample_path_rw2_2ndDeriv_lower)/(h^2), MCI_ARIMA_2nd = mean(sample_path_arima_2ndDeriv_upper - sample_path_arima_2ndDeriv_lower)/(h^2), CR_RW2_2nd = CR_rw2_2nd, CR_ARIMA_2nd = CR_arima_2nd)




plot(diff(compute_g(z))/h ~ z[-1], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "1st deriv", ylim = c(-3,3), main = "RW2 for regular covariate: 1st deriv")
lines(sample_path_rw2_1stDeriv_mean/h ~ z[-1], col = 'red')
lines(sample_path_rw2_1stDeriv_upper/h ~ z[-1], col = 'orange', lty = 2)
lines(sample_path_rw2_1stDeriv_lower/h ~ z[-1], col = 'orange', lty = 2)

for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_rw2_1stDeriv[,i])/h ~ z[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}






plot(diff(compute_g(z))/h ~ z[-1], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "1st deriv", ylim = c(-3,3), main = "ARIMA for regular covariate: 1st deriv")
lines(sample_path_arima_1stDeriv_mean/h ~ z[-1], col = 'red')
lines(sample_path_arima_1stDeriv_upper/h ~ z[-1], col = 'orange', lty = 2)
lines(sample_path_arima_1stDeriv_lower/h ~ z[-1], col = 'orange', lty = 2)

for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_arima_1stDeriv[,i])/h ~ z[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}





plot((diff(compute_g(z), differences = 2)/(h^2)) ~ z[-c(1,2)], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "2nd deriv", ylim = c(-1,1), main = "RW2 for regular covariate: 2nd deriv")
lines(sample_path_rw2_2ndDeriv_mean/(h^2) ~ z[-c(1,2)], col = 'red')
lines(sample_path_rw2_2ndDeriv_upper/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)
lines(sample_path_rw2_2ndDeriv_lower/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)


for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_rw2_2ndDeriv[,i])/(h^2) ~ z[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}



plot(diff(compute_g(z), differences = 2)/(h^2) ~ z[-c(1,2)], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "2nd deriv", ylim = c(-1,1), main = "ARIMA for regular covariate: 2nd deriv")
lines(sample_path_arima_2ndDeriv_mean/(h^2) ~ z[-c(1,2)], col = 'red')
lines(sample_path_arima_2ndDeriv_upper/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)
lines(sample_path_arima_2ndDeriv_lower/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)

for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_arima_2ndDeriv[,i])/(h^2) ~ z[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}


















#############################################################
#############################################################
### Case 3: dense spacing
dis <- 2
z <- round(seq(0.5,100, 0.5),2)
x <- seq(1, 100, dis*0.5)
x_indx <- which(z %in% x)
construct_X <- construct_A(z, x_indx)
d <- diff(z)
H <- compute_H_rue(d,n = length(z))
B <- compute_B(d,n = length(z))
A <- compute_A(d, n = length(z))
k_n <- 1

###### RW2:
Q1 <- t(H) %*% solve(A) %*% H
Q1 <- as(Q1 + Diagonal(length(z), x = 1e-5), "dgTMatrix")
###### ARIMA:
D <- H[-c(1,length(z)),]
R <- B[-c(1,length(z)), -c(1,length(z))]
Q3 <- t(D) %*% solve(R) %*% D
Q3 <- as(as.matrix(Q3 + Diagonal(length(z), x = 1e-5)), "dgTMatrix")


x <- z[x_indx]
compute_g <- function(x){
  5*sin(0.1*x)
}
k <- length(x)
n <- k * k_n
y <- compute_g(rep(x, each = k_n)) + rnorm(n, sd = sqrt(3))
d <- diff(z)
h <- mean(d)
m <- length(z[-x_indx])
X <- as(construct_X, "dgTMatrix")
X <- X[rep(c(1:k), each = k_n),]

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
  W = rep(0, length(z)), # W = c(U); U = B-Spline coefficients
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
gz <- sample_marginal(quad, n_samp)
sample_path_rw2 <- gz$samps
sample_path_rw2 <- as.tibble(as.matrix(sample_path_rw2))
mean_y_rw2 <- apply(sample_path_rw2, 1, mean)
upper_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.95)
lower_y_rw2 <- apply(sample_path_rw2, 1, quantile, p = 0.05)
rIAE_rw2 <- sqrt(mean(abs(compute_g(z) - mean_y_rw2)))
MCI_rw2 <- mean(upper_y_rw2 - lower_y_rw2)
CR_rw2 <- mean(upper_y_rw2 >= compute_g(z) & lower_y_rw2 <= compute_g(z))
sample_path_rw2_x <- construct_X %*% gz$samps
mean_x_rw2 <- sample_path_rw2_x %>% apply(MARGIN = 1, FUN = mean)

plot(mean_y_rw2 ~ z, type = 'l', xlab = "region of interest", ylab = "y", col = "red", ylim = c(-10,10), lty = 1)
lines(upper_y_rw2 ~ z, lty = 2, col = 'orange')
lines(lower_y_rw2 ~ z, lty = 2, col = 'orange')
lines(compute_g(z) ~ z, lty = 3, col = 'black')
points(mean_x_rw2 ~ x, col = 'red')
for (i in sample.int(n_samp,30)) {
  lines(unlist(sample_path_rw2[,i]) ~ z, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "RW2 for dense covariate: g(.)")


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
gz <- sample_marginal(quad, n_samp)
sample_path_arima <- gz$samps



mean_y_arima <- apply(sample_path_arima, 1, mean)
upper_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.95)
lower_y_arima <- apply(sample_path_arima, 1, quantile, p = 0.05)
rIAE_arima <- sqrt(mean(abs(compute_g(z) - mean_y_arima)))
MCI_arima <- mean(upper_y_arima - lower_y_arima)
CR_arima <- mean(upper_y_arima >= compute_g(z) & lower_y_arima <= compute_g(z))
result <- tibble(rIAE_RW2 = rIAE_rw2, rIAE_ARIMA = rIAE_arima, MCI_RW2 = MCI_rw2, MCI_ARIMA = MCI_arima, CR_RW2 = CR_rw2, CR_ARIMA = CR_arima)

sample_path_arima_x <- construct_X %*% gz$samps
mean_x_arima <- sample_path_arima_x %>% apply(MARGIN = 1, FUN = mean)

plot(mean_y_arima ~ z, type = 'l', xlab = "region of interest", ylab = "y", col = "red", ylim = c(-10,10), lty = 1)
lines(upper_y_arima ~ z, lty = 2, col = 'orange')
lines(lower_y_arima ~ z, lty = 2, col = 'orange')
lines(compute_g(z) ~ z, lty = 3, col = 'black')
points(mean_x_arima ~ x, col = 'red')
for (i in sample.int(n_samp,30)) {
  lines(unlist(sample_path_arima[,i]) ~ z, col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}
title(main = "ARIMA for dense covariate: g(.)")





#### First deriv:
sample_path_rw2_1stDeriv <- apply(sample_path_rw2, 2, diff)
sample_path_rw2_1stDeriv_mean <- apply(sample_path_rw2_1stDeriv, 1, mean)
sample_path_rw2_1stDeriv_upper <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.95)
sample_path_rw2_1stDeriv_lower <- apply(sample_path_rw2_1stDeriv, 1, quantile, p = 0.05)
CR_rw2_1st <- mean(sample_path_rw2_1stDeriv_upper >= diff(compute_g(z)) & sample_path_rw2_1stDeriv_lower <= diff(compute_g(z)) )


sample_path_arima_1stDeriv <- apply(sample_path_arima, 2, diff)
sample_path_arima_1stDeriv_mean <- apply(sample_path_arima_1stDeriv, 1, mean)
sample_path_arima_1stDeriv_upper <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.95)
sample_path_arima_1stDeriv_lower <- apply(sample_path_arima_1stDeriv, 1, quantile, p = 0.05)
CR_arima_1st <- mean(sample_path_arima_1stDeriv_upper >= diff(compute_g(z)) & sample_path_arima_1stDeriv_lower <= diff(compute_g(z)) )


result_1st <- tibble(rIAE_RW2_1st = sqrt(mean(abs(diff(compute_g(z))/h - sample_path_rw2_1stDeriv_mean/h))), rIAE_ARIMA_1st = sqrt(mean(abs(diff(compute_g(z))/h - sample_path_arima_1stDeriv_mean/h))), MCI_RW2_1st = mean(sample_path_rw2_1stDeriv_upper - sample_path_rw2_1stDeriv_lower)/h, MCI_ARIMA_1st = mean(sample_path_arima_1stDeriv_upper - sample_path_arima_1stDeriv_lower)/h, CR_RW2_1st = CR_rw2_1st, CR_ARIMA_1st = CR_arima_1st)
result <- cbind(result, result_1st)

### Second deriv
sample_path_rw2_2ndDeriv <- apply(sample_path_rw2, 2, diff, differences = 2)
sample_path_rw2_2ndDeriv_mean <- apply(sample_path_rw2_2ndDeriv, 1, mean)
sample_path_rw2_2ndDeriv_upper <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.95)
sample_path_rw2_2ndDeriv_lower <- apply(sample_path_rw2_2ndDeriv, 1, quantile, p = 0.05)
CR_rw2_2nd <- mean(sample_path_rw2_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_rw2_2ndDeriv_lower <= diff(compute_g(z), differences = 2))

sample_path_arima_2ndDeriv <- apply(sample_path_arima, 2, diff, differences = 2)
sample_path_arima_2ndDeriv_mean <- apply(sample_path_arima_2ndDeriv, 1, mean)
sample_path_arima_2ndDeriv_upper <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.95)
sample_path_arima_2ndDeriv_lower <- apply(sample_path_arima_2ndDeriv, 1, quantile, p = 0.05)
CR_arima_2nd <- mean(sample_path_arima_2ndDeriv_upper >= diff(compute_g(z), differences = 2) & sample_path_arima_2ndDeriv_lower <= diff(compute_g(z), differences = 2))

result_2nd <- tibble(rIAE_RW2_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2)/(h^2) - sample_path_rw2_2ndDeriv_mean/(h^2)))), rIAE_ARIMA_2nd = sqrt(mean(abs(diff(compute_g(z), differences = 2)/(h^2) - sample_path_arima_2ndDeriv_mean/(h^2)))), MCI_RW2_2nd = mean(sample_path_rw2_2ndDeriv_upper - sample_path_rw2_2ndDeriv_lower)/(h^2), MCI_ARIMA_2nd = mean(sample_path_arima_2ndDeriv_upper - sample_path_arima_2ndDeriv_lower)/(h^2), CR_RW2_2nd = CR_rw2_2nd, CR_ARIMA_2nd = CR_arima_2nd)




plot(diff(compute_g(z))/h ~ z[-1], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "1st deriv", ylim = c(-3,3), main = "RW2 for dense covariate: 1st deriv")
lines(sample_path_rw2_1stDeriv_mean/h ~ z[-1], col = 'red')
lines(sample_path_rw2_1stDeriv_upper/h ~ z[-1], col = 'orange', lty = 2)
lines(sample_path_rw2_1stDeriv_lower/h ~ z[-1], col = 'orange', lty = 2)

for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_rw2_1stDeriv[,i])/h ~ z[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}






plot(diff(compute_g(z))/h ~ z[-1], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "1st deriv", ylim = c(-3,3), main = "ARIMA for dense covariate: 1st deriv")
lines(sample_path_arima_1stDeriv_mean/h ~ z[-1], col = 'red')
lines(sample_path_arima_1stDeriv_upper/h ~ z[-1], col = 'orange', lty = 2)
lines(sample_path_arima_1stDeriv_lower/h ~ z[-1], col = 'orange', lty = 2)

for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_arima_1stDeriv[,i])/h ~ z[-1], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}





plot((diff(compute_g(z), differences = 2)/(h^2)) ~ z[-c(1,2)], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "2nd deriv", ylim = c(-1,1), main = "RW2 for dense covariate: 2nd deriv")
lines(sample_path_rw2_2ndDeriv_mean/(h^2) ~ z[-c(1,2)], col = 'red')
lines(sample_path_rw2_2ndDeriv_upper/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)
lines(sample_path_rw2_2ndDeriv_lower/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)


for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_rw2_2ndDeriv[,i])/(h^2) ~ z[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}



plot(diff(compute_g(z), differences = 2)/(h^2) ~ z[-c(1,2)], type = 'l', col = 'black', 
     xlab = "region of interest", 
     ylab = "2nd deriv", ylim = c(-1,1), main = "ARIMA for dense covariate: 2nd deriv")
lines(sample_path_arima_2ndDeriv_mean/(h^2) ~ z[-c(1,2)], col = 'red')
lines(sample_path_arima_2ndDeriv_upper/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)
lines(sample_path_arima_2ndDeriv_lower/(h^2) ~ z[-c(1,2)], col = 'orange', lty = 2)

for (i in sample.int(n_samp,5)) {
  lines(unlist(sample_path_arima_2ndDeriv[,i])/(h^2) ~ z[-c(1,2)], col = rgb(0, 0, 255, max = 255, alpha = 20, names = "grey"))
}








