##### All the relevant functions:

##################### For precision matrix construction:
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




################ For model-based interpolation:
switch_matrix_once <- function(M, from, to){
  L <- diag(1, nrow = nrow(M), ncol = nrow(M))
  original_row <- L[to,]
  L[to,] <- L[from,]
  L[from,] <- original_row
  as(L,"dgCMatrix")%*% as(M, "dgCMatrix") %*% as(t(L), "dgCMatrix")
}
switch_matrix_once_new <- function(M, from, to){
  v <- 1:nrow(M)
  v[from] <- to
  v[to] <- from
  as(M[v,v], "dgCMatrix")
}
switch_matrix <- function(z, x, M){
  all_grid <- sort(c(z,x))
  n <- length(x)
  for (i in 1:length(z)) {
    z_location <- which(all_grid %in% z[i])
    M <- switch_matrix_once_new(M,from = z_location, to = n + i)
    all_grid[z_location] <- all_grid[n+i]
    all_grid[n + i] <- z[i]
  }
  for (j in 1:length(x)) {
    x_location <- which(all_grid %in% x[j])
    M <- switch_matrix_once_new(M,from = x_location, to = j)
    all_grid[x_location] <- all_grid[j]
    all_grid[j] <- x[j]
  }
  M
}
Interpolation_vec_v1 <- function(t, x_grid, gx, GP){
  all_grid <- sort(c(t,x_grid))
  t_location <- which(all_grid %in% t)
  n <- length(all_grid)
  all_d <- diff(all_grid)
  H <- compute_H_rue(all_d, n)
  B <- compute_B(all_d, n)
  A <- compute_A(all_d, n)
  if(GP == "RW2"){
    Q <- as(t(H),"dgCMatrix") %*% as(solve(as(A,"dgCMatrix")),"dgCMatrix") %*% as(H,"dgCMatrix")
  }
  else{
    D <- H[-c(1,length(all_grid)),]
    R <- B[-c(1,length(all_grid)), -c(1,length(all_grid))]
    Q <- t(D) %*% solve(R) %*% D
    Q <- as(t(D),"dgCMatrix") %*% as(solve(as(R,"dgCMatrix")),"dgCMatrix") %*% as(D,"dgCMatrix")
  }
  C <- solve(as(Q + Diagonal(n, x = .Machine$double.eps),"dgCMatrix"))
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
## v2 uses diagonal pertubation to make the matrix invertible
Interpolation_vec_v2 <- function(t, x_grid, quad_samps, GP, mycores = 4){
  gx <- quad_samps$samps[1:ncol(Q1),]
  thetas <- quad_samps$theta$theta1
  sigmas <- exp(-0.5 * thetas)
  all_grid <- sort(c(t,x_grid))
  t_location <- which(all_grid %in% t)
  n <- length(all_grid)
  all_d <- diff(all_grid)
  H <- compute_H_rue(all_d, n)
  B <- compute_B(all_d, n)
  A <- compute_A(all_d, n)
  if(GP == "RW2"){
    Q <- as(t(H),"dgCMatrix") %*% as(solve(as(A,"dgCMatrix")),"dgCMatrix") %*% as(H,"dgCMatrix")
  }
  else{
    D <- H[-c(1,length(all_grid)),]
    R <- B[-c(1,length(all_grid)), -c(1,length(all_grid))]
    Q <- t(D) %*% solve(R) %*% D
    Q <- as(t(D),"dgCMatrix") %*% as(solve(as(R,"dgCMatrix")),"dgCMatrix") %*% as(D,"dgCMatrix")
  }
  C <- solve(as(Q + Diagonal(n, x = 1e-13),"dgCMatrix"))
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
  simulated_gz <- function(i){
    Lt <- chol(((1/sigmas[i])^2) * QAA)
    z <- rnorm(n = length(t))
    v <- solve(a = Lt, b = z)
    v + conditional_mean_list[[i]]
  }
  samples <- mclapply(1:length(conditional_mean_list), simulated_gz, mc.cores = mycores)
  samples
}
Interpolation_vec_v2_updated <- function(t, x_grid, quad_samps, GP, mycores = 4){
  gx <- quad_samps$samps[1:ncol(Q1),]
  thetas <- quad_samps$theta$theta1
  sigmas <- exp(-0.5 * thetas)
  all_grid <- sort(c(t,x_grid))
  t_location <- which(all_grid %in% t)
  n <- length(all_grid)
  all_d <- diff(all_grid)
  H <- compute_H_rue(all_d, n)
  B <- compute_B(all_d, n)
  A <- compute_A(all_d, n)
  if(GP == "RW2"){
    Q <- as(t(H),"dgCMatrix") %*% as(solve(as(A,"dgCMatrix")),"dgCMatrix") %*% as(H,"dgCMatrix")
  }
  else{
    D <- H[-c(1,length(all_grid)),]
    R <- B[-c(1,length(all_grid)), -c(1,length(all_grid))]
    Q <- t(D) %*% solve(R) %*% D
    Q <- as(t(D),"dgCMatrix") %*% as(solve(as(R,"dgCMatrix")),"dgCMatrix") %*% as(D,"dgCMatrix")
  }
  start_time <- Sys.time()
  C <- solve(as(Q + Diagonal(n, x = 1e-13),"dgCMatrix"))
  print(Sys.time() - start_time)
  print("time used to solve Q above")
  start_time <- Sys.time()
  C_trans <- switch_matrix(t,x_grid,C)
  print(Sys.time() - start_time)
  print("time used to switch matrix")
  start_time <- Sys.time()
  Q_trans <- Matrix::forceSymmetric(solve(C_trans))
  print(Sys.time() - start_time)
  print("time used to solve C matrix")
  QAA <- Q_trans[(length(x_grid) + 1):n, (length(x_grid) + 1):n]
  CAA <- solve(QAA)
  QAB <- Q_trans[(length(x_grid) + 1):n, 1:(length(x_grid))]
  updating_mean_matrix <- -CAA %*% QAB
  Lt <- chol(QAA)
  compute_conditional_mean <- function(gx){
    conditional_mean <- updating_mean_matrix %*% (gx)
    conditional_mean
  }
  conditional_mean_list <- apply(gx, 2, compute_conditional_mean)
  simulated_gz <- function(i){
    Lt_i <- ((1/sigmas[i]))*Lt
    z <- rnorm(n = length(t))
    v <- solve(a = Lt_i, b = z)
    v + conditional_mean_list[[i]]
  }
  samples <- mclapply(1:length(conditional_mean_list), simulated_gz, mc.cores = mycores)
  samples
}

## v3 uses eigen-transform to make the matrix invertible: (this method also seems not work)
Interpolation_vec_v3 <- function(t, x_grid, quad_samps, GP, mycores = 4){
  gx <- quad_samps$samps[1:ncol(Q1),]
  thetas <- quad_samps$theta$theta1
  sigmas <- exp(-0.5 * thetas)
  all_grid <- sort(c(t,x_grid))
  t_location <- which(all_grid %in% t)
  n <- length(all_grid)
  all_d <- diff(all_grid)
  H <- compute_H_rue(all_d, n)
  B <- compute_B(all_d, n)
  A <- compute_A(all_d, n)
  if(GP == "RW2"){
    Q <- as(t(H),"dgCMatrix") %*% as(solve(as(A,"dgCMatrix")),"dgCMatrix") %*% as(H,"dgCMatrix")
  }
  else{
    D <- H[-c(1,length(all_grid)),]
    R <- B[-c(1,length(all_grid)), -c(1,length(all_grid))]
    Q <- t(D) %*% solve(R) %*% D
    Q <- as(t(D),"dgCMatrix") %*% as(solve(as(R,"dgCMatrix")),"dgCMatrix") %*% as(D,"dgCMatrix")
  }
  C <- solve(eigenspace_const(Q,all_grid = all_grid, buffer = 1000))
  C_trans <- switch_matrix(t,x_grid,C)
  Q_trans <- Matrix::forceSymmetric(solve(C_trans))
  print(eigen(Q_trans, only.values = T)$values)
  QAA <- Q_trans[(length(x_grid) + 1):n, (length(x_grid) + 1):n]
  QAB <- Q_trans[(length(x_grid) + 1):n, 1:(length(x_grid))]
  QBB <- Q[1:(length(x_grid)), 1:(length(x_grid))]
  compute_conditional_mean <- function(gx){
    conditional_mean <- -solve(QAA)%*% QAB %*% (gx)
    conditional_mean
  }
  conditional_mean_list <- apply(gx, 2, compute_conditional_mean)
  simulated_gz <- function(i){
    Lt <- chol(((1/sigmas[i])^2) * QAA)
    z <- rnorm(n = length(t))
    v <- solve(a = Lt, b = z)
    v + conditional_mean_list[[i]]
  }
  # samples <- mclapply(1:length(conditional_mean_list), simulated_gz, mc.cores = mycores)
  samples <- lapply(1:length(conditional_mean_list), simulated_gz)
  samples
}


## Conditional Kriging:

construct_A <- function(all_grid, x_indx){
  A <- matrix(0, nrow = length(x_indx), ncol = length(all_grid))
  for (i in 1:nrow(A)) {
    A[i,x_indx[i]] <- 1
  }
  A
}

ck_method_1 <- function(t, x_grid, quad_samps, GP, mycores = 4){
  gx <- quad_samps$samps[1:ncol(Q1),]
  thetas <- quad_samps$theta$theta1
  sigmas <- exp(-0.5 * thetas)
  all_grid <- sort(c(t,x_grid))
  t_location <- which(all_grid %in% t)
  n <- length(all_grid)
  all_d <- diff(all_grid)
  H <- compute_H_rue(all_d, n)
  B <- compute_B(all_d, n)
  A <- compute_A(all_d, n)
  if(GP == "RW2"){
    Q <- as(t(H),"dgCMatrix") %*% as(solve(as(A,"dgCMatrix")),"dgCMatrix") %*% as(H,"dgCMatrix")
  }
  else{
    D <- H[-c(1,length(all_grid)),]
    R <- B[-c(1,length(all_grid)), -c(1,length(all_grid))]
    Q <- t(D) %*% solve(R) %*% D
    Q <- as(t(D),"dgCMatrix") %*% as(solve(as(R,"dgCMatrix")),"dgCMatrix") %*% as(D,"dgCMatrix")
  }
  ### Prep Step:
  Q_eigen <- eigen(Q, symmetric = T)
  r <- length(Q_eigen$values) - 2
  Q <- Q + 0.0001 * tcrossprod(Q_eigen$vectors[,(r+1)]) + 0.0001 * tcrossprod(Q_eigen$vectors[,(r+2)])
  C <- solve(as(Q,"dgCMatrix"))
  x_location <- which(all_grid %in% x_grid)
  begin_timeA <- Sys.time()
  A <- as(rbind(construct_A(all_grid, x_location), rep(1, nrow(Q))/sqrt(sum(rep(1, nrow(Q1)))) , all_grid/sqrt(sum(all_grid^2))), "dgCMatrix")
  end_timeA <- Sys.time()
  print(end_timeA - begin_timeA)
  print("Time to construct A")
  begin_time <- Sys.time()
  tCA <- tcrossprod(C,A)
  inv_ACA <- solve(A %*% tCA)
  end_time <- Sys.time()
  print(end_time - begin_time)
  print("Time to compute ACA inverse")
  ### Step 1: sample from unconditional distribution, project to intercept and slope
  begin_time <- Sys.time()
  x_tilde <- matrix(0, nrow = length(all_grid), ncol = ncol(gx))
  int <- matrix(0, nrow = 1, ncol = ncol(gx))
  slope <- matrix(0, nrow = 1, ncol = ncol(gx))
  sd_y <- 1/sqrt(Q_eigen$values[1:r])
  for (i in 1:ncol(gx)) {
    y <- rnorm(n = r, sd = sd_y*sigmas[i])
    x_tilde[,i] <- Q_eigen$vectors[,1:r] %*% y
  }
  for (j in 1:ncol(gx)) {
    int[,j] <- (crossprod(rep(1, nrow(Q1)), gx[,j]))/(sum(rep(1, nrow(Q1))))
    slope[,j] <- (crossprod(x_grid, gx[,j]))/(sum(x_grid^2))
  }
  print(Sys.time() - begin_time)
  print("Time to do step 1")
  
  ### Step 2: correct the samples using posterior
  begin_time <- Sys.time()
  gz <- matrix(0, nrow = length(all_grid), ncol = ncol(gx))
  gx <- rbind(gx,int,slope)
  for (j in 1:ncol(gx)) {
    gz[,j] <- matrix(x_tilde[,j] - tCA %*% inv_ACA %*% (A %*% x_tilde[,j] - gx[,j]), ncol = 1)
  }
  print(Sys.time() - begin_time)
  print("Time to do step 2")
  gz
}

# ck_method_1(t = z_grid, x_grid, quad_samps = gw, "RW2", 4)

# This method currently is quite slow, and does not seem to provide correct answer... maybe speed up the computation
# using matrix operation at the last step. and try to fix it by really sampling from the unconstrained distribution
ck_method_2 <- function(t, x_grid, quad_samps, GP, mycores = 4, buffer = 0.0001){
  gx <- quad_samps$samps[1:ncol(Q1),]
  thetas <- quad_samps$theta$theta1
  sigmas <- exp(-0.5 * thetas)
  all_grid <- sort(c(t,x_grid))
  t_location <- which(all_grid %in% t)
  n <- length(all_grid)
  all_d <- diff(all_grid)
  H <- compute_H_rue(all_d, n)
  B <- compute_B(all_d, n)
  A <- compute_A(all_d, n)
  if(GP == "RW2"){
    Q <- as(t(H),"dgCMatrix") %*% as(solve(as(A,"dgCMatrix")),"dgCMatrix") %*% as(H,"dgCMatrix")
  }
  else{
    D <- H[-c(1,length(all_grid)),]
    R <- B[-c(1,length(all_grid)), -c(1,length(all_grid))]
    Q <- t(D) %*% solve(R) %*% D
    Q <- as(t(D),"dgCMatrix") %*% as(solve(as(R,"dgCMatrix")),"dgCMatrix") %*% as(D,"dgCMatrix")
  }
  ### Prep Step:
  r <- ncol(Q) - 2
  Q_eigen <- eigen(Q, symmetric = T)
  Q <- Q + buffer * tcrossprod(rep(1, nrow(Q))/sqrt(sum(rep(1, nrow(Q1))))) + buffer * tcrossprod(all_grid/sqrt(sum(all_grid^2)))
  # Q <- Q + buffer * tcrossprod(rep(1, nrow(Q))/(sum(rep(1, nrow(Q1))))) + buffer * tcrossprod(all_grid/(sum(all_grid^2)))
  C <- solve(as(Q,"dgCMatrix"))
  x_location <- which(all_grid %in% x_grid)
  begin_timeA <- Sys.time()
  A <- as(rbind(construct_A(all_grid, x_location), rep(1, nrow(Q))/(sum(rep(1, nrow(Q1)))) , all_grid/(sum(all_grid^2))), "dgCMatrix")
  end_timeA <- Sys.time()
  print(end_timeA - begin_timeA)
  print("Time to construct A")
  begin_time <- Sys.time()
  tCA <- tcrossprod(C,A)
  inv_ACA <- solve(A %*% tCA)
  end_time <- Sys.time()
  print(end_time - begin_time)
  print("Time to compute ACA inverse")
  ### Step 1: sample from unconditional distribution, project to intercept and slope
  begin_time <- Sys.time()
  x_tilde <- matrix(0, nrow = length(all_grid), ncol = ncol(gx))
  int <- matrix(0, nrow = 1, ncol = ncol(gx))
  slope <- matrix(0, nrow = 1, ncol = ncol(gx))
  sd_y <- 1/sqrt(Q_eigen$values[1:r])
  for (i in 1:ncol(gx)) {
    y <- rnorm(n = r, sd = sd_y*sigmas[i])
    x_tilde[,i] <- Q_eigen$vectors[,1:r] %*% y + rnorm(1, sd = sigmas[i]/sqrt(buffer)) * rep(1, nrow(Q))/sqrt(sum(rep(1, nrow(Q1)))) + rnorm(1, sd = sigmas[i]/sqrt(buffer)) * all_grid/sqrt(sum(all_grid^2))
  }
  for (j in 1:ncol(gx)) {
    int[,j] <- (crossprod(rep(1, nrow(Q1)), gx[,j]))/(sum(rep(1, nrow(Q1))))
    slope[,j] <- (crossprod(x_grid, gx[,j]))/(sum(x_grid^2))
  }
  print(mean(int))
  print(mean(slope))
  
  print(Sys.time() - begin_time)
  print("Time to do step 1")
  
  ### Step 2: correct the samples using posterior
  begin_time <- Sys.time()
  # gz <- matrix(0, nrow = length(all_grid), ncol = ncol(gx))
  gx <- rbind(gx,int,slope)
  gz <- x_tilde - tCA %*% inv_ACA %*% (A %*% x_tilde - gx)
  # for (j in 1:ncol(gx)) {
  #   gz[,j] <- matrix(x_tilde[,j] - tCA %*% inv_ACA %*% (A %*% x_tilde[,j] - gx[,j]), ncol = 1)
  # }
  print(Sys.time() - begin_time)
  print("Time to do step 2")
  gz
}







##################### For Full Rank constraining:
eigenspace_const <- function(Q, all_grid, buffer = 10*.Machine$double.eps){
  eigen_decom <- eigen(Q, symmetric = TRUE)
  Q_new <- Q + (abs(eigen_decom$values[(ncol(Q)-1)]) + buffer)*(rep(1,ncol(Q))) %*% t(rep(1,ncol(Q))) + (abs(eigen_decom$values[ncol(Q)]) + buffer) * eigen_decom$values[(ncol(Q))]* all_grid %*% t(all_grid)
  as(Q_new, "dgTMatrix")
}





##############
compute_deriv <- function(vec, order){diff(vec, differences = order)}


