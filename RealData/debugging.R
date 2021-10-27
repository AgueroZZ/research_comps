
ck_method_2 <- function(t, x_grid, quad_samps, GP, mycores = 4, buffer = 1){
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
  Q_eigen$values[(r+1):ncol(Q)] <- 0
  # round_tol <- min(abs(Q_eigen$values[1:r]))
  Q <- Q_eigen$vectors %*% Diagonal(x = Q_eigen$values) %*% t(Q_eigen$vectors) + buffer * tcrossprod(Q_eigen$vectors[,(r+1)]) + buffer * tcrossprod(Q_eigen$vectors[,(r+2)])
  
  # Q <- Q + buffer * tcrossprod(rep(1, nrow(Q))/(sum(rep(1, nrow(Q1))))) + buffer * tcrossprod(all_grid/(sum(all_grid^2)))

  C <- solve(as(Q,"dgCMatrix"))
  
  x_location <- which(all_grid %in% x_grid)
  begin_timeA <- Sys.time()
  A <- as(rbind(construct_A(all_grid, x_location), Q_eigen$vectors[,(r+1)], Q_eigen$vectors[,(r+2)]), "dgCMatrix")
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
    x_tilde[,i] <- Q_eigen$vectors[,1:r] %*% y + rnorm(1, sd = sigmas[i]/sqrt(buffer)) * Q_eigen$vectors[,(r+1)] + rnorm(1, sd = sigmas[i]/sqrt(buffer)) * Q_eigen$vectors[,(r+2)]
  }
  # return(list(samples = x_tilde, Q = Q))
  # for (j in 1:ncol(gx)) {
  #   int[,j] <- (crossprod(rep(1, nrow(Q1)), gx[,j]))/(sum(rep(1, nrow(Q1))))
  #   slope[,j] <- (crossprod(x_grid, gx[,j]))/(sum(x_grid^2))
  # }
  # print(mean(int))
  # print(mean(slope))
  
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




start_time <- Sys.time()
gz <- ck_method_2(t = z_grid, x_grid, quad_samps = gw, "RW2", 4)
Sys.time() - start_time


gz_mean <- apply(gz, 1, mean)
plot(gz_mean ~ sort(c(x_grid,z_grid)), type = 'l')
