safe_xlogx <- function(x) {
  out <- numeric(length(x))
  idx <- x > 0
  out[idx] <- x[idx] * log(x[idx])
  out
}

log_det_spd <- function(A) {
  R <- chol(A)
  2 * sum(log(diag(R)))
}

E_quad_ih <- function(Yi, Bi, Omega_inv_i, mu_h, Sigma_h) {
  resid_i <- Yi - as.vector(Bi %*% mu_h)
  
  mean_part <- as.numeric(crossprod(resid_i, Omega_inv_i %*% resid_i))
  A_i <- crossprod(Bi, Omega_inv_i %*% Bi)
  trace_part <- sum(diag(A_i %*% Sigma_h))
  
  mean_part + trace_part
}

make_Omega_i <- function(t_i, delta, jitter = 1e-8) {
  D_i <- abs(outer(t_i, t_i, "-"))
  Omega_i <- exp(-delta * D_i)
  diag(Omega_i) <- diag(Omega_i) + jitter
  Omega_i
}

make_Omega_objects <- function(t_list, delta, jitter = 1e-8) {
  N <- length(t_list)
  
  Omega_list <- vector("list", N)
  Omega_inv_list <- vector("list", N)
  log_det_Omega <- numeric(N)
  
  for (i in seq_len(N)) {
    Omega_i <- make_Omega_i(t_list[[i]], delta, jitter)
    R_i <- chol(Omega_i)
    
    Omega_list[[i]] <- Omega_i
    Omega_inv_list[[i]] <- chol2inv(R_i)
    log_det_Omega[i] <- 2 * sum(log(diag(R_i)))
  }
  
  list(
    Omega_list = Omega_list,
    Omega_inv_list = Omega_inv_list,
    log_det_Omega = log_det_Omega
  )
}


############################################################
## OU covariance matrix
############################################################

make_Omega_from_grid <- function(t_grid, delta, sigma = 1) {
  D <- abs(outer(t_grid, t_grid, "-"))
  sigma^2 * exp(-delta * D)
}

