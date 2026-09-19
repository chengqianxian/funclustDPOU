###### This R script provides functions for each update

### for q(phi_h)
update_Sigma_h <- function(S0, tau_h_mean, r_h, B_list, Omega_inv_list) {
  # S0: M x M prior covariance matrix
  # tau_h_mean: E_q[tau_h] = a_tilde_h / b_tilde_h
  # r_h: length-N vector of responsibilities r_ih for cluster h
  # B_list: list of B_i matrices, each n_i x M
  # Omega_inv_list: list of Omega_i(delta)^(-1), each n_i x n_i
  
  N <- length(B_list)
  M <- ncol(B_list[[1]])
  
  if (length(r_h) != N) {
    stop("length(r_h) must equal length(B_list).")
  }
  if (length(Omega_inv_list) != N) {
    stop("Omega_inv_list must have the same length as B_list.")
  }
  
  S0_inv <- solve(S0)
  precision_h <- S0_inv
  
  for (i in 1:N) {
    Bi <- B_list[[i]]
    Omega_inv_i <- Omega_inv_list[[i]]
    
    precision_h <- precision_h +
      tau_h_mean * r_h[i] * crossprod(Bi, Omega_inv_i %*% Bi)
  }
  
  Sigma_h <- solve(precision_h)
  
  return(Sigma_h)
}

update_mu_h <- function(Sigma_h, S0, m0, tau_h_mean, r_h, B_list, Omega_inv_list, Y_list) {
  # Sigma_h: M x M updated covariance matrix for cluster h
  # S0: M x M prior covariance matrix
  # m0: length-M prior mean vector
  # tau_h_mean: E_q[tau_h] = a_tilde_h / b_tilde_h
  # r_h: length-N vector of responsibilities r_ih for cluster h
  # B_list: list of B_i matrices, each n_i x M
  # Omega_inv_list: list of Omega_i(delta)^(-1), each n_i x n_i
  # Y_list: list of Y_i vectors, each length n_i
  
  N <- length(Y_list)
  
  if (length(B_list) != N) {
    stop("B_list and Y_list must have the same length.")
  }
  if (length(Omega_inv_list) != N) {
    stop("Omega_inv_list and Y_list must have the same length.")
  }
  if (length(r_h) != N) {
    stop("length(r_h) must equal length(Y_list).")
  }
  
  rhs_h <- solve(S0, m0)
  
  for (i in seq_len(N)) {
    Bi <- B_list[[i]]
    Omega_inv_i <- Omega_inv_list[[i]]
    Yi <- Y_list[[i]]
    
    rhs_h <- rhs_h +
      tau_h_mean * r_h[i] * crossprod(Bi, Omega_inv_i %*% Yi)
  }
  
  mu_h <- as.vector(Sigma_h %*% rhs_h)
  
  return(mu_h)
}

### for q(tau_h)
update_a_tilde_h <- function(a0, r_h, n_vec) {
  # a0: prior shape
  # r_h: length-N responsibility vector for cluster h
  # n_vec: length-N vector, where n_vec[i] = length(Y_i)
  
  if (length(r_h) != length(n_vec)) {
    stop("length(r_h) must equal length(n_vec).")
  }
  
  return(a0 + 0.5 * sum(r_h * n_vec))
}
update_b_tilde_h <- function(b0, r_h, expected_quad_h) {
  # b0: prior rate
  # r_h: length-N responsibility vector for cluster h
  # expected_quad_h: length-N vector of Q_ih values for fixed h
  
  if (length(r_h) != length(expected_quad_h)) {
    stop("length(r_h) must equal length(expected_quad_h).")
  }
  
  return(b0 + 0.5 * sum(r_h * expected_quad_h))
}

### for q(c_i)
update_r_i <- function(E_log_pi,
                       E_log_tau,
                       E_tau,
                       expected_quad_i,
                       log_det_Omega_i,
                       n_i) {
  # E_log_pi: length-H vector, E_q[log pi_h]
  # E_log_tau: length-H vector, E_q[log tau_h]
  # E_tau: length-H vector, E_q[tau_h]
  # expected_quad_i: length-H vector, Q_ih for fixed curve i
  # log_det_Omega_i: scalar, log |Omega_i(delta)|
  # n_i: number of observations for curve i
  
  H <- length(E_log_pi)
  
  if (length(E_log_tau) != H) {
    stop("E_log_tau must have the same length as E_log_pi.")
  }
  if (length(E_tau) != H) {
    stop("E_tau must have the same length as E_log_pi.")
  }
  if (length(expected_quad_i) != H) {
    stop("expected_quad_i must have the same length as E_log_pi.")
  }
  
  alpha_i <- E_log_pi +
    0.5 * n_i * E_log_tau -
    0.5 * log_det_Omega_i -
    0.5 * E_tau * expected_quad_i
  
  alpha_i_shifted <- alpha_i - max(alpha_i)
  r_i <- exp(alpha_i_shifted)
  r_i <- r_i / sum(r_i)
  
  return(r_i)
}

update_r <- function(E_log_pi,
                     E_log_tau,
                     E_tau,
                     expected_quad,
                     log_det_Omega,
                     n_vec) {
  # expected_quad: N x H matrix, entry [i, h] = Q_ih
  # log_det_Omega: length-N vector, log |Omega_i(delta)|
  # n_vec: length-N vector of n_i
  
  N <- nrow(expected_quad)
  H <- ncol(expected_quad)
  
  if (length(E_log_pi) != H) stop("length(E_log_pi) must equal ncol(expected_quad).")
  if (length(E_log_tau) != H) stop("length(E_log_tau) must equal ncol(expected_quad).")
  if (length(E_tau) != H) stop("length(E_tau) must equal ncol(expected_quad).")
  if (length(log_det_Omega) != N) stop("length(log_det_Omega) must equal nrow(expected_quad).")
  if (length(n_vec) != N) stop("length(n_vec) must equal nrow(expected_quad).")
  
  r <- matrix(0, nrow = N, ncol = H)
  
  for (i in seq_len(N)) {
    r[i, ] <- update_r_i(
      E_log_pi = E_log_pi,
      E_log_tau = E_log_tau,
      E_tau = E_tau,
      expected_quad_i = expected_quad[i, ],
      log_det_Omega_i = log_det_Omega[i],
      n_i = n_vec[i]
    )
  }
  
  return(r)
}

### for q(v_h), h = 1, ..., H - 1
update_gamma_h1 <- function(r_h) {
  # r_h: length-N vector of responsibilities for stick-breaking component h
  
  return(1 + sum(r_h))
}

update_gamma_h2 <- function(alpha, r, h) {
  # alpha: DP concentration parameter
  # r: N x H responsibility matrix
  # h: stick-breaking index, 1 <= h <= H - 1
  
  H <- ncol(r)
  
  if (H < 2) {
    stop("The truncation level H must be at least 2.")
  }
  if (h < 1 || h > H - 1) {
    stop("h must be between 1 and H - 1, where H = ncol(r).")
  }
  
  # gamma_{h2} = alpha + sum_i sum_{ell > h} r_{i ell};
  # the sum includes component H even though v_H itself is fixed at 1.
  alpha + sum(r[, (h + 1):H, drop = FALSE])
}


