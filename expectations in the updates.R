### This R manuscript provides calculations for all the expectations
### all are one-time calculations, that is, result is a vector or a matrix
E_tau <- function(a_tilde, b_tilde) {
  if (length(a_tilde) != length(b_tilde)) {
    stop("a_tilde and b_tilde must have the same length.")
  }
  
  a_tilde / b_tilde
}

E_log_tau <- function(a_tilde, b_tilde) {
  if (length(a_tilde) != length(b_tilde)) {
    stop("a_tilde and b_tilde must have the same length.")
  }
  
  digamma(a_tilde) - log(b_tilde)
}

## for all curves and all clusters, N \times H
E_quad <- function(Y_list, B_list, Omega_inv_list, mu_list, Sigma_list) {
  # Y_list: list of Y_i
  # B_list: list of B_i
  # Omega_inv_list: list of Omega_i^{-1}
  # mu_list: list of mu_h
  # Sigma_list: list of Sigma_h
  
  N <- length(Y_list)
  H <- length(mu_list)
  
  if (length(B_list) != N) {
    stop("B_list and Y_list must have the same length.")
  }
  if (length(Omega_inv_list) != N) {
    stop("Omega_inv_list and Y_list must have the same length.")
  }
  if (length(Sigma_list) != H) {
    stop("Sigma_list and mu_list must have the same length.")
  }
  
  Q <- matrix(0, nrow = N, ncol = H)
  
  for (i in 1:N) {
    for (h in 1:H) {
      Q[i, h] <- E_quad_ih(
        Yi = Y_list[[i]],
        Bi = B_list[[i]],
        Omega_inv_i = Omega_inv_list[[i]],
        mu_h = mu_list[[h]],
        Sigma_h = Sigma_list[[h]]
      )
    }
  }
  
  Q
}

E_log_v <- function(gamma_1, gamma_2) {
  if (length(gamma_1) != length(gamma_2)) {
    stop("gamma_1 and gamma_2 must have the same length.")
  }
  
  digamma(gamma_1) - digamma(gamma_1 + gamma_2)
}

E_log_1_minus_v <- function(gamma_1, gamma_2) {
  if (length(gamma_1) != length(gamma_2)) {
    stop("gamma_1 and gamma_2 must have the same length.")
  }
  
  digamma(gamma_2) - digamma(gamma_1 + gamma_2)
}

E_log_pi <- function(gamma_1, gamma_2) {
  H <- length(gamma_1)
  
  if (length(gamma_2) != H) {
    stop("gamma_1 and gamma_2 must have the same length.")
  }
  
  Elogv <- E_log_v(gamma_1, gamma_2)
  Elog1mv <- E_log_1_minus_v(gamma_1, gamma_2)
  
  Elogpi <- numeric(H)
  
  for (h in 1:H) {
    if (h == 1) {
      Elogpi[h] <- Elogv[h]
    } else {
      Elogpi[h] <- Elogv[h] + sum(Elog1mv[1:(h - 1)])
    }
  }
  
  Elogpi
}

