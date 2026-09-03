elbo_loglik <- function(r,
                        n_vec,
                        log_det_Omega,
                        E_log_tau,
                        E_tau,
                        E_quad) {
  # r: N x H responsibility matrix
  # n_vec: length-N vector
  # log_det_Omega: length-N vector, log |Omega_i|
  # E_log_tau: length-H vector
  # E_tau: length-H vector
  # E_quad: N x H matrix, Q_ih
  
  N <- nrow(r)
  H <- ncol(r)
  
  if (length(n_vec) != N) stop("length(n_vec) must equal nrow(r).")
  if (length(log_det_Omega) != N) stop("length(log_det_Omega) must equal nrow(r).")
  if (length(E_log_tau) != H) stop("length(E_log_tau) must equal ncol(r).")
  if (length(E_tau) != H) stop("length(E_tau) must equal ncol(r).")
  if (!all(dim(E_quad) == dim(r))) stop("E_quad must have the same dimension as r.")
  
  val <- 0
  
  for (i in seq_len(N)) {
    for (h in seq_len(H)) {
      term_ih <- -0.5 * n_vec[i] * log(2 * pi) +
        0.5 * n_vec[i] * E_log_tau[h] -
        0.5 * log_det_Omega[i] -
        0.5 * E_tau[h] * E_quad[i, h]
      
      val <- val + r[i, h] * term_ih
    }
  }
  
  val
}


elbo_log_p_c <- function(r, E_log_pi) {
  H <- ncol(r)
  
  if (length(E_log_pi) != H) {
    stop("length(E_log_pi) must equal ncol(r).")
  }
  
  sum(r * matrix(E_log_pi, nrow = nrow(r), ncol = H, byrow = TRUE))
}

elbo_log_p_v <- function(alpha, E_log_1_minus_v) {
  H <- length(E_log_1_minus_v)
  
  sum(
    -lbeta(1, alpha) +
      (alpha - 1) * E_log_1_minus_v
  )
}

elbo_log_p_phi <- function(mu_list, Sigma_list, m0, S0) {
  H <- length(mu_list)
  M <- length(m0)
  
  if (length(Sigma_list) != H) {
    stop("Sigma_list and mu_list must have the same length.")
  }
  
  S0_inv <- solve(S0)
  log_det_S0 <- log_det_spd(S0)
  
  val <- 0
  
  for (h in seq_len(H)) {
    mu_h <- mu_list[[h]]
    Sigma_h <- Sigma_list[[h]]
    
    diff_h <- mu_h - m0
    
    quad_h <- as.numeric(crossprod(diff_h, S0_inv %*% diff_h))
    trace_h <- sum(diag(S0_inv %*% Sigma_h))
    
    val <- val +
      -0.5 * M * log(2 * pi) -
      0.5 * log_det_S0 -
      0.5 * (trace_h + quad_h)
  }
  
  val
}

elbo_log_p_tau <- function(a0, b0, E_tau, E_log_tau) {
  H <- length(E_tau)
  
  if (length(E_log_tau) != H) {
    stop("E_tau and E_log_tau must have the same length.")
  }
  
  sum(
    a0 * log(b0) -
      lgamma(a0) +
      (a0 - 1) * E_log_tau -
      b0 * E_tau
  )
}


elbo_log_q_c <- function(r) {
  sum(safe_xlogx(as.vector(r)))
}

elbo_log_q_v <- function(gamma_1,
                         gamma_2,
                         E_log_v,
                         E_log_1_minus_v) {
  H <- length(gamma_1)
  
  if (length(gamma_2) != H) stop("gamma_2 must have same length as gamma_1.")
  if (length(E_log_v) != H) stop("E_log_v must have same length as gamma_1.")
  if (length(E_log_1_minus_v) != H) {
    stop("E_log_1_minus_v must have same length as gamma_1.")
  }
  
  sum(
    -lbeta(gamma_1, gamma_2) +
      (gamma_1 - 1) * E_log_v +
      (gamma_2 - 1) * E_log_1_minus_v
  )
}

elbo_log_q_phi <- function(Sigma_list) {
  H <- length(Sigma_list)
  M <- nrow(Sigma_list[[1]])
  
  val <- 0
  
  for (h in seq_len(H)) {
    Sigma_h <- Sigma_list[[h]]
    
    val <- val +
      -0.5 * M * log(2 * pi) -
      0.5 * log_det_spd(Sigma_h) -
      0.5 * M
  }
  
  val
}

elbo_log_q_tau <- function(a_tilde,
                           b_tilde,
                           E_tau,
                           E_log_tau) {
  H <- length(a_tilde)
  
  if (length(b_tilde) != H) stop("b_tilde must have same length as a_tilde.")
  if (length(E_tau) != H) stop("E_tau must have same length as a_tilde.")
  if (length(E_log_tau) != H) stop("E_log_tau must have same length as a_tilde.")
  
  sum(
    a_tilde * log(b_tilde) -
      lgamma(a_tilde) +
      (a_tilde - 1) * E_log_tau -
      b_tilde * E_tau
  )
}

compute_elbo <- function(r,
                         n_vec,
                         log_det_Omega,
                         E_quad,
                         E_tau,
                         E_log_tau,
                         E_log_pi,
                         E_log_v,
                         E_log_1_minus_v,
                         mu_list,
                         Sigma_list,
                         a_tilde,
                         b_tilde,
                         gamma_1,
                         gamma_2,
                         m0,
                         S0,
                         a0,
                         b0,
                         alpha) {
  
  term_loglik <- elbo_loglik(
    r = r,
    n_vec = n_vec,
    log_det_Omega = log_det_Omega,
    E_log_tau = E_log_tau,
    E_tau = E_tau,
    E_quad = E_quad
  )
  
  term_log_p_c <- elbo_log_p_c(
    r = r,
    E_log_pi = E_log_pi
  )
  
  term_log_p_v <- elbo_log_p_v(
    alpha = alpha,
    E_log_1_minus_v = E_log_1_minus_v
  )
  
  term_log_p_phi <- elbo_log_p_phi(
    mu_list = mu_list,
    Sigma_list = Sigma_list,
    m0 = m0,
    S0 = S0
  )
  
  term_log_p_tau <- elbo_log_p_tau(
    a0 = a0,
    b0 = b0,
    E_tau = E_tau,
    E_log_tau = E_log_tau
  )
  
  term_log_q_c <- elbo_log_q_c(r)
  
  term_log_q_v <- elbo_log_q_v(
    gamma_1 = gamma_1,
    gamma_2 = gamma_2,
    E_log_v = E_log_v,
    E_log_1_minus_v = E_log_1_minus_v
  )
  
  term_log_q_phi <- elbo_log_q_phi(
    Sigma_list = Sigma_list
  )
  
  term_log_q_tau <- elbo_log_q_tau(
    a_tilde = a_tilde,
    b_tilde = b_tilde,
    E_tau = E_tau,
    E_log_tau = E_log_tau
  )
  
  elbo <- term_loglik +
    term_log_p_c +
    term_log_p_v +
    term_log_p_phi +
    term_log_p_tau -
    term_log_q_c -
    term_log_q_v -
    term_log_q_phi -
    term_log_q_tau
  
  list(
    elbo = elbo,
    components = list(
      loglik = term_loglik,
      log_p_c = term_log_p_c,
      log_p_v = term_log_p_v,
      log_p_phi = term_log_p_phi,
      log_p_tau = term_log_p_tau,
      log_q_c = term_log_q_c,
      log_q_v = term_log_q_v,
      log_q_phi = term_log_q_phi,
      log_q_tau = term_log_q_tau
    )
  )
}

