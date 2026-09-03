variational_EM_DP_OU <- function(Y_list,
                                 B_list,
                                 t_list,
                                 H,
                                 m0,
                                 S0,
                                 a0,
                                 b0,
                                 alpha,
                                 init = list(),
                                 max_iter = 200,
                                 tol = 1e-6,
                                 delta_bounds = c(1e-4, 100),
                                 update_delta = TRUE,
                                 jitter = 1e-8,
                                 verbose = TRUE) {
  ############################################################
  ## Basic dimensions
  ############################################################
  
  N <- length(Y_list)
  n_vec <- vapply(Y_list, length, numeric(1))
  M <- length(m0)
  
  if (length(B_list) != N) stop("B_list and Y_list must have the same length.")
  if (length(t_list) != N) stop("t_list and Y_list must have the same length.")
  if (ncol(B_list[[1]]) != M) stop("ncol(B_i) must equal length(m0).")
  if (!all(dim(S0) == c(M, M))) stop("S0 must be M x M.")
  
  ############################################################
  ## Initialization
  ############################################################
  
  delta <- if (!is.null(init$delta)) init$delta else 1
  
  r <- if (!is.null(init$r)) {
    init$r
  } else {
    tmp <- matrix(runif(N * H), nrow = N, ncol = H)
    tmp / rowSums(tmp)
  }
  
  if (!all(dim(r) == c(N, H))) {
    stop("init$r must be an N x H matrix.")
  }
  
  # No initialization of mu_list and Sigma_list is needed because
  # q(phi_h) is updated before mu_list and Sigma_list are used.
  mu_list <- vector("list", H)
  Sigma_list <- vector("list", H)
  
  a_tilde <- if (!is.null(init$a_tilde)) {
    init$a_tilde
  } else {
    rep(a0 + 0.5 * mean(n_vec) * N / H, H)
  }
  
  b_tilde <- if (!is.null(init$b_tilde)) {
    init$b_tilde
  } else {
    rep(b0 + 1, H)
  }
  
  if (length(a_tilde) != H) stop("init$a_tilde must have length H.")
  if (length(b_tilde) != H) stop("init$b_tilde must have length H.")
  
  gamma_1 <- if (!is.null(init$gamma_1)) {
    init$gamma_1
  } else {
    vapply(seq_len(H), function(h) update_gamma_h1(r[, h]), numeric(1))
  }
  
  gamma_2 <- if (!is.null(init$gamma_2)) {
    init$gamma_2
  } else {
    vapply(seq_len(H), function(h) update_gamma_h2(alpha, r, h), numeric(1))
  }
  
  if (length(gamma_1) != H) stop("init$gamma_1 must have length H.")
  if (length(gamma_2) != H) stop("init$gamma_2 must have length H.")
  
  elbo_history <- numeric(max_iter)
  delta_history <- numeric(max_iter)
  
  ############################################################
  ## M-step objective for delta
  ############################################################
  
  objective_delta <- function(delta_candidate,
                              r,
                              E_tau_vec,
                              Y_list,
                              B_list,
                              t_list,
                              mu_list,
                              Sigma_list,
                              jitter) {
    omega_obj <- make_Omega_objects(
      t_list = t_list,
      delta = delta_candidate,
      jitter = jitter
    )
    
    Q <- E_quad(
      Y_list = Y_list,
      B_list = B_list,
      Omega_inv_list = omega_obj$Omega_inv_list,
      mu_list = mu_list,
      Sigma_list = Sigma_list
    )
    
    val <- 0
    
    for (i in seq_len(N)) {
      for (h in seq_len(H)) {
        val <- val + r[i, h] * (
          -0.5 * omega_obj$log_det_Omega[i] -
            0.5 * E_tau_vec[h] * Q[i, h]
        )
      }
    }
    
    val
  }
  
  ############################################################
  ## Variational EM iterations
  ############################################################
  
  elbo_obj <- NULL
  Q <- NULL
  E_tau_vec <- NULL
  E_log_tau_vec <- NULL
  E_log_v_vec <- NULL
  E_log_1_minus_v_vec <- NULL
  E_log_pi_vec <- NULL
  
  for (iter in seq_len(max_iter)) {
    ##########################################################
    ## Build OU covariance objects for current delta
    ##########################################################
    
    omega_obj <- make_Omega_objects(
      t_list = t_list,
      delta = delta,
      jitter = jitter
    )
    
    Omega_inv_list <- omega_obj$Omega_inv_list
    log_det_Omega <- omega_obj$log_det_Omega
    
    ##########################################################
    ## Current expectations
    ##########################################################
    
    E_tau_vec <- E_tau(a_tilde, b_tilde)
    E_log_tau_vec <- E_log_tau(a_tilde, b_tilde)
    E_log_v_vec <- E_log_v(gamma_1, gamma_2)
    E_log_1_minus_v_vec <- E_log_1_minus_v(gamma_1, gamma_2)
    E_log_pi_vec <- E_log_pi(gamma_1, gamma_2)
    
    ##########################################################
    ## Update q(phi_h): Sigma_h and mu_h
    ##########################################################
    
    for (h in seq_len(H)) {
      Sigma_list[[h]] <- update_Sigma_h(
        S0 = S0,
        tau_h_mean = E_tau_vec[h],
        r_h = r[, h],
        B_list = B_list,
        Omega_inv_list = Omega_inv_list
      )
      
      mu_list[[h]] <- update_mu_h(
        Sigma_h = Sigma_list[[h]],
        S0 = S0,
        m0 = m0,
        tau_h_mean = E_tau_vec[h],
        r_h = r[, h],
        B_list = B_list,
        Omega_inv_list = Omega_inv_list,
        Y_list = Y_list
      )
    }
    
    ##########################################################
    ## Update q(tau_h): a_tilde_h and b_tilde_h
    ##########################################################
    
    Q <- E_quad(
      Y_list = Y_list,
      B_list = B_list,
      Omega_inv_list = Omega_inv_list,
      mu_list = mu_list,
      Sigma_list = Sigma_list
    )
    
    for (h in seq_len(H)) {
      a_tilde[h] <- update_a_tilde_h(
        a0 = a0,
        r_h = r[, h],
        n_vec = n_vec
      )
      
      b_tilde[h] <- update_b_tilde_h(
        b0 = b0,
        r_h = r[, h],
        expected_quad_h = Q[, h]
      )
    }
    
    ##########################################################
    ## Update q(v_h): gamma_1 and gamma_2
    ##########################################################
    
    for (h in seq_len(H)) {
      gamma_1[h] <- update_gamma_h1(r[, h])
      gamma_2[h] <- update_gamma_h2(alpha, r, h)
    }
    
    ##########################################################
    ## Refresh expectations after tau and v updates
    ##########################################################
    
    E_tau_vec <- E_tau(a_tilde, b_tilde)
    E_log_tau_vec <- E_log_tau(a_tilde, b_tilde)
    E_log_v_vec <- E_log_v(gamma_1, gamma_2)
    E_log_1_minus_v_vec <- E_log_1_minus_v(gamma_1, gamma_2)
    E_log_pi_vec <- E_log_pi(gamma_1, gamma_2)
    
    ##########################################################
    ## Update q(c_i): responsibilities r_ih
    ##########################################################
    
    r <- update_r(
      E_log_pi = E_log_pi_vec,
      E_log_tau = E_log_tau_vec,
      E_tau = E_tau_vec,
      expected_quad = Q,
      log_det_Omega = log_det_Omega,
      n_vec = n_vec
    )
    
    ##########################################################
    ## M-step: update delta
    ##########################################################
    
    if (update_delta) {
      opt <- optimize(
        f = function(d) {
          objective_delta(
            delta_candidate = d,
            r = r,
            E_tau_vec = E_tau_vec,
            Y_list = Y_list,
            B_list = B_list,
            t_list = t_list,
            mu_list = mu_list,
            Sigma_list = Sigma_list,
            jitter = jitter
          )
        },
        interval = delta_bounds,
        maximum = TRUE
      )
      
      delta <- opt$maximum
    }
    
    ##########################################################
    ## Recompute objects after delta update
    ##########################################################
    
    omega_obj <- make_Omega_objects(
      t_list = t_list,
      delta = delta,
      jitter = jitter
    )
    
    Omega_inv_list <- omega_obj$Omega_inv_list
    log_det_Omega <- omega_obj$log_det_Omega
    
    Q <- E_quad(
      Y_list = Y_list,
      B_list = B_list,
      Omega_inv_list = Omega_inv_list,
      mu_list = mu_list,
      Sigma_list = Sigma_list
    )
    
    ##########################################################
    ## Compute ELBO
    ##########################################################
    
    elbo_obj <- compute_elbo(
      r = r,
      n_vec = n_vec,
      log_det_Omega = log_det_Omega,
      E_quad = Q,
      E_tau = E_tau_vec,
      E_log_tau = E_log_tau_vec,
      E_log_pi = E_log_pi_vec,
      E_log_v = E_log_v_vec,
      E_log_1_minus_v = E_log_1_minus_v_vec,
      mu_list = mu_list,
      Sigma_list = Sigma_list,
      a_tilde = a_tilde,
      b_tilde = b_tilde,
      gamma_1 = gamma_1,
      gamma_2 = gamma_2,
      m0 = m0,
      S0 = S0,
      a0 = a0,
      b0 = b0,
      alpha = alpha
    )
    
    elbo_history[iter] <- elbo_obj$elbo
    delta_history[iter] <- delta
    
    ##########################################################
    ## Convergence check
    ##########################################################
    
    if (verbose) {
      cat(
        sprintf(
          "iter = %d, ELBO = %.6f, delta = %.6f\n",
          iter,
          elbo_history[iter],
          delta
        )
      )
    }
    
    if (iter > 1) {
      elbo_diff <- abs(elbo_history[iter] - elbo_history[iter - 1])
      rel_diff <- elbo_diff / (abs(elbo_history[iter - 1]) + 1e-8)
      
      if (rel_diff < tol) {
        elbo_history <- elbo_history[seq_len(iter)]
        delta_history <- delta_history[seq_len(iter)]
        break
      }
    }
  }
  
  ############################################################
  ## Return fitted object
  ############################################################
  
  cluster <- max.col(r, ties.method = "first")
  
  list(
    r = r,
    cluster = cluster,
    mu_list = mu_list,
    Sigma_list = Sigma_list,
    a_tilde = a_tilde,
    b_tilde = b_tilde,
    gamma_1 = gamma_1,
    gamma_2 = gamma_2,
    delta = delta,
    elbo_history = elbo_history,
    delta_history = delta_history,
    final_expectations = list(
      E_tau = E_tau_vec,
      E_log_tau = E_log_tau_vec,
      E_log_v = E_log_v_vec,
      E_log_1_minus_v = E_log_1_minus_v_vec,
      E_log_pi = E_log_pi_vec,
      E_quad = Q
    ),
    final_elbo_components = elbo_obj$components
  )
}