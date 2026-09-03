############################################################
## Repeat Scenario 1 simulation for multiple seeds, delta = 5
## MCMC version
############################################################

run_one_scenario1_mcmc <- function(seed,
                                   K_true = 3,
                                   H = 8,
                                   n_per_cluster = 50,
                                   n_grid = 100,
                                   M = 6,
                                   sigma_true = 0.4,
                                   delta_true = 5,
                                   a0 = 2,
                                   b0 = 1,
                                   alpha = 1,
                                   MCMC_iter = 5000,
                                   burn_in = 1000,
                                   s_delta = 0.1,
                                   delta_bounds = c(1e-4, 30),
                                   active_threshold = 5,
                                   verbose = FALSE) {
  ############################################################
  ## Set seed: controls data generation and initialization
  ############################################################
  
  set.seed(seed)
  
  ############################################################
  ## Simulation settings
  ############################################################
  
  N <- K_true * n_per_cluster
  t_grid <- seq(0, 1, length.out = n_grid)
  
  phi_true <- rbind(
    c(1.5, 1.0, 1.6, 1.8, 1.0, 1.5),
    c(1.8, 0.6, 0.4, 2.6, 2.8, 1.6),
    c(1.2, 1.8, 2.2, 0.8, 0.6, 1.8)
  )
  
  B <- splines::bs(
    x = t_grid,
    df = M,
    degree = 3,
    intercept = TRUE
  )
  
  B <- as.matrix(B)
  
  ############################################################
  ## OU covariance matrix
  ############################################################
  
  Omega_true <- make_Omega_from_grid(
    t_grid = t_grid,
    delta = delta_true,
    sigma = sigma_true
  )
  
  ############################################################
  ## Generate data
  ############################################################
  
  Y_list <- vector("list", N)
  B_list <- vector("list", N)
  t_list <- vector("list", N)
  true_label <- integer(N)
  
  idx <- 1
  
  for (k in seq_len(K_true)) {
    mean_k <- as.vector(B %*% phi_true[k, ])
    
    for (i in seq_len(n_per_cluster)) {
      eps_i <- as.vector(MASS::mvrnorm(
        n = 1,
        mu = rep(0, n_grid),
        Sigma = Omega_true
      ))
      
      Y_list[[idx]] <- mean_k + eps_i
      B_list[[idx]] <- B
      t_list[[idx]] <- t_grid
      true_label[idx] <- k
      
      idx <- idx + 1
    }
  }
  
  Y_mat <- do.call(rbind, Y_list)
  
  ############################################################
  ## Hyperparameters
  ############################################################
  
  m0 <- rep(0, M)
  S0 <- diag(100, M)
  
  ############################################################
  ## Initialization by k-means
  ############################################################
  
  km <- kmeans(Y_mat, centers = K_true, nstart = 20)
  
  c_init <- km$cluster
  
  # Optional but useful: initialize tau around prior mean
  tau_init <- rep(a0 / b0, H)
  
  # Random initial delta, same spirit as your VBEM code
  delta_init <- sample(3:7, 1)
  
  ############################################################
  ## Fit MCMC model
  ############################################################
  
  fit <- mcmc_DP_OU(
    Y_list = Y_list,
    B_list = B_list,
    t_list = t_list,
    H = H,
    m0 = m0,
    S0 = S0,
    a0 = a0,
    b0 = b0,
    alpha = alpha,
    MCMC_iter = MCMC_iter,
    burn_in = burn_in,
    init = list(
      c = c_init,
      tau = tau_init,
      delta = delta_init
    ),
    s_delta = s_delta,
    delta_bounds = delta_bounds,
    jitter = 1e-8,
    verbose = verbose
  )
  
  ############################################################
  ## Posterior clustering estimate
  ############################################################
  
  posterior_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  estimated_label <- apply(fit$c_samples, 2, posterior_mode)
  
  ############################################################
  ## Metrics and diagnostics
  ############################################################
  
  metrics <- clustering_metrics(
    true_label = true_label,
    estimated_label = estimated_label
  )
  
  cluster_size <- tabulate(estimated_label, nbins = H)
  n_active <- sum(cluster_size > active_threshold)
  
  delta_est <- mean(fit$delta_samples)
  delta_sd <- sd(fit$delta_samples)
  
  out <- data.frame(
    seed = seed,
    metrics,
    delta_est = delta_est,
    delta_sd = delta_sd,
    delta_true = delta_true,
    delta_acceptance_rate = fit$delta_acceptance_rate,
    MCMC_iter = MCMC_iter,
    burn_in = burn_in,
    n_active = n_active
  )
  
  out
}

############################################################
## Run 50 MCMC replications
############################################################

seeds <- 1:50

results_list_mcmc_1.2 <- vector("list", length(seeds))

start_time_mcmc_1.2 <- Sys.time()

for (s in seq_along(seeds)) {
  cat("Running seed", seeds[s], "...\n")
  
  results_list_mcmc_1.2[[s]] <- run_one_scenario1_mcmc(
    seed = seeds[s],
    MCMC_iter = 5000,
    burn_in = 1000,
    s_delta = 0.1,
    verbose = FALSE
  )
}

end_time_mcmc_1.2 <- Sys.time()

simulation_results_mcmc_1.2 <- do.call(rbind, results_list_mcmc_1.2)

time_used_mcmc_1.2 <- end_time_mcmc_1.2 - start_time_mcmc_1.2

simulation_results_mcmc_1.2
time_used_mcmc_1.2

############################################################
## Summary over 50 MCMC replications
############################################################

metric_cols <- c(
  "Accuracy",
  "V_measure",
  "Rand",
  "Adjusted_Rand",
  "Jaccard"
)

summary_results_mcmc_1.2 <- data.frame(
  Metric = metric_cols,
  Mean = round(sapply(simulation_results_mcmc_1.2[metric_cols], mean), 4),
  SD   = round(sapply(simulation_results_mcmc_1.2[metric_cols], sd), 4),
  Min  = round(sapply(simulation_results_mcmc_1.2[metric_cols], min), 4),
  Max  = round(sapply(simulation_results_mcmc_1.2[metric_cols], max), 4)
)

summary_results_mcmc_1.2


############################################################
## Delta summary
############################################################

delta_true <- 5

delta_est <- simulation_results_mcmc_1.2$delta_est

delta_summary_mcmc_1.2 <- data.frame(
  Mean = round(mean(delta_est), 4),
  SD = round(sd(delta_est), 4),
  Bias = round(mean(delta_est - delta_true), 4),
  MSE = round(mean((delta_est - delta_true)^2), 4)
)

delta_summary_mcmc_1.2

write.csv(
  simulation_results_mcmc_1.2,
  file = "scenario1_50_replications_results_mcmc_1.2.csv",
  row.names = FALSE
)
