##### --------------------------------------------------------------------------
############################################################
## Repeat Scenario 1 simulation for multiple seeds, delta = 3
############################################################

############################################################
## Source H-1 truncated stick-breaking implementation
############################################################

source("Helpers.R")
source("expectations.R")
source("update_equations.R")
source("ELBO.R")
source("performance_metrics.R")
source("main_function.R")


run_one_scenario1 <- function(seed,
                              K_true = 3,
                              H = 8,
                              n_per_cluster = 50,
                              n_grid = 100,
                              M = 6,
                              sigma_true = 0.4,
                              delta_true = 3,
                              a0 = 2,
                              b0 = 1,
                              alpha = 1,
                              max_iter = 200,
                              tol = 1e-6,
                              delta_bounds = c(1e-4, 30),
                              active_threshold = 5,
                              verbose = FALSE) {
  ############################################################
  ## Set seed: controls data generation and k-means init
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
  
  r_init <- matrix(1e-3, nrow = N, ncol = H)
  
  for (i in seq_len(N)) {
    r_init[i, km$cluster[i]] <- 1
  }
  
  r_init <- r_init / rowSums(r_init)
  
  gamma_1_init <- vapply(seq_len(H - 1), function(h) {
    update_gamma_h1(r_init[, h])
  }, numeric(1))
  
  gamma_2_init <- vapply(seq_len(H - 1), function(h) {
    update_gamma_h2(alpha, r_init, h)
  }, numeric(1))
  
  a_tilde_init <- rep(a0 + 0.5 * n_grid * N / H, H)
  b_tilde_init <- rep(b0 + 1, H)
  delta_init <- sample(1:5, 1)
  
  ############################################################
  ## Fit model
  ############################################################
  
  fit <- variational_EM_DP_OU(
    Y_list = Y_list,
    B_list = B_list,
    t_list = t_list,
    H = H,
    m0 = m0,
    S0 = S0,
    a0 = a0,
    b0 = b0,
    alpha = alpha,
    init = list(
      r = r_init,
      a_tilde = a_tilde_init,
      b_tilde = b_tilde_init,
      gamma_1 = gamma_1_init,
      gamma_2 = gamma_2_init,
      delta = delta_init
    ),
    max_iter = max_iter,
    tol = tol,
    delta_bounds = delta_bounds,
    update_delta = TRUE,
    jitter = 1e-8,
    verbose = verbose
  )
  
  ############################################################
  ## Metrics and diagnostics
  ############################################################
  
  estimated_label <- fit$cluster
  
  metrics <- clustering_metrics(
    true_label = true_label,
    estimated_label = estimated_label
  )
  
  cluster_size <- colSums(fit$r)
  n_active <- sum(cluster_size > active_threshold)
  
  out <- data.frame(
    seed = seed,
    metrics,
    final_ELBO = tail(fit$elbo_history, 1),
    delta_est = fit$delta,
    delta_true = delta_true,
    n_iter = length(fit$elbo_history),
    n_active = n_active
  )
  
  out
}


############################################################
## Run 50 replications
############################################################

seeds <- 1:50

results_list_1.1 <- vector("list", length(seeds))

start_time.1.1 <- Sys.time()

for (s in seq_along(seeds)) {
  cat("Running seed", seeds[s], "...\n")
  
  results_list_1.1[[s]] <- run_one_scenario1(
    seed = seeds[s],
    verbose = FALSE
  )
}

end_time.1.1 <- Sys.time()

simulation_results_1.1 <- do.call(rbind, results_list_1.1)

time_used.1.1 <- end_time.1.1 - start_time.1.1

simulation_results_1.1

time_used.1.1

############################################################
## Summary over 50 replications 
# (  "final_ELBO",
# "delta_est",
# "n_iter",
# "n_active")
############################################################

metric_cols <- c(
  "Accuracy",
  "V_measure",
  "Rand",
  "Adjusted_Rand",
  "Jaccard"
)

summary_results_1.1 <- data.frame(
  Metric = metric_cols,
  Mean = round(sapply(simulation_results_1.1[metric_cols], mean), 4),
  SD   = round(sapply(simulation_results_1.1[metric_cols], sd), 4),
  Min  = round(sapply(simulation_results_1.1[metric_cols], min), 4),
  Max  = round(sapply(simulation_results_1.1[metric_cols], max), 4)
)

summary_results_1.1

write.csv(
  simulation_results_1.1,
  file = "scenario1_50_replications_results_1.1.csv",
  row.names = FALSE
)


delta_true <- 3

delta_est <- simulation_results_1.1$delta_est

delta_summary <- data.frame(
  Mean = round(mean(delta_est), 4),
  SD = round(sd(delta_est), 4),
  Bias = round(mean(delta_est - delta_true), 4),
  MSE = round(mean((delta_est - delta_true)^2), 4)
)

delta_summary



