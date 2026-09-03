############################################################
## Performance metric functions for clustering
############################################################

choose2 <- function(x) {
  x * (x - 1) / 2
}

contingency_table <- function(true_label, estimated_label) {
  true_label <- as.factor(true_label)
  estimated_label <- as.factor(estimated_label)
  table(true_label, estimated_label)
}

############################################################
## 1. Accuracy after optimal label matching
############################################################

max_assignment_sum <- function(W) {
  # Maximizes sum of assigned entries in a square or rectangular matrix W.
  # Pads W to square and uses dynamic programming over subsets.
  
  W <- as.matrix(W)
  nr <- nrow(W)
  nc <- ncol(W)
  n <- max(nr, nc)
  
  W_pad <- matrix(0, nrow = n, ncol = n)
  W_pad[seq_len(nr), seq_len(nc)] <- W
  
  # dp[[mask + 1]] stores best value after assigning some rows to columns in mask
  dp <- rep(-Inf, 2^n)
  dp[1] <- 0
  
  for (i in seq_len(n)) {
    new_dp <- rep(-Inf, 2^n)
    
    for (mask in 0:(2^n - 1)) {
      if (!is.finite(dp[mask + 1])) next
      
      used <- as.logical(intToBits(mask))[seq_len(n)]
      available_cols <- which(!used)
      
      for (j in available_cols) {
        new_mask <- bitwOr(mask, bitwShiftL(1, j - 1))
        new_dp[new_mask + 1] <- max(
          new_dp[new_mask + 1],
          dp[mask + 1] + W_pad[i, j]
        )
      }
    }
    
    dp <- new_dp
  }
  
  dp[2^n]
}

clustering_accuracy <- function(true_label, estimated_label) {
  tab <- contingency_table(true_label, estimated_label)
  correct <- max_assignment_sum(tab)
  correct / length(true_label)
}

############################################################
## 2. Pair-count quantities: TP, TN, FP, FN
############################################################

pair_counts <- function(true_label, estimated_label) {
  tab <- contingency_table(true_label, estimated_label)
  
  n <- sum(tab)
  
  TP <- sum(choose2(tab))
  true_pairs <- sum(choose2(rowSums(tab)))
  estimated_pairs <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(n)
  
  FP <- estimated_pairs - TP
  FN <- true_pairs - TP
  TN <- total_pairs - TP - FP - FN
  
  list(
    TP = as.numeric(TP),
    TN = as.numeric(TN),
    FP = as.numeric(FP),
    FN = as.numeric(FN),
    total_pairs = as.numeric(total_pairs)
  )
}

############################################################
## 3. Rand index
############################################################

rand_index <- function(true_label, estimated_label) {
  pc <- pair_counts(true_label, estimated_label)
  
  (pc$TP + pc$TN) / pc$total_pairs
}

############################################################
## 4. Adjusted Rand index
############################################################

adjusted_rand_index <- function(true_label, estimated_label) {
  tab <- contingency_table(true_label, estimated_label)
  
  n <- sum(tab)
  
  sum_comb <- sum(choose2(tab))
  row_comb <- sum(choose2(rowSums(tab)))
  col_comb <- sum(choose2(colSums(tab)))
  total_comb <- choose2(n)
  
  expected_index <- row_comb * col_comb / total_comb
  max_index <- 0.5 * (row_comb + col_comb)
  
  denom <- max_index - expected_index
  
  if (abs(denom) < .Machine$double.eps) {
    return(0)
  }
  
  (sum_comb - expected_index) / denom
}

############################################################
## 5. Jaccard index
############################################################

jaccard_index <- function(true_label, estimated_label) {
  pc <- pair_counts(true_label, estimated_label)
  
  denom <- pc$TP + pc$FP + pc$FN
  
  if (denom == 0) {
    return(0)
  }
  
  pc$TP / denom
}

############################################################
## 6. V-measure
############################################################

entropy_from_counts <- function(counts) {
  counts <- counts[counts > 0]
  probs <- counts / sum(counts)
  -sum(probs * log(probs))
}

v_measure <- function(true_label, estimated_label, beta = 1) {
  tab <- contingency_table(true_label, estimated_label)
  n <- sum(tab)
  
  # Entropy of true classes and estimated clusters
  H_true <- entropy_from_counts(rowSums(tab))
  H_est <- entropy_from_counts(colSums(tab))
  
  # H(true | estimated)
  H_true_given_est <- 0
  for (j in seq_len(ncol(tab))) {
    n_j <- sum(tab[, j])
    if (n_j > 0) {
      probs <- tab[, j] / n_j
      probs <- probs[probs > 0]
      H_true_given_est <- H_true_given_est -
        (n_j / n) * sum(probs * log(probs))
    }
  }
  
  # H(estimated | true)
  H_est_given_true <- 0
  for (i in seq_len(nrow(tab))) {
    n_i <- sum(tab[i, ])
    if (n_i > 0) {
      probs <- tab[i, ] / n_i
      probs <- probs[probs > 0]
      H_est_given_true <- H_est_given_true -
        (n_i / n) * sum(probs * log(probs))
    }
  }
  
  homogeneity <- if (H_true == 0) 1 else 1 - H_true_given_est / H_true
  completeness <- if (H_est == 0) 1 else 1 - H_est_given_true / H_est
  
  if (homogeneity + completeness == 0) {
    V <- 0
  } else {
    V <- (1 + beta) * homogeneity * completeness /
      (beta * homogeneity + completeness)
  }
  
  list(
    v_measure = V,
    homogeneity = homogeneity,
    completeness = completeness
  )
}

############################################################
## 7. Wrapper: all metrics
############################################################

clustering_metrics <- function(true_label, estimated_label) {
  pc <- pair_counts(true_label, estimated_label)
  vm <- v_measure(true_label, estimated_label)
  
  data.frame(
    Accuracy = clustering_accuracy(true_label, estimated_label),
    V_measure = vm$v_measure,
    Rand = rand_index(true_label, estimated_label),
    Adjusted_Rand = adjusted_rand_index(true_label, estimated_label),
    Jaccard = jaccard_index(true_label, estimated_label)
  )
}
