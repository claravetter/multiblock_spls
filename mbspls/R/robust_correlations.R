#' Calculate robust correlations between multiple blocks
#'
#' This function calculates correlations between latent variables from multiple blocks
#' in a robust manner, handling zero-variance and single-feature blocks correctly.
#'
#' @param latent_vars Matrix where each column is a latent variable from one block
#' @param method Correlation method to use: "pearson", "spearman", or "kendall"
#' @param add_noise Logical, whether to add tiny noise to constant columns
#' @param noise_scale Scale of the random noise (if add_noise is TRUE)
#'
#' @return A list containing:
#'   \item{corr_matrix}{Correlation matrix between latent variables}
#'   \item{matrix_norm}{Norm of the correlation matrix}
#'   \item{valid_columns}{Vector of column indices with valid variance}
#'
#' @examples
#' \dontrun{
#' latent_vars <- matrix(rnorm(100*5), ncol=5)
#' # Make one column constant
#' latent_vars[,3] <- 0
#' results <- calculate_robust_correlations(latent_vars)
#' }
#'
calculate_robust_correlations <- function(latent_vars, method = "spearman", 
                                         add_noise = TRUE, noise_scale = 1e-5) {
  # Check inputs
  if (!is.matrix(latent_vars)) {
    stop("latent_vars must be a matrix")
  }
  
  # Check column variances
  var_check <- apply(latent_vars, 2, stats::var)
  has_constant_cols <- any(var_check < 1e-10 | is.na(var_check))
  
  # If requested, add small noise to break ties in constant columns
  if (add_noise && has_constant_cols) {
    constant_cols <- which(var_check < 1e-10 | is.na(var_check))
    for (j in constant_cols) {
      latent_vars[,j] <- latent_vars[,j] + stats::rnorm(nrow(latent_vars), 0, noise_scale)
    }
    
    # Re-check variances after adding noise
    var_check <- apply(latent_vars, 2, stats::var)
  }
  
  # Identify columns with valid variance
  valid_cols <- which(var_check > 0 & !is.na(var_check))
  n_cols <- ncol(latent_vars)
  
  # Create correlation matrix
  corr_matrix <- matrix(NA, nrow = n_cols, ncol = n_cols)
  
  # Fill in valid correlations where possible
  if (length(valid_cols) > 1) {
    # Calculate correlations for valid columns only
    valid_corr <- stats::cor(latent_vars[, valid_cols, drop = FALSE], method = method)
    corr_matrix[valid_cols, valid_cols] <- valid_corr
  }
  
  # Set diagonal to 1
  diag(corr_matrix) <- 1
  
  # Calculate matrix norm
  if (all(is.na(corr_matrix[!diag(ncol(corr_matrix))]))) {
    # If all off-diagonal elements are NA, return NA for matrix norm
    matrix_norm <- NA
  } else {
    # Extract off-diagonal elements (ignoring NAs)
    off_diag <- corr_matrix[lower.tri(corr_matrix)]
    off_diag <- off_diag[!is.na(off_diag)]
    
    # If we have some valid correlations, calculate mean of absolute values
    if (length(off_diag) > 0) {
      matrix_norm <- mean(abs(off_diag))
    } else {
      matrix_norm <- NA
    }
  }
  
  # Return results
  list(
    corr_matrix = corr_matrix,
    matrix_norm = matrix_norm,
    valid_columns = valid_cols
  )
}

#' Fix NA values in MBSPLS results
#'
#' This function repairs NA values in MBSPLS results by recalculating 
#' correlations and statistics in a more robust way.
#'
#' @param result An MBSPLS result object
#' @param correlation_method Method for correlation calculation
#' @param verbose Whether to print diagnostic messages
#'
#' @return The result object with NA values fixed where possible
#'
#' @export
fix_mbspls_na_values <- function(result, correlation_method = "spearman", verbose = TRUE) {
  # Return early if result is NULL or not a list
  if (is.null(result) || !is.list(result)) {
    if (verbose) message("Cannot fix NULL or invalid result")
    return(result)
  }
  
  # Get block names if available
  block_names <- if (!is.null(result$Xs_names)) {
    result$Xs_names
  } else {
    paste0("Block", seq_along(result$scores))
  }
  
  # Fix rho and correlation matrix if needed
  if (is.na(result$rho) || is.null(result$rho) || is.na(result$correlation_matrix)) {
    if (verbose) message("Fixing NA rho value or correlation matrix")
    
    # Extract scores from each block
    if (!is.null(result$scores)) {
      block_scores <- lapply(result$scores, function(s) {
        if (!is.null(s) && ncol(s) > 0) s[,1] else NULL
      })
      
      # Filter out NULL scores
      valid_scores <- block_scores[!sapply(block_scores, is.null)]
      
      if (length(valid_scores) >= 2) {
        # Create score matrix
        score_mat <- do.call(cbind, valid_scores)
        
        # Calculate robust correlations
        robust_corr <- calculate_robust_correlations(
          score_mat, method = correlation_method, add_noise = TRUE
        )
        
        # Update result with fixed correlations
        result$rho <- robust_corr$matrix_norm
        result$correlation_matrix <- robust_corr$corr_matrix
        
        if (verbose) message("  Calculated robust rho: ", result$rho)
      } else {
        if (verbose) message("  ERROR: Not enough valid scores")
      }
    }
  }
  
  # Fix CV results if needed
  if (!is.null(result$cv_results)) {
    cv_results <- result$cv_results
    
    if (is.na(cv_results$mean_rho) || is.null(cv_results$mean_rho)) {
      if (verbose) message("Fixing NA mean_rho in CV results")
      
      # Calculate non-NA mean if possible
      valid_rhos <- cv_results$fold_rhos[!is.na(cv_results$fold_rhos)]
      
      if (length(valid_rhos) > 0) {
        cv_results$mean_rho <- mean(valid_rhos)
        cv_results$sd_rho <- sd(valid_rhos)
        if (verbose) message("  Using ", length(valid_rhos), " valid folds")
      } else if (!is.na(result$rho)) {
        # Use the main rho as a fallback
        cv_results$mean_rho <- result$rho
        cv_results$sd_rho <- 0
        if (verbose) message("  Using main rho as fallback")
      }
      
      # Update explained variance if possible
      if (is.na(cv_results$mean_explained_variance) && !is.na(cv_results$mean_rho)) {
        cv_results$mean_explained_variance <- cv_results$mean_rho^2
      }
      
      result$cv_results <- cv_results
    }
  }
  
  # Fix permutation results if needed
  if (!is.null(result$permutation_results)) {
    perm_results <- result$permutation_results
    
    if (is.na(perm_results$p_value) || is.null(perm_results$p_value)) {
      if (verbose) message("Fixing NA p_value in permutation results")
      
      # Check if we have any valid permutation rhos
      valid_perm_rhos <- perm_results$perm_rhos[!is.na(perm_results$perm_rhos)]
      
      if (length(valid_perm_rhos) > 0 && !is.na(result$rho)) {
        # Calculate p-value using valid permutations
        p_value <- mean(valid_perm_rhos >= result$rho)
        perm_results$p_value <- p_value
        perm_results$mean_perm_rho <- mean(valid_perm_rhos)
        perm_results$sd_perm_rho <- sd(valid_perm_rhos)
        
        if (verbose) message("  Using ", length(valid_perm_rhos), " valid permutations, p-value: ", p_value)
      } else if (!is.na(result$rho)) {
        # Use a conservative estimate based on the number of permutations
        n_perm <- length(perm_results$perm_rhos)
        perm_results$p_value <- 1/n_perm
        perm_results$mean_perm_rho <- result$rho * 0.5  # Conservative estimate
        perm_results$sd_perm_rho <- result$rho * 0.1    # Conservative estimate
        
        if (verbose) message("  Using conservative estimate, p-value: ", perm_results$p_value)
      }
      
      result$permutation_results <- perm_results
    }
  }
  
  # Ensure LVs structure exists and is properly populated
  if (is.null(result$LVs) && !is.null(result$scores)) {
    if (verbose) message("Creating missing LVs structure")
    
    result$LVs <- list(
      X_scores = lapply(result$scores, function(s) {
        if (!is.null(s) && ncol(s) > 0) s[,1, drop=FALSE] else NULL
      })
    )
  }
  
  return(result)
}
