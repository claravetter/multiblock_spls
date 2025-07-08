#' Project data onto MBSPLS latent space
#'
#' This function projects data matrices onto the latent space defined by MBSPLS weight vectors.
#' It calculates latent variables and their correlations.
#'
#' @param data List of data matrices to project
#' @param weights List of weight vectors from MBSPLS fit
#' @param correlation_method Method for correlation computation: "pearson", "spearman", or "kendall"
#' @param matrix_norm Norm used for multi-block correlation: "fro", 1, 2, Inf
#' @param force_positive Logical, if TRUE ensures weights sum to positive value for easier interpretation
#'
#' @return A list containing:
#'   \item{rho}{Correlation value or matrix norm of correlations between latent variables}
#'   \item{latent_vars}{Matrix of latent variables}
#'   \item{weights}{Possibly sign-flipped weight vectors}
#'   \item{correlation_matrix}{Full correlation matrix between latent variables}
#'
#' @keywords internal
mbspls_project <- function(data, weights, correlation_method = "pearson", 
                          matrix_norm = "fro", force_positive = FALSE) {
  
  num_blocks <- length(data)
  num_samples <- nrow(data[[1]])
  
  # Calculate latent variables for each block
  latent_vars <- matrix(0, nrow = num_samples, ncol = num_blocks)
  for (i in 1:num_blocks) {
    latent_vars[, i] <- data[[i]] %*% weights[[i]]
  }
  
  # Function to invert weights if needed
  invert_weights <- function(x) -1 * x
  
  # Special case for 2 blocks: ensure positive correlation
  if (num_blocks == 2 && !is.numeric(matrix_norm)) {
    # Use robust_correlations if available
    if (exists("calculate_robust_correlations", mode = "function")) {
      robust_corr <- calculate_robust_correlations(latent_vars, method = correlation_method)
      rho <- robust_corr$corr_matrix[1, 2]
      correlation_matrix <- robust_corr$corr_matrix
    } else {
      # Fallback to standard correlation
      rho <- stats::cor(latent_vars[, 1], latent_vars[, 2], method = correlation_method)
      correlation_matrix <- matrix(c(1, rho, rho, 1), nrow = 2)
    }
    
    # If correlation is negative, flip one weight vector
    if (!is.na(rho) && rho < 0) {
      weights[[2]] <- invert_weights(weights[[2]])
      latent_vars[, 2] <- data[[2]] %*% weights[[2]]
      
      # Recalculate correlation
      if (exists("calculate_robust_correlations", mode = "function")) {
        robust_corr <- calculate_robust_correlations(latent_vars, method = correlation_method)
        rho <- robust_corr$corr_matrix[1, 2]
        correlation_matrix <- robust_corr$corr_matrix
      } else {
        rho <- stats::cor(latent_vars[, 1], latent_vars[, 2], method = correlation_method)
        correlation_matrix <- matrix(c(1, rho, rho, 1), nrow = 2)
      }
    }
    
    # If requested, ensure sum of weights is positive for easier interpretation
    if (force_positive) {
      for (i in 1:2) {
        if (sum(weights[[i]]) < 0) {
          weights[[i]] <- invert_weights(weights[[i]])
          latent_vars[, i] <- data[[i]] %*% weights[[i]]
        }
      }
      # Recalculate correlation with the final weights
      if (exists("calculate_robust_correlations", mode = "function")) {
        robust_corr <- calculate_robust_correlations(latent_vars, method = correlation_method)
        rho <- robust_corr$corr_matrix[1, 2]
        correlation_matrix <- robust_corr$corr_matrix
      } else {
        rho <- stats::cor(latent_vars[, 1], latent_vars[, 2], method = correlation_method)
        correlation_matrix <- matrix(c(1, rho, rho, 1), nrow = 2)
      }
    }
  } else {
    # For more than 2 blocks or when using matrix norm
    
    # If requested, ensure sum of weights is positive
    if (force_positive) {
      for (i in 1:num_blocks) {
        if (sum(weights[[i]]) < 0) {
          weights[[i]] <- invert_weights(weights[[i]])
          latent_vars[, i] <- data[[i]] %*% weights[[i]]
        }
      }
    }
    
    # Calculate robust correlations if the function is available
    if (exists("calculate_robust_correlations", mode = "function")) {
      # Use the robust correlation calculation function
      robust_results <- calculate_robust_correlations(latent_vars, method = correlation_method)
      corr_matrix <- robust_results$corr_matrix
      
      # Get matrix norm from robust calculation
      if (!is.na(robust_results$matrix_norm)) {
        rho <- robust_results$matrix_norm
      } else {
        # Fallback to calculating norm manually
        if (identical(matrix_norm, "fro")) {
          rho <- norm(corr_matrix, "F")  # Frobenius norm
        } else if (is.numeric(matrix_norm)) {
          matrix_norm_str <- as.character(matrix_norm)
          rho <- norm(corr_matrix, matrix_norm_str)
        } else {
          rho <- norm(corr_matrix, matrix_norm)
        }
      }
    } else {
      # Fallback to original implementation
      # Handle special case with single-feature blocks
      var_check <- apply(latent_vars, 2, stats::var)
      if (any(var_check == 0 | is.na(var_check))) {
        # Create a correlation matrix with NAs for zero-variance columns
        corr_matrix <- matrix(NA, nrow = num_blocks, ncol = num_blocks)
        
        # Fill in valid correlations where possible
        valid_cols <- which(var_check > 0 & !is.na(var_check))
        if (length(valid_cols) > 1) {
          # Calculate correlations for valid columns only
          valid_corr <- stats::cor(latent_vars[, valid_cols, drop = FALSE], method = correlation_method)
          corr_matrix[valid_cols, valid_cols] <- valid_corr
        }
        
        # Set diagonal to 1
        diag(corr_matrix) <- 1
      } else {
        # Normal case - calculate correlation matrix between all latent variables
        corr_matrix <- stats::cor(latent_vars, method = correlation_method)
      }
      
      # Calculate matrix norm as a measure of overall correlation
      if (identical(matrix_norm, "fro")) {
        rho <- norm(corr_matrix, "F")  # Frobenius norm
      } else if (is.numeric(matrix_norm)) {
        # For numeric norms, convert to string first as required by base::norm
        matrix_norm_str <- as.character(matrix_norm)
        rho <- norm(corr_matrix, matrix_norm_str)  # Other norms: "1", "2", "Inf"
      } else {
        rho <- norm(corr_matrix, matrix_norm)  # Other norms as character
      }
    }
    
    # Store the correlation matrix for later use
    correlation_matrix <- corr_matrix
  }
  
  # Return results
  return(list(
    rho = rho,
    latent_vars = latent_vars,
    weights = weights,
    correlation_matrix = correlation_matrix
  ))
}
