#' Core MBSPLS algorithm implementation
#'
#' This function implements the Multi-Block Sparse Partial Least Squares algorithm.
#' It extracts weight vectors that maximize the covariance between multiple data blocks
#' subject to sparsity constraints.
#'
#' @param matrices List of data matrices, each representing a block.
#' @param sparsity Numeric vector of sparsity parameters for each block.
#' @param g_factors Matrix specifying the type of optimization.
#' @param convergence_threshold Numeric threshold for convergence. Default: 1e-5.
#' @param max_iterations Maximum number of iterations for the algorithm. Default: 1000.
#' @param print_level Integer controlling verbosity. Default: 0.
#' @param initial_weights Optional list of initial weight vectors.
#' @param ensure_nonzero_blocks Logical, if TRUE, ensures each block has at least one non-zero feature. Default: FALSE.
#'
#' @return A list containing:
#'   \item{weights}{List of weight vectors for each matrix}
#'   \item{covariances}{List of covariance matrices between each pair of matrices}
#'   \item{latent_vars}{Latent variable matrices}
#'   \item{converged}{Logical indicating whether the algorithm converged}
#'   \item{iterations}{Number of iterations until convergence}
#'   \item{convergence_diff}{Final convergence difference value}
#'   \item{dimensions}{Dimensions of each block}
#'
#' @keywords internal
mbspls_fit_core <- function(matrices, 
                           sparsity, 
                           g_factors, 
                           convergence_threshold = 1e-5, 
                           max_iterations = 1000,
                           print_level = 0,
                           ensure_nonzero_blocks = FALSE,
                           initial_weights = NULL) {
  
  # Initial checks
  num_matrices <- length(matrices)
  
  # Initialize cell arrays
  weights <- vector("list", num_matrices)
  covariances <- vector("list", num_matrices * num_matrices)
  dim(covariances) <- c(num_matrices, num_matrices)
  
  # Check sparsity parameters and number of features in each matrix
  num_weights <- numeric(num_matrices)
  no_sparse_matrix <- logical(num_matrices)
  
  for (i in 1:num_matrices) {
    num_weights[i] <- ncol(matrices[[i]])
    if (sparsity[i] < 1 || sparsity[i] > sqrt(num_weights[i])) {
      no_sparse_matrix[i] <- TRUE
      if (print_level > -1) {
        warning(sprintf("Sparsity parameter for matrix %d is out of interval: 1 <= c <= sqrt(number_features). Not using sparsity on this matrix.", i))
      }
    } else {
      no_sparse_matrix[i] <- FALSE
    }
  }
  
  # Check the dimensions of all matrices (samples should match)
  dims <- t(sapply(matrices, dim))
  if (length(unique(dims[, 1])) > 1) {
    stop("All matrices must have the same number of samples (rows).")
  }
  
  # Initialize weight vectors (random or provided)
  if (is.null(initial_weights)) {
    # Random initialization (set seed for reproducibility)
    set.seed(42)
    for (i in 1:num_matrices) {
      weights[[i]] <- rnorm(num_weights[i])
      # Normalize the weight vector
      weights[[i]] <- weights[[i]] / sqrt(sum(weights[[i]]^2))
    }
  } else {
    # Use provided initial weights
    weights <- initial_weights
    # Check and normalize weights if needed
    for (i in 1:num_matrices) {
      if (length(weights[[i]]) != num_weights[i]) {
        stop(sprintf("Initial weight vector for matrix %d has incorrect length.", i))
      }
      norm_w <- sqrt(sum(weights[[i]]^2))
      if (abs(norm_w - 1) > 1e-10) {
        weights[[i]] <- weights[[i]] / norm_w
      }
    }
  }
  
  # Pre-compute cross-products to avoid repeated calculations
  cross_products <- vector("list", num_matrices * num_matrices)
  dim(cross_products) <- c(num_matrices, num_matrices)
  
  for (i in 1:num_matrices) {
    for (j in 1:num_matrices) {
      if (i != j) {
        cross_products[[i, j]] <- t(matrices[[i]]) %*% matrices[[j]]
      }
    }
  }
  
  # Main iteration loop
  k <- 0
  diff <- Inf
  prev_weights <- weights
  latent_vars <- vector("list", num_matrices)
  
  while ((diff > convergence_threshold) && (k < max_iterations)) {
    k <- k + 1
    
    # Compute latent variables
    for (i in 1:num_matrices) {
      latent_vars[[i]] <- matrices[[i]] %*% weights[[i]]
    }
    
    # Update weights
    for (i in 1:num_matrices) {
      # Compute gradient based on cross-covariance with other blocks
      new_weight <- numeric(num_weights[i])
      
      for (j in 1:num_matrices) {
        if (i != j) {
          new_weight <- new_weight + g_factors[i, j] * (cross_products[[i, j]] %*% weights[[j]])
        }
      }
      
      # Save the original weight before thresholding (for ensure_nonzero_blocks)
      original_weight <- new_weight
      
      # Apply soft-thresholding if sparsity is requested
      if (!no_sparse_matrix[i]) {
        # Soft thresholding
        new_weight <- soft_threshold(new_weight, sparsity[i])
      }
      
      # If ensure_nonzero_blocks is TRUE and all weights became zero or nearly zero,
      # keep the feature with the largest absolute original weight
      if (ensure_nonzero_blocks && all(abs(new_weight) < 1e-8)) {
        # Find the index of the feature with the highest absolute coefficient
        max_idx <- which.max(abs(original_weight))
        
        # Keep only the largest feature
        if (print_level >= 0) {
          message(sprintf("Block %d has all zero weights after thresholding. Ensuring at least one non-zero feature (feature %d)", i, max_idx))
        }
        
        new_weight <- rep(0, length(new_weight))
        new_weight[max_idx] <- sign(original_weight[max_idx]) * 0.1  # Small but non-zero value
      }
      
      # Normalize the weight vector
      norm_w <- sqrt(sum(new_weight^2))
      if (norm_w > .Machine$double.eps) {
        weights[[i]] <- new_weight / norm_w
      } else {
        if (print_level >= 0) {
          warning(sprintf("Weight vector %d is all zeros after thresholding. Using previous weights.", i))
        }
        weights[[i]] <- prev_weights[[i]]
      }
    }
    
    # Check convergence: compute difference from previous weights
    diff_weights <- lapply(1:num_matrices, function(i) {
      abs(sum(weights[[i]] * prev_weights[[i]])) - 1
    })
    diff <- max(abs(unlist(diff_weights)))
    
    # Store current weights for next iteration
    prev_weights <- weights
  }
  
  # Final computation of latent variables
  for (i in 1:num_matrices) {
    latent_vars[[i]] <- matrices[[i]] %*% weights[[i]]
  }
  
  # Compute final covariances
  for (i in 1:num_matrices) {
    for (j in 1:num_matrices) {
      covariances[[i, j]] <- t(latent_vars[[i]]) %*% latent_vars[[j]] / (nrow(matrices[[i]]) - 1)
    }
  }
  
  # Check convergence success
  converged <- diff <= convergence_threshold
  
  if (!converged && print_level >= 0) {
    warning(sprintf("Algorithm did not converge after %d iterations. Final diff: %g", 
                   k, diff))
  } else if (print_level >= 1) {
    message(sprintf("Algorithm converged after %d iterations. Final diff: %g", 
                   k, diff))
  }
  
  # Return results
  return(list(
    weights = weights,
    covariances = covariances,
    latent_vars = latent_vars,
    converged = converged,
    iterations = k,
    convergence_diff = diff,
    dimensions = dims
  ))
}

#' Soft thresholding function
#'
#' Applies soft thresholding to a vector based on the specified lambda.
#'
#' @param w Numeric vector to apply soft thresholding to
#' @param lambda Threshold value
#'
#' @return Soft-thresholded vector
#' @keywords internal
soft_threshold <- function(w, lambda) {
  # First, ensure the vector is properly formatted for our C++ function
  w <- as.vector(w)
  
  # Use the pure R implementation (C++ version will be used if available at build time)
  # Sort values by absolute magnitude
  abs_w <- abs(w)
  ord <- order(abs_w, decreasing = TRUE)
  abs_sorted <- abs_w[ord]
  
  # Compute the threshold
  cumsum_abs <- cumsum(abs_sorted)
  k <- which(abs_sorted > (cumsum_abs - lambda) / seq_along(w))[1]
  
  if (is.na(k)) {
    # All weights below threshold
    return(numeric(length(w)))
  }
  
  # Compute the threshold value
  threshold <- (cumsum_abs[k] - lambda) / k
  
  # Apply soft-thresholding
  result <- pmax(abs_w - threshold, 0) * sign(w)
  
  return(result)
}
