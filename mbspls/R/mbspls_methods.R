#' Print method for mbspls objects
#'
#' @param x An mbspls object
#' @param ... Further arguments passed to or from other methods
#'
#' @return Invisibly returns the input object
#' @export
print.mbspls <- function(x, ...) {
  cat("Multi-Block Sparse PLS Results\n")
  cat("------------------------------\n")
  
  # Basic information
  cat("Number of blocks:", length(x$weights), "\n")
  cat("Correlation/norm:", sprintf("%.4f", x$rho), "\n")
  
  # CV information if available
  if (!is.null(x$cv_results)) {
    cat("\nCross-validation results:\n")
    cat("  Number of folds:", length(x$cv_results$fold_rhos), "\n")
    cat("  Mean correlation:", sprintf("%.4f", x$cv_results$mean_rho), "\n")
    cat("  SD correlation:", sprintf("%.4f", x$cv_results$sd_rho), "\n")
  }
  
  # Permutation test information
  if (!is.null(x$permutation_results)) {
    cat("\nPermutation test results:\n")
    cat("  Number of permutations:", length(x$permutation_results$perm_rhos), "\n")
    cat("  p-value:", sprintf("%.4f", x$permutation_results$p_value), "\n")
    cat("  Mean perm. correlation:", sprintf("%.4f", x$permutation_results$mean_perm_rho), "\n")
  }
  
  # Bootstrap information
  if (!is.null(x$bootstrap_results)) {
    cat("\nBootstrap results:\n")
    cat("  Number of bootstrap samples:", length(x$bootstrap_results$boot_rhos), "\n")
    cat("  Mean bootstrap correlation:", sprintf("%.4f", x$bootstrap_results$mean_boot_rho), "\n")
    cat("  SD bootstrap correlation:", sprintf("%.4f", x$bootstrap_results$sd_boot_rho), "\n")
  }
  
  invisible(x)
}

#' Summary method for mbspls objects
#'
#' @param object An mbspls object
#' @param ... Further arguments passed to or from other methods
#'
#' @return A list with summary statistics
#' @export
summary.mbspls <- function(object, ...) {
  # Count non-zero loadings in each block
  nonzero_counts <- lapply(object$weights, function(w) sum(abs(w) > 1e-10))
  names(nonzero_counts) <- names(object$weights)
  
  # Calculate proportion of variance explained by each block
  if (!is.null(object$latent_vars)) {
    if (is.matrix(object$latent_vars)) {
      # Calculate R-squared between latent variables
      R2 <- cor(object$latent_vars)^2
    } else if (is.list(object$latent_vars)) {
      # Extract latent variables from list form
      lv_mat <- do.call(cbind, object$latent_vars)
      R2 <- cor(lv_mat)^2
    } else {
      R2 <- NA
    }
  } else {
    R2 <- NA
  }
  
  # Compile summary statistics
  result <- list(
    rho = object$rho,
    nonzero_loadings = nonzero_counts,
    R_squared = R2
  )
  
  # Add CV results if available
  if (!is.null(object$cv_results)) {
    result$cv <- list(
      mean_rho = object$cv_results$mean_rho,
      sd_rho = object$cv_results$sd_rho,
      q2 = 1 - (1 - object$cv_results$mean_rho^2)  # Simple Q2 approximation
    )
  }
  
  # Add permutation test results if available
  if (!is.null(object$permutation_results)) {
    result$permutation <- list(
      p_value = object$permutation_results$p_value,
      mean_rho = object$permutation_results$mean_perm_rho,
      sd_rho = object$permutation_results$sd_perm_rho
    )
  }
  
  # Add bootstrap results if available
  if (!is.null(object$bootstrap_results)) {
    # Calculate consistency of feature selection
    if (length(object$bootstrap_results$boot_weights) > 0) {
      feature_consistency <- lapply(seq_along(object$weights), function(i) {
        orig_nonzero <- abs(object$weights[[i]]) > 1e-10
        boot_weights <- object$bootstrap_results$boot_weights[[i]]
        
        # For each bootstrap, check which features are selected
        boot_nonzero <- apply(boot_weights, 2, function(w) abs(w) > 1e-10)
        
        # Calculate selection frequency for each feature
        selection_freq <- rowMeans(boot_nonzero)
        
        # Calculate consistency - proportion of bootstraps that select the same features
        consistency <- mean(orig_nonzero == (selection_freq > 0.5))
        
        # Return selection frequencies and consistency
        list(
          selection_freq = selection_freq,
          consistency = consistency
        )
      })
      names(feature_consistency) <- names(object$weights)
      
      result$bootstrap <- list(
        mean_rho = object$bootstrap_results$mean_boot_rho,
        sd_rho = object$bootstrap_results$sd_boot_rho,
        feature_consistency = feature_consistency
      )
    }
  }
  
  class(result) <- "summary.mbspls"
  return(result)
}

#' Print method for summary.mbspls objects
#'
#' @param x A summary.mbspls object
#' @param ... Further arguments passed to or from other methods
#'
#' @return Invisibly returns the input object
#' @export
print.summary.mbspls <- function(x, ...) {
  cat("Summary of Multi-Block Sparse PLS Results\n")
  cat("----------------------------------------\n\n")
  
  cat("Correlation/norm:", sprintf("%.4f", x$rho), "\n\n")
  
  cat("Non-zero loadings per block:\n")
  for (i in seq_along(x$nonzero_loadings)) {
    block_name <- names(x$nonzero_loadings)[i]
    count <- x$nonzero_loadings[[i]]
    cat(sprintf("  %s: %d\n", block_name, count))
  }
  
  cat("\nR-squared between blocks:\n")
  if (!is.list(x$R_squared) && !is.matrix(x$R_squared)) {
    cat("  Not available\n")
  } else if (is.matrix(x$R_squared)) {
    # Pretty print the matrix
    print(round(x$R_squared, 4))
  }
  
  if (!is.null(x$cv)) {
    cat("\nCross-validation:\n")
    cat("  Mean correlation:", sprintf("%.4f", x$cv$mean_rho), "\n")
    cat("  SD correlation:", sprintf("%.4f", x$cv$sd_rho), "\n")
    cat("  Q² approximation:", sprintf("%.4f", x$cv$q2), "\n")
  }
  
  if (!is.null(x$permutation)) {
    cat("\nPermutation test:\n")
    cat("  p-value:", sprintf("%.4f", x$permutation$p_value), "\n")
    cat("  Mean perm. correlation:", sprintf("%.4f", x$permutation$mean_rho), "\n")
    cat("  SD perm. correlation:", sprintf("%.4f", x$permutation$sd_rho), "\n")
  }
  
  if (!is.null(x$bootstrap)) {
    cat("\nBootstrap results:\n")
    cat("  Mean bootstrap correlation:", sprintf("%.4f", x$bootstrap$mean_rho), "\n")
    cat("  SD bootstrap correlation:", sprintf("%.4f", x$bootstrap$sd_rho), "\n")
    cat("\nFeature selection consistency:\n")
    for (i in seq_along(x$bootstrap$feature_consistency)) {
      block_name <- names(x$bootstrap$feature_consistency)[i]
      consistency <- x$bootstrap$feature_consistency[[i]]$consistency
      cat(sprintf("  %s: %.2f%%\n", block_name, consistency * 100))
    }
  }
  
  invisible(x)
}

#' Extract model coefficients from an mbspls object
#'
#' @param object An mbspls object
#' @param ... Further arguments (not used)
#'
#' @return A list with the weight vectors (coefficients) for each block
#' @export
coef.mbspls <- function(object, ...) {
  # Return the weight vectors
  return(object$weights)
}

#' Predict method for mbspls objects
#'
#' Project new data onto the MBSPLS latent space.
#'
#' @param object An mbspls object
#' @param newdata A list of matrices (new data blocks)
#' @param ... Further arguments (not used)
#'
#' @return A list containing projected latent variables and correlations
#' @export
predict.mbspls <- function(object, newdata, ...) {
  # Check input data
  if (!is.list(newdata) || !all(sapply(newdata, is.matrix))) {
    stop("newdata must be a list of matrices")
  }
  
  if (length(newdata) != length(object$weights)) {
    stop("Number of matrices in newdata must match the number of blocks in the model")
  }
  
  # Project new data onto the latent space
  projection <- get("mbspls_project", envir = asNamespace("mbspls"))(
    data = newdata,
    weights = object$weights,
    correlation_method = object$params$correlation_method,
    matrix_norm = object$params$matrix_norm
  )
  
  # Return projection results
  return(list(
    latent_vars = projection$latent_vars,
    rho = projection$rho
  ))
}

#' Plot method for mbspls objects
#'
#' Creates plots for MBSPLS results, including latent variables, loadings, and cross-validation results.
#'
#' @param x An mbspls object
#' @param what What to plot: "latent_vars", "loadings", "cv", "permutation", or "bootstrap"
#' @param blocks Which blocks to include in the plot (for latent_vars and loadings)
#' @param ... Further arguments passed to the plotting functions
#'
#' @return A ggplot object (if ggplot2 is available) or NULL (with base graphics)
#' @export
plot.mbspls <- function(x, what = c("latent_vars", "loadings", "cv", "permutation", "bootstrap"),
                        blocks = NULL, ...) {
  what <- match.arg(what)
  
  # Check if ggplot2 is available
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  
  if (what == "latent_vars") {
    # Plot latent variables
    if (is.null(blocks) || length(blocks) != 2) {
      # Default to first two blocks if not specified
      blocks <- seq_len(min(2, length(x$latent_vars)))
    }
    
    # Extract latent variables
    if (is.matrix(x$latent_vars)) {
      lv1 <- x$latent_vars[, blocks[1]]
      lv2 <- x$latent_vars[, blocks[2]]
      block_names <- colnames(x$latent_vars)[blocks]
      if (is.null(block_names)) {
        block_names <- paste("Block", blocks)
      }
    } else if (is.list(x$latent_vars)) {
      lv1 <- x$latent_vars[[blocks[1]]]
      lv2 <- x$latent_vars[[blocks[2]]]
      block_names <- names(x$latent_vars)[blocks]
      if (is.null(block_names)) {
        block_names <- paste("Block", blocks)
      }
    } else {
      stop("Latent variables not available in the expected format")
    }
    
    # Create plot
    if (has_ggplot) {
      df <- data.frame(LV1 = lv1, LV2 = lv2)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = LV1, y = LV2)) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue") +
        ggplot2::labs(
          x = paste("Latent Variable:", block_names[1]),
          y = paste("Latent Variable:", block_names[2]),
          title = paste("MBSPLS Latent Variables"),
          subtitle = paste("Correlation:", round(x$rho, 4))
        ) +
        ggplot2::theme_minimal()
      return(p)
    } else {
      # Base R plot
      plot(lv1, lv2, main = "MBSPLS Latent Variables",
           xlab = paste("Latent Variable:", block_names[1]),
           ylab = paste("Latent Variable:", block_names[2]),
           pch = 16)
      abline(lm(lv2 ~ lv1), col = "blue")
      text(max(lv1) * 0.8, max(lv2) * 0.9, 
           paste("Correlation:", round(x$rho, 4)))
      return(invisible(NULL))
    }
  } else if (what == "loadings") {
    # Plot loadings
    if (is.null(blocks)) {
      blocks <- seq_along(x$weights)
    }
    
    if (has_ggplot) {
      # Create a data frame for plotting
      df_list <- lapply(blocks, function(i) {
        block_name <- names(x$weights)[i]
        if (is.null(block_name)) block_name <- paste("Block", i)
        
        w <- x$weights[[i]]
        names_w <- names(w)
        if (is.null(names_w)) names_w <- paste0("Var", seq_along(w))
        
        data.frame(
          Variable = names_w,
          Loading = w,
          Block = block_name,
          stringsAsFactors = FALSE
        )
      })
      
      df <- do.call(rbind, df_list)
      
      # Create plot
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Loading, fill = Block)) +
        ggplot2::geom_col() +
        ggplot2::facet_wrap(~ Block, scales = "free_x") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(
          title = "MBSPLS Loadings",
          x = "Variable",
          y = "Loading Weight"
        )
      
      return(p)
    } else {
      # Base R plot - split screen for multiple blocks
      n_blocks <- length(blocks)
      par(mfrow = c(ceiling(n_blocks/2), min(n_blocks, 2)))
      
      for (i in blocks) {
        block_name <- names(x$weights)[i]
        if (is.null(block_name)) block_name <- paste("Block", i)
        
        w <- x$weights[[i]]
        names_w <- names(w)
        if (is.null(names_w)) names_w <- paste0("Var", seq_along(w))
        
        barplot(w, main = paste("Loadings -", block_name),
                ylab = "Loading Weight", las = 2)
        abline(h = 0, lty = 2)
      }
      
      return(invisible(NULL))
    }
  } else if (what == "cv" && !is.null(x$cv_results)) {
    # Plot cross-validation results
    if (has_ggplot) {
      df <- data.frame(
        Fold = seq_along(x$cv_results$fold_rhos),
        Correlation = x$cv_results$fold_rhos
      )
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Fold, y = Correlation)) +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = x$cv_results$mean_rho, color = "red", linetype = "dashed") +
        ggplot2::labs(
          title = "Cross-validation Results",
          subtitle = paste("Mean correlation:", round(x$cv_results$mean_rho, 4)),
          x = "Fold",
          y = "Correlation"
        ) +
        ggplot2::theme_minimal()
      
      return(p)
    } else {
      # Base R plot
      plot(x$cv_results$fold_rhos, type = "b", 
           main = "Cross-validation Results",
           xlab = "Fold", ylab = "Correlation")
      abline(h = x$cv_results$mean_rho, col = "red", lty = 2)
      text(length(x$cv_results$fold_rhos) * 0.8, x$cv_results$mean_rho * 1.1,
           paste("Mean:", round(x$cv_results$mean_rho, 4)))
      return(invisible(NULL))
    }
  } else if (what == "permutation" && !is.null(x$permutation_results)) {
    # Plot permutation test results
    if (has_ggplot) {
      df <- data.frame(
        Permutation = seq_along(x$permutation_results$perm_rhos),
        Correlation = x$permutation_results$perm_rhos
      )
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Correlation)) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::geom_vline(xintercept = x$rho, color = "red") +
        ggplot2::labs(
          title = "Permutation Test Results",
          subtitle = paste("p-value:", round(x$permutation_results$p_value, 4)),
          x = "Correlation",
          y = "Count"
        ) +
        ggplot2::theme_minimal()
      
      return(p)
    } else {
      # Base R plot
      hist(x$permutation_results$perm_rhos, 
           main = "Permutation Test Results",
           xlab = "Correlation", ylab = "Count")
      abline(v = x$rho, col = "red")
      text(x$rho * 1.1, length(x$permutation_results$perm_rhos) * 0.1,
           paste("p-value:", round(x$permutation_results$p_value, 4)))
      return(invisible(NULL))
    }
  } else if (what == "bootstrap" && !is.null(x$bootstrap_results)) {
    # Plot bootstrap results
    if (is.null(blocks)) {
      blocks <- 1
    }
    block_idx <- blocks[1]  # Just use the first requested block
    
    if (has_ggplot && length(x$bootstrap_results$boot_weights) > 0) {
      # Extract bootstrap weights for the selected block
      boot_weights <- x$bootstrap_results$boot_weights[[block_idx]]
      
      # Calculate mean and standard deviation for each variable
      means <- rowMeans(boot_weights)
      sds <- apply(boot_weights, 1, sd)
      
      # Create data frame for plotting
      df <- data.frame(
        Variable = seq_along(means),
        Mean = means,
        Lower = means - sds,
        Upper = means + sds
      )
      
      # Add original weights
      df$Original <- x$weights[[block_idx]]
      
      # Create plot
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Variable)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), 
                              width = 0.2, alpha = 0.5) +
        ggplot2::geom_point(ggplot2::aes(y = Mean), color = "blue") +
        ggplot2::geom_point(ggplot2::aes(y = Original), color = "red", shape = 4) +
        ggplot2::labs(
          title = paste("Bootstrap Results -", names(x$weights)[block_idx]),
          subtitle = "Blue dots: bootstrap mean; Red crosses: original weights",
          x = "Variable",
          y = "Weight"
        ) +
        ggplot2::theme_minimal()
      
      return(p)
    } else {
      # Base R plot
      if (length(x$bootstrap_results$boot_weights) > 0) {
        boot_weights <- x$bootstrap_results$boot_weights[[block_idx]]
        
        means <- rowMeans(boot_weights)
        sds <- apply(boot_weights, 1, sd)
        
        # Plot means with error bars
        plot(seq_along(means), means, 
             main = paste("Bootstrap Results -", names(x$weights)[block_idx]),
             xlab = "Variable", ylab = "Weight",
             pch = 16, col = "blue")
        
        # Add error bars
        for (i in seq_along(means)) {
          lines(c(i, i), c(means[i] - sds[i], means[i] + sds[i]), col = "blue")
        }
        
        # Add original weights
        points(seq_along(x$weights[[block_idx]]), x$weights[[block_idx]], 
               col = "red", pch = 4)
        
        legend("topright", 
               legend = c("Bootstrap mean", "Original weights"),
               pch = c(16, 4),
               col = c("blue", "red"))
      } else {
        plot(x$bootstrap_results$boot_rhos, 
             main = "Bootstrap Results",
             xlab = "Bootstrap Sample", ylab = "Correlation")
        abline(h = mean(x$bootstrap_results$boot_rhos), col = "red", lty = 2)
      }
      
      return(invisible(NULL))
    }
  } else {
    stop("Invalid 'what' argument or requested data not available")
  }
}
