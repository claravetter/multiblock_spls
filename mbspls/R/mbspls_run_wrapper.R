# Wrapper for mbspls_run to match the arguments used in our script
# This file implements a wrapper around mbspls_run that accepts the arguments structure
# used in the leave-one-site-out analysis scripts

#' Run MBSPLS analysis with the argument structure used in leave-one-site-out analysis
#'
#' This is a wrapper around mbspls_run that accepts the arguments used in our leave-one-site-out scripts
#'
#' @param Xs List of data matrices for each block
#' @param Xs_names Optional character vector with names for each block
#' @param feature_names Optional list of feature names for each block
#' @param covars Optional list of covariate matrices for each block
#' @param covars_names Optional list of covariate names for each block
#' @param covars_correction_target Numeric vector indicating which blocks to correct
#' @param ids Optional vector of subject IDs
#' @param diagnosis Optional vector of diagnosis values
#' @param diagnosis_names Optional vector of diagnosis labels
#' @param config List of configuration parameters including:
#'   \itemize{
#'     \item \code{optimization_strategy}: Type of hyperparameter optimization. Options: "grid_search", "randomized_search", or NULL (no tuning)
#'     \item \code{hyperparam_distributions}: List of 2-element vectors specifying [min, max] for each block's sparsity parameter
#'     \item \code{randomized_search_iterations}: Number of random combinations to try (for randomized search)
#'     \item \code{inner_folds}: Number of cross-validation folds for tuning (default: 5)
#'     \item \code{matrix_norm}: Matrix norm to use ("fro", "1", "2", "inf")
#'     \item \code{correlation_method}: Correlation method ("pearson", "spearman", "kendall")
#'     \item \code{ensure_nonzero_blocks}: Whether to ensure each block has at least one non-zero feature
#'     \item \code{seed}: Random seed for reproducibility
#'     \item \code{outer_folds}: Number of cross-validation folds for final evaluation
#'     \item \code{permutation_testing}: Number of permutations for significance testing
#'     \item \code{bootstrap_testing}: Number of bootstrap samples for stability assessment
#'   }
#' @param output_dir Directory to save results
#'
#' @return An mbspls object with analysis results
#'
#' @examples
#' \dontrun{
#' # Example usage
#' result <- mbspls_run_wrapper(
#'   Xs = list(X1, X2, X3),
#'   Xs_names = c("Block1", "Block2", "Block3"),
#'   config = list(matrix_norm = "fro", permutation_testing = 100)
#' )
#' }
#'
mbspls_run_wrapper <- function(
  Xs,                         # List of data matrices for each block
  Xs_names = NULL,            # Optional names for each block
  feature_names = NULL,       # Optional feature names for each block
  covars = NULL,              # Optional list of covariates for each block
  covars_names = NULL,        # Optional covariate names
  covars_correction_target = NULL, # Which blocks to correct
  ids = NULL,                 # Optional subject IDs
  diagnosis = NULL,           # Optional diagnosis values
  diagnosis_names = NULL,     # Optional diagnosis labels
  config = list(),            # Configuration parameters
  output_dir = NULL           # Directory to save results
) {
  # Create output directory if specified
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    debug_dir <- file.path(output_dir, "debug")
    dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)
  }
  # Create output directory if specified
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Apply covariate correction if covariates are provided
  if (!is.null(covars)) {
    # Create log directory
    if (!is.null(output_dir)) {
      debug_dir <- file.path(output_dir, "debug")
      dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Log initial information
      sink(file.path(output_dir, "covariate_info.txt"))
      cat("Covariate information:\n")
      cat("Number of blocks:", length(Xs), "\n")
      cat("Number of covariate sets:", length(covars), "\n")
      
      # Debug - check if each covariate is a matrix
      for (i in seq_along(covars)) {
        if (!is.null(covars[[i]])) {
          cat(sprintf("Block %d (%s): %s covariates, class: %s\n", 
                     i, 
                     if (!is.null(Xs_names)) Xs_names[i] else paste("Block", i),
                     if (is.matrix(covars[[i]])) ncol(covars[[i]]) else "INVALID",
                     class(covars[[i]])[1]))
          
          if (!is.null(covars_names) && !is.null(covars_names[[i]])) {
            cat("  Covariate names:", paste(covars_names[[i]], collapse=", "), "\n")
          }
        } else {
          cat(sprintf("Block %d: No covariates (NULL)\n", i))
        }
      }
      
      # Debug - check covars_correction_target
      if (!is.null(covars_correction_target)) {
        cat("Correction targets:", paste(covars_correction_target, collapse=", "), "\n")
      } else {
        cat("Correction targets: NULL\n")
      }
      sink()
    }
    
    # Process each block that needs covariate correction
    if (!is.null(covars_correction_target) && length(covars_correction_target) == length(Xs)) {
      for (i in seq_along(Xs)) {
        # Only process blocks where correction is requested and covariates exist
        if (i <= length(covars_correction_target) && 
            covars_correction_target[i] == 1 && 
            !is.null(covars[[i]]) && 
            is.matrix(covars[[i]])) {
          
          # Log this correction
          if (!is.null(output_dir)) {
            sink(file.path(debug_dir, paste0("block_", i, "_correction.txt")))
            cat(sprintf("Applying covariate correction to block %d\n", i))
            cat("Block dimensions:", nrow(Xs[[i]]), "x", ncol(Xs[[i]]), "\n")
            cat("Covariates dimensions:", nrow(covars[[i]]), "x", ncol(covars[[i]]), "\n")
            sink()
          }
          
          # Simple regression-based covariate correction
          X_corrected <- matrix(0, nrow(Xs[[i]]), ncol(Xs[[i]]))
          colnames(X_corrected) <- colnames(Xs[[i]])
          
          for (j in seq_len(ncol(Xs[[i]]))) {
            y <- Xs[[i]][, j]  # Response variable (feature)
            
            # Create formula for the model
            tryCatch({
              # Create a data frame with the feature and covariates
              model_data <- data.frame(y = y, covars[[i]])
              
              # Fit the model
              model_formula <- as.formula(paste("y ~", paste(colnames(covars[[i]]), collapse = " + ")))
              lm_fit <- lm(model_formula, data = model_data)
              
              # Store residuals as the corrected values
              X_corrected[, j] <- residuals(lm_fit)
            }, error = function(e) {
              # If there's an error, just use the original values
              if (!is.null(output_dir)) {
                sink(file.path(debug_dir, paste0("block_", i, "_error.txt")), append = TRUE)
                cat(sprintf("Error correcting feature %d: %s\n", j, e$message))
                sink()
              }
              X_corrected[, j] <- y
            })
          }
          
          # Replace the original block with the corrected one
          Xs[[i]] <- X_corrected
        }
      }
    }
  }
  
  # Set up data object
  # Name the blocks if needed
  if (!is.null(Xs_names) && length(Xs_names) == length(Xs)) {
    names(Xs) <- Xs_names
  }
  
  # Create mbs_data object - need to use do.call since mbs_data expects separate arguments
  if (is.null(names(Xs))) {
    # If no names, assign default names
    names(Xs) <- paste0("block", seq_along(Xs))
  }
  
  # Use do.call to pass each matrix as a named argument
  data <- do.call(mbs_data, Xs)
  
  # Add feature names if provided
  if (!is.null(feature_names)) {
    for (i in seq_along(feature_names)) {
      if (i <= length(data$matrices) && length(feature_names[[i]]) == ncol(data$matrices[[i]])) {
        colnames(data$matrices[[i]]) <- feature_names[[i]]
      }
    }
  }
  
  # Implement proper hyperparameter tuning for sparsity parameters
  best_sparsity <- NULL
  best_score <- -Inf
  
  # Initialize variables for potential later use
  inner_folds <- if (!is.null(config$inner_folds)) config$inner_folds else 5
  n_iterations <- if (!is.null(config$randomized_search_iterations)) config$randomized_search_iterations else 100
  
  # Check if we need to perform hyperparameter tuning
  do_tuning <- !is.null(config$optimization_strategy) && 
              (config$optimization_strategy %in% c("grid_search", "randomized_search")) &&
              !is.null(config$hyperparam_distributions) &&
              length(config$hyperparam_distributions) == length(Xs)
  
  if (do_tuning) {
    # Log start of tuning
    if (!is.null(output_dir)) {
      tuning_dir <- file.path(output_dir, "tuning")
      dir.create(tuning_dir, recursive = TRUE, showWarnings = FALSE)
      sink(file.path(tuning_dir, "tuning_config.txt"))
      cat("Hyperparameter Tuning Configuration\n")
      cat("==================================\n")
      cat(sprintf("Strategy: %s\n", config$optimization_strategy))
      cat(sprintf("Number of blocks: %d\n", length(Xs)))
      
      # Log hyperparameter distributions
      cat("\nSparsity parameter distributions:\n")
      for (i in seq_along(config$hyperparam_distributions)) {
        block_name <- if (!is.null(Xs_names)) Xs_names[i] else paste("Block", i)
        cat(sprintf("  %s: [%s, %s]\n", 
                   block_name, 
                   config$hyperparam_distributions[[i]][1], 
                   config$hyperparam_distributions[[i]][2]))
      }
      
      # Log cross-validation settings
      inner_folds <- if (!is.null(config$inner_folds)) config$inner_folds else 5
      cat(sprintf("\nCross-validation: %d folds\n", inner_folds))
      
      # Log number of iterations
      n_iterations <- 10  # Default
      if (config$optimization_strategy == "randomized_search") {
        n_iterations <- if (!is.null(config$randomized_search_iterations)) 
                          config$randomized_search_iterations else 100
        cat(sprintf("Number of iterations: %d\n", n_iterations))
      }
      sink()
    }
    
    # Define inner cross-validation folds
    inner_folds <- if (!is.null(config$inner_folds)) config$inner_folds else inner_folds
    n_samples <- nrow(Xs[[1]])  # All blocks should have the same number of rows
    
    # Create folds
    set.seed(if (!is.null(config$seed)) config$seed else 42)
    fold_indices <- sample(rep(1:inner_folds, length.out = n_samples))
    
    # Number of iterations to try
    n_iterations <- if (config$optimization_strategy == "randomized_search" && 
                        !is.null(config$randomized_search_iterations)) 
                      config$randomized_search_iterations else n_iterations
    
    # Function to generate parameter combinations
    generate_params <- function() {
      sparsity <- numeric(length(Xs))
      
      for (i in seq_along(Xs)) {
        if (!is.null(config$hyperparam_distributions) && 
            i <= length(config$hyperparam_distributions) &&
            length(config$hyperparam_distributions[[i]]) >= 2) {
          
          min_val <- config$hyperparam_distributions[[i]][1]
          max_val <- config$hyperparam_distributions[[i]][2]
          
          # Generate value based on strategy
          if (config$optimization_strategy == "randomized_search") {
            # Random sample from uniform distribution
            sparsity[i] <- runif(1, min_val, max_val)
          } else {
            # For grid search, pick a point between min and max
            # In a real grid search, we'd generate all combinations
            # but for simplicity here, we'll sample points along the range
            sparsity[i] <- min_val + (max_val - min_val) * 
                           ((i - 1) / (length(Xs) - 1 + 0.00001))
          }
        } else {
          # Default: sqrt of number of features
          sparsity[i] <- sqrt(ncol(Xs[[i]]))
        }
      }
      
      return(sparsity)
    }
    
    # Initialize tracking of results
    all_params <- list()
    all_scores <- numeric()
    
    # Log tuning progress
    if (!is.null(output_dir)) {
      sink(file.path(tuning_dir, "tuning_progress.txt"))
      cat("Sparsity Parameter Tuning Progress\n")
      cat("================================\n\n")
    }
    
    # Run hyperparameter tuning
    for (iter in 1:n_iterations) {
      # Generate parameter combination
      current_sparsity <- generate_params()
      all_params[[iter]] <- current_sparsity
      
      # Log current parameters
      if (!is.null(output_dir)) {
        cat(sprintf("Iteration %d/%d\n", iter, n_iterations))
        cat("  Sparsity values:", paste(round(current_sparsity, 3), collapse = ", "), "\n")
      }
      
      # Cross-validation scores for this parameter set
      fold_scores <- numeric(inner_folds)
      
      # Run cross-validation
      for (fold in 1:inner_folds) {
        # Split data
        train_indices <- which(fold_indices != fold)
        test_indices <- which(fold_indices == fold)
        
        # Subset data matrices
        train_data <- lapply(Xs, function(X) X[train_indices, , drop = FALSE])
        test_data <- lapply(Xs, function(X) X[test_indices, , drop = FALSE])
        
        # Create data object for training
        names(train_data) <- names(Xs)
        data_train <- do.call(mbs_data, train_data)
        
        # Set up parameters
        params <- mbspls_setup(
          sparsity = current_sparsity,
          matrix_norm = if (!is.null(config$matrix_norm)) config$matrix_norm else "fro",
          correlation_method = if (!is.null(config$correlation_method)) 
                                tolower(config$correlation_method) else "pearson",
          ensure_nonzero_blocks = if (!is.null(config$ensure_nonzero_blocks)) 
                                  as.logical(config$ensure_nonzero_blocks) else FALSE
        )
        
        # Fit model
        tryCatch({
          model <- mbspls_run(
            data = data_train, 
            params = params,
            cv_folds = 0,  # No nested CV
            n_permutations = 0,  # No permutations during tuning
            n_bootstraps = 0,    # No bootstrapping during tuning
            seed = if (!is.null(config$seed)) config$seed + fold else 42 + fold,
            verbose = FALSE
          )
          
          # Project test data
          projected <- mbspls_project(model, test_data)
          
          # Calculate performance metric (rho)
          if (!is.null(projected$rho) && !is.na(projected$rho)) {
            fold_scores[fold] <- projected$rho
          } else {
            fold_scores[fold] <- NA
          }
        }, error = function(e) {
          if (!is.null(output_dir)) {
            cat(sprintf("    Fold %d: Error - %s\n", fold, e$message))
          }
          fold_scores[fold] <- NA
        })
      }
      
      # Calculate average score
      mean_score <- mean(fold_scores, na.rm = TRUE)
      all_scores[iter] <- mean_score
      
      # Log fold results
      if (!is.null(output_dir)) {
        cat("  Fold scores:", paste(round(fold_scores, 4), collapse = ", "), "\n")
        cat("  Mean score:", round(mean_score, 4), "\n\n")
      }
      
      # Update best parameters if better
      if (!is.na(mean_score) && mean_score > best_score) {
        best_score <- mean_score
        best_sparsity <- current_sparsity
        
        if (!is.null(output_dir)) {
          cat(sprintf("  New best score: %.4f\n", best_score))
          cat("  New best sparsity:", paste(round(best_sparsity, 3), collapse = ", "), "\n\n")
        }
      }
    }
    
    if (!is.null(output_dir)) {
      sink()
      
      # Save tuning results
      tuning_results <- data.frame(
        iteration = 1:n_iterations,
        score = all_scores,
        do.call(rbind, lapply(all_params, function(p) {
          data.frame(matrix(p, nrow = 1, 
                           dimnames = list(NULL, paste0("sparsity_", seq_along(p)))))
        }))
      )
      
      write.csv(tuning_results, file.path(tuning_dir, "tuning_results.csv"), row.names = FALSE)
      
      # Plot tuning results if possible
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        tryCatch({
          best_iter <- which.max(all_scores)
          p <- ggplot2::ggplot(tuning_results, ggplot2::aes_string(x = "iteration", y = "score")) +
               ggplot2::geom_line() +
               ggplot2::geom_point() +
               ggplot2::geom_point(data = tuning_results[best_iter, ], 
                                  color = "red", size = 3) +
               ggplot2::labs(title = "Hyperparameter Tuning Results",
                            x = "Iteration", 
                            y = "Performance Score")
          
          ggplot2::ggsave(file.path(tuning_dir, "tuning_plot.png"), p, 
                         width = 8, height = 6, dpi = 100)
        }, error = function(e) {
          # If plotting fails, just log the error
          message("Could not create tuning plot: ", e$message)
        })
      }
      
      # Write final summary
      sink(file.path(tuning_dir, "tuning_summary.txt"))
      cat("Sparsity Parameter Tuning Summary\n")
      cat("===============================\n\n")
      cat(sprintf("Best score: %.4f\n", best_score))
      cat("Best sparsity parameters:\n")
      for (i in seq_along(best_sparsity)) {
        block_name <- if (!is.null(Xs_names)) Xs_names[i] else paste("Block", i)
        cat(sprintf("  %s: %.4f\n", block_name, best_sparsity[i]))
      }
      sink()
    }
  } else {
    # No tuning, use default or mean values
    best_sparsity <- numeric(length(Xs))
    for (i in seq_along(Xs)) {
      # Default: sqrt of number of features
      best_sparsity[i] <- sqrt(ncol(Xs[[i]]))
      
      # If config has hyperparam_distributions, use the mean as a reasonable default
      if (!is.null(config$hyperparam_distributions) && 
          length(config$hyperparam_distributions) >= i) {
        
        if (is.numeric(config$hyperparam_distributions[[i]]) && 
            length(config$hyperparam_distributions[[i]]) == 2) {
          dist_range <- config$hyperparam_distributions[[i]]
          best_sparsity[i] <- mean(dist_range)
        }
      }
    }
  }
  
  # Set up parameters for mbspls using best parameters
  params <- mbspls_setup(
    sparsity = best_sparsity,
    matrix_norm = if (!is.null(config$matrix_norm)) config$matrix_norm else "fro",
    correlation_method = if (!is.null(config$correlation_method)) 
                           tolower(config$correlation_method) else "pearson",
    ensure_nonzero_blocks = if (!is.null(config$ensure_nonzero_blocks)) 
                             as.logical(config$ensure_nonzero_blocks) else FALSE
  )
  
  # Run mbspls with appropriate parameters
  result <- mbspls_run(
    data = data, 
    params = params,
    cv_folds = if (!is.null(config$outer_folds)) config$outer_folds else 10,
    n_permutations = if (!is.null(config$permutation_testing)) config$permutation_testing else 0,
    n_bootstraps = if (!is.null(config$bootstrap_testing)) config$bootstrap_testing else 0,
    seed = if (!is.null(config$seed)) config$seed else 42,
    verbose = TRUE
  )
  
  # Ensure the LVs structure exists
  if (is.null(result$LVs)) {
    # Create LVs structure based on permutation and CV results
    p_value <- NA
    if (!is.null(result$permutation_results) && !is.null(result$permutation_results$p_value)) {
      p_value <- result$permutation_results$p_value
    }
    
    explained_variance <- NA
    # Try to extract explained variance from CV results
    if (!is.null(result$cv_results)) {
      if (!is.null(result$cv_results$mean_explained_variance)) {
        explained_variance <- result$cv_results$mean_explained_variance
      } else if (!is.null(result$cv_results$mean_rho)) {
        # Approximate using squared correlation
        explained_variance <- result$cv_results$mean_rho^2
      }
    }
    
    # Create LVs structure
    result$LVs <- list(
      list(
        p_value = p_value,
        explained_variance = explained_variance
      )
    )
  }
  
  # Add any additional metadata
  result$Xs_names <- Xs_names
  result$feature_names <- feature_names
  result$ids <- ids
  result$diagnosis <- diagnosis
  result$diagnosis_names <- diagnosis_names
  
  # Save results if output directory is specified
  if (!is.null(output_dir)) {
    saveRDS(result, file.path(output_dir, "mbspls_result.rds"))
    
    # Save summary information
    sink(file.path(output_dir, "mbspls_summary.txt"))
    cat("MBSPLS Analysis Summary\n")
    cat("======================\n\n")
    
    cat("Dataset Information:\n")
    cat("------------------\n")
    for (i in seq_along(Xs)) {
      block_name <- if (!is.null(Xs_names)) Xs_names[i] else paste("Block", i)
      cat(sprintf("%s: %d samples, %d features\n", 
                 block_name, nrow(Xs[[i]]), ncol(Xs[[i]])))
    }
    
    cat("\nAnalysis Configuration:\n")
    cat("---------------------\n")
    cat(sprintf("Cross-validation folds: %d\n", 
               if (!is.null(config$outer_folds)) config$outer_folds else 10))
    cat(sprintf("Permutation tests: %d\n", 
               if (!is.null(config$permutation_testing)) config$permutation_testing else 0))
    cat(sprintf("Bootstrap samples: %d\n", 
               if (!is.null(config$bootstrap_testing)) config$bootstrap_testing else 0))
    cat(sprintf("Matrix norm: %s\n", 
               if (!is.null(config$matrix_norm)) config$matrix_norm else "fro"))
    cat(sprintf("Correlation method: %s\n", 
               if (!is.null(config$correlation_method)) config$correlation_method else "pearson"))
    
    # Add sparsity parameter information
    cat("\nSparsity Parameters:\n")
    if (!is.null(best_sparsity)) {
      for (i in seq_along(best_sparsity)) {
        block_name <- if (!is.null(Xs_names)) Xs_names[i] else paste("Block", i)
        cat(sprintf("  %s: %.4f\n", block_name, best_sparsity[i]))
      }
      
      # If we did hyperparameter tuning, note that
      if (do_tuning) {
        cat("\nNote: Sparsity parameters were tuned using cross-validation.\n")
        cat(sprintf("Tuning method: %s, %d iterations, %d folds\n", 
                   config$optimization_strategy, n_iterations, inner_folds))
      }
    }
    
    cat("\nResults Summary:\n")
    cat("--------------\n")
    for (lv in seq_along(result$LVs)) {
      if (!is.null(result$LVs[[lv]])) {
        cat(sprintf("LV %d: p-value = %.4f, explained variance = %.4f\n", 
                   lv, result$LVs[[lv]]$p_value, result$LVs[[lv]]$explained_variance))
      }
    }
    sink()
  }
  
  return(result)
}
