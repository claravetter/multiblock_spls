#' Run MBSPLS analysis
#'
#' This is the main function for Multi-Block Sparse PLS analysis. It handles
#' cross-validation, permutation testing, and bootstrapping.
#'
#' @param data A list of matrices or an mbs_data object
#' @param params An mbspls_params object from mbspls_setup()
#' @param cv_folds Number of cross-validation folds. Default: 10
#' @param n_permutations Number of permutations for significance testing. Default: 0
#' @param n_bootstraps Number of bootstrap samples for stability assessment. Default: 0
#' @param seed Random seed for reproducibility. Default: 42
#' @param verbose Logical; if TRUE, prints progress. Default: TRUE
#'
#' @return An mbspls object with the following components:
#'   \item{weights}{Final weight vectors for each block}
#'   \item{latent_vars}{Latent variables for each block}
#'   \item{rho}{Correlation values or matrix norms}
#'   \item{cv_results}{Cross-validation results}
#'   \item{permutation_results}{Permutation test results (if requested)}
#'   \item{bootstrap_results}{Bootstrap results (if requested)}
#'   \item{params}{Parameters used for the analysis}
#'
#' @importFrom stats cor sd
#' @importFrom utils packageVersion
#' @useDynLib mbspls, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#'
#' @examples
#' \dontrun{
#' # Create synthetic multi-block data
#' X1 <- matrix(rnorm(100*20), 100, 20)
#' X2 <- matrix(rnorm(100*30), 100, 30)
#' data <- mbs_data(block1 = X1, block2 = X2)
#'
#' # Set up parameters
#' params <- mbspls_setup(sparsity = c(2, 2))
#'
#' # Run MBSPLS with cross-validation
#' result <- mbspls_run(data, params, cv_folds = 5)
#' }
mbspls_run <- function(data, params, cv_folds = 10, 
                      n_permutations = 0, n_bootstraps = 0,
                      seed = 42, verbose = TRUE) {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Check and convert input data if needed
  if (!inherits(data, "mbs_data")) {
    if (is.list(data) && all(sapply(data, is.matrix))) {
      data <- get("mbs_data", envir = asNamespace("mbspls"))(matrices = data)
    } else {
      stop("Data must be an mbs_data object or a list of matrices")
    }
  }
  
  # Extract matrices for easier handling
  matrices <- data$matrices
  
  # Check that parameters match data
  if (length(params$sparsity) != length(matrices)) {
    stop("Number of sparsity parameters must match number of data blocks")
  }
  
  # Set up parallelization if needed
  if (n_permutations > 0 || n_bootstraps > 0 || cv_folds > 1) {
    if (requireNamespace("future", quietly = TRUE) && 
        requireNamespace("furrr", quietly = TRUE)) {
      
      if (verbose) {
        message("Setting up parallel processing with ", 
                params$parallel$workers, " workers")
      }
      
      # Save current plan and restore it when done
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      
      # Set up parallel plan according to parameters
      if (params$parallel$strategy == "multisession") {
        future::plan(future::multisession, workers = params$parallel$workers)
      } else if (params$parallel$strategy == "multicore") {
        future::plan(future::multicore, workers = params$parallel$workers)
      } else if (params$parallel$strategy == "cluster") {
        # Assumes a pre-configured cluster
        if (!exists("cl", envir = globalenv())) {
          warning("Cluster 'cl' not found in global environment, falling back to multisession")
          future::plan(future::multisession, workers = params$parallel$workers)
        } else {
          future::plan(future::cluster, workers = get("cl", envir = globalenv()))
        }
      } else {
        # Default to sequential for unknown strategies
        future::plan(future::sequential)
      }
    } else {
      warning("Packages 'future' and 'furrr' are recommended for parallel processing")
    }
  }
  
  # Initialize results
  results <- list(
    weights = list(),
    latent_vars = list(),
    rho = NULL,
    cv_results = NULL,
    permutation_results = NULL,
    bootstrap_results = NULL,
    params = params,
    call = match.call(),
    version = utils::packageVersion("mbspls"),
    LVs = NULL  # Initialize LVs list that will store p-values and explained variance
  )
  
  # Run full dataset analysis
  if (verbose) message("Fitting MBSPLS on full dataset")
  
  # Use temporary placeholder as test data since we're using all data
  full_fit <- get("mbspls_fit_fold", envir = asNamespace("mbspls"))(
    training_data = matrices,
    test_data = matrices,
    params = params
  )
  
  results$weights <- full_fit$weights
  results$latent_vars <- full_fit$latent_vars
  results$rho <- full_fit$rho
  
  # Run cross-validation if requested
  if (cv_folds > 1) {
    if (verbose) message("Performing ", cv_folds, "-fold cross-validation")
    
    # Create fold indices
    n_samples <- nrow(matrices[[1]])
    fold_indices <- split(sample(n_samples), rep(1:cv_folds, length.out = n_samples))
    
    # Function to run a single CV fold
    run_cv_fold <- function(fold_idx) {
      # Create training and test indices
      test_idx <- fold_indices[[fold_idx]]
      train_idx <- setdiff(1:n_samples, test_idx)
      
      # Split data into training and test sets
      train_data <- lapply(matrices, function(mat) mat[train_idx, , drop = FALSE])
      test_data <- lapply(matrices, function(mat) mat[test_idx, , drop = FALSE])
      
      # Run MBSPLS on this fold
      fold_fit <- get("mbspls_fit_fold", envir = asNamespace("mbspls"))(
        training_data = train_data,
        test_data = test_data,
        params = params
      )
      
      # Return fold results
      return(list(
        rho = fold_fit$rho,
        weights = fold_fit$weights,
        test_idx = test_idx
      ))
    }
    
    # Run CV folds (in parallel if possible)
    if (requireNamespace("furrr", quietly = TRUE) && 
        !inherits(future::plan(), "sequential")) {
      cv_results <- furrr::future_map(1:cv_folds, run_cv_fold, 
                                     .options = furrr::furrr_options(seed = params$parallel$seed))
    } else {
      cv_results <- lapply(1:cv_folds, run_cv_fold)
    }
    
    # Extract and summarize CV results
    cv_rhos <- sapply(cv_results, function(x) x$rho)
    
    results$cv_results <- list(
      folds = cv_results,
      mean_rho = mean(cv_rhos),
      sd_rho = sd(cv_rhos),
      fold_rhos = cv_rhos
    )
  }
  
  # Run permutation testing if requested
  if (n_permutations > 0) {
    if (verbose) message("Running ", n_permutations, " permutations")
    
    # Function to run a single permutation
    run_permutation <- function(perm_idx) {
      # Permute the samples in each block independently
      perm_matrices <- lapply(matrices, function(mat) {
        mat[sample(nrow(mat)), , drop = FALSE]
      })
      
      # Run MBSPLS on permuted data
      perm_fit <- get("mbspls_fit_fold", envir = asNamespace("mbspls"))(
        training_data = perm_matrices,
        test_data = perm_matrices,
        params = params
      )
      
      # Return permutation results
      return(list(
        rho = perm_fit$rho,
        converged = perm_fit$converged
      ))
    }
    
    # Run permutations (in parallel if possible)
    if (requireNamespace("furrr", quietly = TRUE) && 
        !inherits(future::plan(), "sequential")) {
      perm_results <- furrr::future_map(1:n_permutations, run_permutation,
                                       .options = furrr::furrr_options(seed = params$parallel$seed))
    } else {
      perm_results <- lapply(1:n_permutations, run_permutation)
    }
    
    # Extract permutation rhos and calculate p-value
    perm_rhos <- sapply(perm_results, function(x) x$rho)
    perm_converged <- sapply(perm_results, function(x) x$converged)
    
    # Calculate p-value (proportion of permutations with rho >= observed rho)
    # Handle NA values
    if (is.na(results$rho) || all(is.na(perm_rhos))) {
      p_value <- NA
    } else {
      valid_perms <- !is.na(perm_rhos)
      if (sum(valid_perms) > 0) {
        p_value <- mean(perm_rhos[valid_perms] >= results$rho, na.rm = TRUE)
      } else {
        p_value <- NA
      }
    }
    
    results$permutation_results <- list(
      perm_rhos = perm_rhos,
      converged = perm_converged,
      p_value = p_value,
      mean_perm_rho = mean(perm_rhos, na.rm = TRUE),
      sd_perm_rho = sd(perm_rhos, na.rm = TRUE)
    )
    
    # Calculate explained variance for CV results
    # For now, we'll use a simple measure: (rho^2) which represents the squared correlation
    # This is a simplified approximation of explained variance
    mean_explained_variance <- results$rho^2
    
    # Add to CV results
    if (!is.null(results$cv_results)) {
      results$cv_results$mean_explained_variance <- mean_explained_variance
    } else {
      results$cv_results <- list(
        mean_explained_variance = mean_explained_variance
      )
    }
    
    # Create LVs structure with p-value and explained variance
    results$LVs <- list(
      list(
        p_value = p_value,
        explained_variance = mean_explained_variance
      )
    )
  }
  
  # Run bootstrap if requested
  if (n_bootstraps > 0) {
    if (verbose) message("Running ", n_bootstraps, " bootstrap samples")
    
    # Function to run a single bootstrap
    run_bootstrap <- function(boot_idx) {
      # Create bootstrap sample (with replacement)
      n_samples <- nrow(matrices[[1]])
      boot_idx <- sample(n_samples, replace = TRUE)
      
      # Sample data with replacement
      boot_matrices <- lapply(matrices, function(mat) {
        mat[boot_idx, , drop = FALSE]
      })
      
      # Run MBSPLS on bootstrap sample
      boot_fit <- get("mbspls_fit_fold", envir = asNamespace("mbspls"))(
        training_data = boot_matrices,
        test_data = boot_matrices,
        params = params
      )
      
      # Return bootstrap results
      return(list(
        weights = boot_fit$weights,
        rho = boot_fit$rho,
        converged = boot_fit$converged
      ))
    }
    
    # Run bootstraps (in parallel if possible)
    if (requireNamespace("furrr", quietly = TRUE) && 
        !inherits(future::plan(), "sequential")) {
      boot_results <- furrr::future_map(1:n_bootstraps, run_bootstrap,
                                       .options = furrr::furrr_options(seed = params$parallel$seed))
    } else {
      boot_results <- lapply(1:n_bootstraps, run_bootstrap)
    }
    
    # Extract and summarize bootstrap results
    boot_rhos <- sapply(boot_results, function(x) x$rho)
    boot_converged <- sapply(boot_results, function(x) x$converged)
    
    # Calculate stability of weights
    boot_weights <- lapply(seq_along(matrices), function(i) {
      # Extract weights for block i from all bootstraps
      all_weights <- sapply(boot_results, function(x) x$weights[[i]])
      return(all_weights)
    })
    
    results$bootstrap_results <- list(
      boot_weights = boot_weights,
      boot_rhos = boot_rhos,
      converged = boot_converged,
      mean_boot_rho = mean(boot_rhos),
      sd_boot_rho = sd(boot_rhos)
    )
  }
  
  # Set class for custom methods
  class(results) <- c("mbspls", "list")
  
  return(results)
}
