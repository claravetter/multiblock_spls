#' Fit MBSPLS on training data and project onto test data
#'
#' This function trains MBSPLS on training data and projects the test data
#' onto the resulting latent space.
#'
#' @param training_data List of training data matrices, each representing a block
#' @param test_data List of test data matrices, each representing a block
#' @param params List of MBSPLS parameters from mbspls_setup()
#' @param initial_weights Optional list of initial weight vectors
#'
#' @return A list containing:
#'   \item{rho}{Correlation value or matrix norm between latent variables}
#'   \item{weights}{Weight vectors for each block}
#'   \item{latent_vars}{Latent variables for test data}
#'   \item{training_latent_vars}{Latent variables for training data}
#'
#' @keywords internal
mbspls_fit_fold <- function(training_data, test_data, params, initial_weights = NULL) {
  
  # Extract parameters
  sparsity <- params$sparsity
  g_factors <- params$gs
  convergence_threshold <- params$convergence_threshold
  max_iterations <- params$max_iterations
  correlation_method <- params$correlation_method
  matrix_norm <- params$matrix_norm
  print_level <- params$print_level
  ensure_nonzero_blocks <- params$ensure_nonzero_blocks
  
  # Perform MBSPLS on training data
  fit <- get("mbspls_fit_core", envir = asNamespace("mbspls"))(
    matrices = training_data,
    sparsity = sparsity,
    g_factors = g_factors,
    convergence_threshold = convergence_threshold,
    max_iterations = max_iterations,
    print_level = print_level,
    ensure_nonzero_blocks = ensure_nonzero_blocks,
    initial_weights = initial_weights
  )
  
  # Project test data using the trained weights
  projection <- get("mbspls_project", envir = asNamespace("mbspls"))(
    data = test_data,
    weights = fit$weights,
    correlation_method = correlation_method,
    matrix_norm = matrix_norm,
    force_positive = TRUE
  )
  
  # Return results
  return(list(
    rho = projection$rho,
    weights = projection$weights,
    latent_vars = projection$latent_vars,
    training_latent_vars = fit$latent_vars,
    converged = fit$converged,
    iterations = fit$iterations
  ))
}
