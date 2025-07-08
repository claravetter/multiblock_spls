#' Set up MBSPLS parameters
#'
#' This function configures all parameters needed for Multi-Block Sparse PLS analysis.
#' It sets defaults and validates user inputs.
#'
#' @param sparsity Numeric vector of sparsity parameters for each block.
#'        Values should be between 1 and sqrt(number of features) in each block.
#' @param g_factors Matrix or numeric value specifying the type of optimization.
#'        Default is (1/(num_blocks-1)) * (matrix of 1s - identity matrix).
#' @param convergence_threshold Numeric threshold for convergence. Default: 1e-5.
#' @param max_iterations Maximum number of iterations for the algorithm. Default: 1000.
#' @param correlation_method Method for correlation computation. Default: "pearson".
#'        Other options: "spearman" or "kendall".
#' @param num_components Number of latent variables/components to extract. Default: 2.
#' @param matrix_norm Norm used for multi-block correlation. Default: "fro".
#'        Options: "fro" (Frobenius), 1, 2, Inf.
#' @param print_level Integer controlling verbosity. Higher means more output.
#'        Default: 0 (minimal output).
#' @param ensure_nonzero_blocks Logical, if TRUE, ensures each block has at least 
#'        one non-zero feature, preventing blocks from dropping out. Default: FALSE.
#' @param parallel_settings List of parallel configuration options. Default: NULL.
#'
#' @return A list with all parameters needed for MBSPLS analysis.
#' @export
#'
#' @examples
#' # Setup with default parameters for 3 blocks
#' params <- mbspls_setup(sparsity = c(2, 2, 2))
#'
#' # Custom setup
#' params <- mbspls_setup(
#'   sparsity = c(3, 2, 4),
#'   convergence_threshold = 1e-6,
#'   max_iterations = 2000,
#'   correlation_method = "spearman"
#' )
mbspls_setup <- function(sparsity,
                         g_factors = NULL,
                         convergence_threshold = 1e-5,
                         max_iterations = 1000,
                         correlation_method = c("pearson", "spearman", "kendall"),
                         num_components = 2,
                         matrix_norm = "fro",
                         print_level = 0,
                         ensure_nonzero_blocks = FALSE,
                         parallel_settings = NULL) {
  
  # Validate inputs
  correlation_method <- match.arg(correlation_method)
  
  # Handle matrix_norm - allow both character and numeric inputs
  if (is.numeric(matrix_norm)) {
    # Convert numeric to character for validation
    matrix_norm_char <- as.character(matrix_norm)
    matrix_norm_val <- matrix_norm  # Keep original numeric value
  } else {
    # For character inputs
    valid_norms <- c("fro", "1", "2", "inf")
    if (!matrix_norm %in% valid_norms) {
      warning("Invalid matrix_norm '", matrix_norm, "', using 'fro' instead")
      matrix_norm <- "fro"
    }
    matrix_norm_char <- matrix_norm
    
    # Convert string norms to appropriate values
    matrix_norm_val <- switch(matrix_norm,
                             "fro" = "fro",
                             "inf" = Inf,
                             "1" = 1,
                             "2" = 2,
                             matrix_norm)
  }
  
  # Make sure we have a valid value
  if (is.null(matrix_norm_val)) {
    matrix_norm_val <- "fro"
    warning("Invalid matrix_norm specified, using 'fro' instead")
  }
  
  # Create default g_factors if not provided
  num_blocks <- length(sparsity)
  if (is.null(g_factors)) {
    g_factors <- (1/(num_blocks - 1)) * (matrix(1, num_blocks, num_blocks) - diag(num_blocks))
  } else if (length(g_factors) == 1) {
    # Use a constant value for all off-diagonal elements
    temp_g <- matrix(0, num_blocks, num_blocks)
    temp_g[lower.tri(temp_g) | upper.tri(temp_g)] <- g_factors
    g_factors <- temp_g
  }
  
  # Set up parallel configuration
  if (is.null(parallel_settings)) {
    # Default parallel settings
    available_cores <- min(parallel::detectCores(logical = FALSE), 8)
    parallel_settings <- list(
      workers = available_cores,
      strategy = "multisession",
      seed = TRUE   # for reproducibility
    )
  }
  
  # Return the parameter list
  params <- list(
    sparsity = sparsity,
    gs = g_factors,
    convergence_threshold = convergence_threshold,
    max_iterations = max_iterations,
    correlation_method = correlation_method,
    num_components = num_components,
    matrix_norm = matrix_norm_val,
    ensure_nonzero_blocks = ensure_nonzero_blocks,
    print_level = print_level,
    parallel = parallel_settings
  )
  
  # Add validation info and class
  class(params) <- c("mbspls_params", "list")
  return(params)
}

#' Print MBSPLS parameters
#'
#' @param x An mbspls_params object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the input object.
#' @export
print.mbspls_params <- function(x, ...) {
  cat("MBSPLS parameters:\n")
  cat("  Number of blocks:", length(x$sparsity), "\n")
  cat("  Sparsity parameters:", paste(x$sparsity, collapse = ", "), "\n")
  cat("  Convergence threshold:", x$convergence_threshold, "\n")
  cat("  Maximum iterations:", x$max_iterations, "\n")
  cat("  Correlation method:", x$correlation_method, "\n")
  cat("  Number of components:", x$num_components, "\n")
  cat("  Matrix norm:", x$matrix_norm, "\n")
  cat("  Parallel workers:", x$parallel$workers, "\n")
  
  invisible(x)
}
