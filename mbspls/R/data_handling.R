#' Read MATLAB .mat file with MBSPLS input
#'
#' This function reads a MATLAB .mat file containing input data for MBSPLS analysis
#' and converts it to an R-friendly format.
#'
#' @param path Character string with the path to the .mat file
#' @param drop_singleton Logical; if TRUE, singleton dimensions are dropped. Default: TRUE.
#'
#' @return A list containing the MBSPLS input data with components:
#'   \item{matrices}{List of data matrices (blocks)}
#'   \item{setup}{List of setup parameters}
#'   \item{metadata}{Additional metadata if available}
#'
#' @importFrom R.matlab readMat
#' @export
#'
#' @examples
#' \dontrun{
#' # Read MBSPLS data from a MATLAB file
#' data <- read_matlab_input("path/to/data.mat")
#' }
read_matlab_input <- function(path, drop_singleton = TRUE) {
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("Package 'R.matlab' is required for reading MATLAB files")
  }
  
  # Read the .mat file
  mat_data <- R.matlab::readMat(path)
  
  # Check for key components and convert structure
  result <- list(
    matrices = list(),
    setup = list(),
    metadata = list()
  )
  
  # Check for and process 'input' structure
  if ("input" %in% names(mat_data)) {
    input_struct <- mat_data$input
    
    # Check for X and Y matrices in input
    if (!is.null(input_struct$X)) {
      # Convert X to an R matrix, dropping singleton dimensions if requested
      X <- if (drop_singleton) drop_singleton_dims(input_struct$X) else input_struct$X
      result$matrices[[1]] <- X
      
      # Set name if available, otherwise use default
      if (!is.null(input_struct$X.name)) {
        names(result$matrices)[1] <- input_struct$X.name
      } else {
        names(result$matrices)[1] <- "block1"
      }
    }
    
    # Check for Y matrix
    if (!is.null(input_struct$Y)) {
      Y <- if (drop_singleton) drop_singleton_dims(input_struct$Y) else input_struct$Y
      result$matrices[[2]] <- Y
      
      # Set name if available, otherwise use default
      if (!is.null(input_struct$Y.name)) {
        names(result$matrices)[2] <- input_struct$Y.name
      } else {
        names(result$matrices)[2] <- "block2"
      }
    }
    
    # Check for multiple blocks (Xs) in input
    if (!is.null(input_struct$Xs)) {
      Xs <- input_struct$Xs
      
      # Handle cell array of matrices
      for (i in seq_along(Xs)) {
        block <- if (drop_singleton) drop_singleton_dims(Xs[[i]]) else Xs[[i]]
        result$matrices[[i]] <- block
        
        # Set name if available, otherwise use default
        if (!is.null(input_struct$Xs.names) && i <= length(input_struct$Xs.names)) {
          names(result$matrices)[i] <- input_struct$Xs.names[[i]]
        } else {
          names(result$matrices)[i] <- paste0("block", i)
        }
      }
    }
    
    # Extract other important parameters from input
    param_fields <- c("permutation_testing", "bootstrap_testing", "hyperopt_params")
    for (field in param_fields) {
      if (!is.null(input_struct[[field]])) {
        result$metadata[[field]] <- input_struct[[field]]
      }
    }
  }
  
  # Check for 'setup' structure
  if ("setup" %in% names(mat_data)) {
    setup_struct <- mat_data$setup
    
    # Extract setup parameters
    setup_fields <- c("analysis_folder", "scratch_space", "n_cores")
    for (field in setup_fields) {
      if (!is.null(setup_struct[[field]])) {
        result$setup[[field]] <- setup_struct[[field]]
      }
    }
  }
  
  # Check for 'matrices' structure (sometimes data is stored here directly)
  if ("matrices" %in% names(mat_data) && length(result$matrices) == 0) {
    matrices_struct <- mat_data$matrices
    
    # Extract matrices - specific handling depends on your MATLAB structure
    # This is a simplified example:
    for (name in names(matrices_struct)) {
      if (is.matrix(matrices_struct[[name]]) || is.array(matrices_struct[[name]])) {
        block <- if (drop_singleton) drop_singleton_dims(matrices_struct[[name]]) else matrices_struct[[name]]
        n_blocks <- length(result$matrices) + 1
        result$matrices[[n_blocks]] <- block
        names(result$matrices)[n_blocks] <- name
      }
    }
  }
  
  # Add class for easier handling
  class(result) <- c("mbs_data", "list")
  
  return(result)
}

#' Write MBSPLS data to RDS file
#'
#' This function saves MBSPLS data as an RDS file for future use in R.
#'
#' @param data An mbs_data object or list with matrices and settings
#' @param path Character string with the path to save the RDS file
#'
#' @return Invisibly returns the path where the data was saved
#' @export
#'
#' @examples
#' \dontrun{
#' # Save MBSPLS data to an RDS file
#' write_rds_input(mbspls_data, "path/to/data.rds")
#' }
write_rds_input <- function(data, path) {
  # Make sure the path ends with .rds
  if (!grepl("\\.rds$", path, ignore.case = TRUE)) {
    path <- paste0(path, ".rds")
  }
  
  # Save the data
  saveRDS(data, file = path)
  
  # Return the path invisibly
  invisible(path)
}

#' Drop singleton dimensions from arrays
#'
#' Helper function to remove singleton dimensions from arrays or matrices.
#'
#' @param x An array or matrix
#'
#' @return The input with singleton dimensions removed
#' @keywords internal
drop_singleton_dims <- function(x) {
  if (is.array(x) || is.matrix(x)) {
    # Find non-singleton dimensions
    non_singleton <- which(dim(x) > 1)
    
    if (length(non_singleton) < length(dim(x))) {
      # Keep only non-singleton dimensions
      return(array(x, dim = dim(x)[non_singleton]))
    }
  }
  return(x)
}

#' Create an MBSPLS data object from a list of matrices
#'
#' @param ... Named matrices to include in the MBSPLS analysis
#' @param metadata Optional list with additional metadata
#'
#' @return An object of class 'mbs_data'
#' @export
#'
#' @examples
#' # Create MBSPLS data from matrices
#' X1 <- matrix(rnorm(100*20), 100, 20)
#' X2 <- matrix(rnorm(100*30), 100, 30)
#' X3 <- matrix(rnorm(100*10), 100, 10)
#' data <- mbs_data(proteome = X1, metabolome = X2, clinical = X3)
mbs_data <- function(..., metadata = list()) {
  # Collect matrices passed as named arguments
  matrices <- list(...)
  
  # Check that all arguments are matrices with the same number of rows
  if (length(matrices) == 0) {
    stop("At least one matrix must be provided")
  }
  
  for (i in seq_along(matrices)) {
    if (!is.matrix(matrices[[i]])) {
      stop("All arguments must be matrices")
    }
    
    if (i > 1 && nrow(matrices[[i]]) != nrow(matrices[[1]])) {
      stop("All matrices must have the same number of rows (samples)")
    }
  }
  
  # Create the result object
  result <- list(
    matrices = matrices,
    metadata = metadata
  )
  
  # Add class
  class(result) <- c("mbs_data", "list")
  
  return(result)
}

#' Print method for mbs_data objects
#'
#' @param x An mbs_data object
#' @param ... Further arguments passed to or from other methods
#'
#' @return Invisibly returns the input object
#' @export
print.mbs_data <- function(x, ...) {
  cat("Multi-Block Data:\n")
  cat("  Number of samples:", nrow(x$matrices[[1]]), "\n")
  cat("  Number of blocks:", length(x$matrices), "\n\n")
  
  cat("Blocks:\n")
  for (i in seq_along(x$matrices)) {
    cat(sprintf("  %s: %d variables\n", 
                names(x$matrices)[i], 
                ncol(x$matrices[[i]])))
  }
  
  # Print metadata summary if available
  if (length(x$metadata) > 0) {
    cat("\nMetadata fields:", paste(names(x$metadata), collapse = ", "), "\n")
  }
  
  invisible(x)
}

#' Subset method for mbs_data objects
#'
#' @param x An mbs_data object
#' @param i Indices for subsetting rows (samples)
#' @param j Block names or indices for subsetting columns (blocks)
#' @param ... Not used
#'
#' @return A subset mbs_data object
#' @export
`[.mbs_data` <- function(x, i, j, ...) {
  # Handle missing indices
  if (missing(i)) i <- seq_len(nrow(x$matrices[[1]]))
  if (missing(j)) j <- seq_along(x$matrices)
  
  # If j is character, match with block names
  if (is.character(j)) {
    j <- match(j, names(x$matrices))
    if (any(is.na(j))) {
      stop("Some block names not found: ", 
           paste(j[is.na(j)], collapse = ", "))
    }
  }
  
  # Subset matrices
  matrices <- x$matrices[j]
  for (k in seq_along(matrices)) {
    matrices[[k]] <- matrices[[k]][i, , drop = FALSE]
  }
  
  # Create new mbs_data object
  result <- list(
    matrices = matrices,
    metadata = x$metadata
  )
  
  class(result) <- c("mbs_data", "list")
  return(result)
}
