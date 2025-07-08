# Define global variables used in ggplot2 calls to avoid R CMD check notes
utils::globalVariables(c(
  ".data", "LV1", "LV2", "Variable", "Loading", "Block", "Fold", "Correlation", 
  "Lower", "Upper", "Mean", "Original"
))

.onLoad <- function(libname, pkgname) {
  # Make sure internal functions are linked properly
  mbspls_env <- asNamespace(pkgname)
  
  # Register S3 methods directly (without roxygen)
  s3_register <- function(generic, class, method = NULL) {
    if (is.null(method)) {
      method <- paste0(generic, ".", class)
    }
    
    method_fn <- get(method, envir = mbspls_env)
    registerS3method(generic, class, method_fn, envir = mbspls_env)
  }
  
  # Register all needed S3 methods
  s3_register("print", "mbspls")
  s3_register("summary", "mbspls")
  s3_register("coef", "mbspls")
  s3_register("predict", "mbspls")
  s3_register("plot", "mbspls")
  s3_register("print", "mbspls_params")
  s3_register("print", "summary.mbspls")
  s3_register("print", "mbs_data")
  s3_register("[", "mbs_data")
  
  invisible()
}

.onUnload <- function(libpath) {
  library.dynam.unload("mbspls", libpath)
}
