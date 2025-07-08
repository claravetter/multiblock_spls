library(devtools)
# Load the package
load_all()

# Create synthetic data to test
set.seed(123)
X1 <- matrix(rnorm(50*10), 50, 10)
X2 <- matrix(rnorm(50*15), 50, 15)
data <- mbs_data(block1 = X1, block2 = X2)

# Set up parameters
params <- mbspls_setup(sparsity = c(2, 2))

# Run MBSPLS with cross-validation
cat("Starting MBSPLS run...\n")
result <- tryCatch({
  mbspls_run(data, params, cv_folds = 2)
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

# Check if successful
if (!is.null(result)) {
  cat("Success! mbspls_run completed without errors\n")
  cat("Result is a:", class(result)[1], "object\n")
  cat("Number of blocks:", length(result$weights), "\n")
  cat("Correlation:", result$rho, "\n")
  if (!is.null(result$cv_results)) {
    cat("CV mean correlation:", result$cv_results$mean_rho, "\n")
  }
}

# Save workspace for debugging
save.image("mbspls_test.RData")
cat("Test complete. Workspace saved to 'mbspls_test.RData'\n")
