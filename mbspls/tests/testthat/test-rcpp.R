test_that("Rcpp functions are correctly loaded and working", {
  # Skip if we can't access the Rcpp functions (e.g., during CRAN checks)
  if (!exists("soft_threshold_cpp", envir = asNamespace("mbspls"))) {
    skip("Rcpp functions not available in this environment")
  }
  
  # Create a simple test vector
  test_vec <- c(1.5, 0.3, -2.5, 0.8, -0.1)
  lambda <- 1.0
  
  # Test the function
  result <- soft_threshold_cpp(test_vec, lambda)
  
  # Basic expectations
  expect_is(result, "numeric")
  expect_length(result, length(test_vec))
  
  # Compare with R implementation for small vectors
  soft_thresh_r <- function(w, lambda) {
    abs_w <- abs(w)
    ord <- order(abs_w, decreasing = TRUE)
    abs_sorted <- abs_w[ord]
    
    cumsum_abs <- cumsum(abs_sorted)
    k <- which(abs_sorted > (cumsum_abs - lambda) / seq_along(w))[1]
    
    if (is.na(k)) {
      return(numeric(length(w)))
    }
    
    threshold <- (cumsum_abs[k] - lambda) / k
    pmax(abs_w - threshold, 0) * sign(w)
  }
  
  r_result <- soft_thresh_r(test_vec, lambda)
  expect_equal(result, r_result, tolerance = 1e-6)
})
