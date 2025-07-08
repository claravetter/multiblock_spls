test_that("mbspls_setup creates valid parameter object", {
  # Test with minimal parameters
  params <- mbspls_setup(sparsity = c(2, 2))
  expect_s3_class(params, "mbspls_params")
  expect_equal(params$sparsity, c(2, 2))
  expect_equal(params$convergence_threshold, 1e-5)
  
  # Test with custom parameters
  params <- mbspls_setup(
    sparsity = c(3, 4, 5),
    convergence_threshold = 1e-6,
    correlation_method = "spearman"
  )
  expect_s3_class(params, "mbspls_params")
  expect_equal(params$sparsity, c(3, 4, 5))
  expect_equal(params$convergence_threshold, 1e-6)
  expect_equal(params$correlation_method, "spearman")
})

test_that("mbs_data creates valid data object", {
  # Create test matrices
  X1 <- matrix(rnorm(50*10), 50, 10)
  X2 <- matrix(rnorm(50*15), 50, 15)
  
  # Test creation
  data <- mbs_data(block1 = X1, block2 = X2)
  expect_s3_class(data, "mbs_data")
  expect_equal(length(data$matrices), 2)
  expect_equal(dim(data$matrices$block1), c(50, 10))
  expect_equal(dim(data$matrices$block2), c(50, 15))
  
  # Test subsetting
  subset_data <- data[1:25, ]
  expect_s3_class(subset_data, "mbs_data")
  expect_equal(dim(subset_data$matrices$block1), c(25, 10))
  
  subset_data <- data[, "block1"]
  expect_s3_class(subset_data, "mbs_data")
  expect_equal(length(subset_data$matrices), 1)
  expect_equal(names(subset_data$matrices), "block1")
})
