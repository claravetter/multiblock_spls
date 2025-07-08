#' mbspls: Multi-Block Sparse Partial Least Squares
#'
#' The mbspls package implements Multi-Block Sparse Partial Least Squares (MBSPLS)
#' for integrative analysis of multiple data blocks. The package provides functionality
#' for dimensionality reduction, feature selection, and correlation analysis between
#' multiple high-dimensional datasets.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{mbspls_setup}}: Set up parameters for MBSPLS analysis
#'   \item \code{\link{mbspls_run}}: Run MBSPLS analysis with cross-validation,
#'         permutation testing, and bootstrapping
#'   \item \code{\link{mbs_data}}: Create an MBSPLS data object from a list of matrices
#'   \item \code{\link{read_matlab_input}}: Read MATLAB .mat file with MBSPLS input
#'   \item \code{\link{write_rds_input}}: Write MBSPLS data to RDS file
#' }
#'
#' @section S3 Methods:
#' \itemize{
#'   \item \code{print.mbspls}: Print method for mbspls objects
#'   \item \code{summary.mbspls}: Summary method for mbspls objects
#'   \item \code{coef.mbspls}: Extract model coefficients from an mbspls object
#'   \item \code{predict.mbspls}: Project new data onto the MBSPLS latent space
#'   \item \code{plot.mbspls}: Create plots for MBSPLS results
#' }
#'
#' @useDynLib mbspls, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats cor sd lm
#' @importFrom utils packageVersion
#' @importFrom graphics abline barplot hist legend lines par plot points text
#' @docType package
#' @name mbspls-package
NULL
