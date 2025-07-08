// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::vec soft_threshold_cpp(const arma::vec& w, const double lambda) {
  int n = w.n_elem;
  
  // Sort absolute values in descending order
  arma::vec abs_w = arma::abs(w);
  arma::uvec idx = arma::sort_index(abs_w, "descend");
  arma::vec abs_sorted = abs_w(idx);
  
  // Compute cumulative sum
  arma::vec cumsum_abs = arma::cumsum(abs_sorted);
  
  // Find threshold point
  int k = 0;
  for (int i = 0; i < n; i++) {
    if (abs_sorted(i) > (cumsum_abs(i) - lambda) / (i + 1)) {
      k = i + 1;
    } else {
      break;
    }
  }
  
  // If all weights are below threshold, return zeros
  if (k == 0) {
    return arma::zeros<arma::vec>(n);
  }
  
  // Compute threshold value
  double threshold = (cumsum_abs(k-1) - lambda) / k;
  
  // Apply soft-thresholding
  arma::vec result = arma::max(abs_w - threshold, 0.0) % arma::sign(w);
  
  return result;
}

// [[Rcpp::export]]
arma::mat soft_thresh_update_cpp(const arma::mat& X, const arma::mat& Y, 
                                 const arma::vec& w_Y, const double lambda) {
  // Compute the gradient X'Y * w_Y
  arma::vec new_weight = X.t() * (Y * w_Y);
  
  // Apply soft thresholding
  arma::vec thresholded = soft_threshold_cpp(new_weight, lambda);
  
  // Normalize
  double norm_val = arma::norm(thresholded);
  if (norm_val > 1e-10) {
    thresholded = thresholded / norm_val;
  }
  
  return thresholded;
}
