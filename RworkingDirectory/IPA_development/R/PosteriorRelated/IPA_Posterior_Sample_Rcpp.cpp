#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sampleRcpp( NumericVector x,
                          int size,
                          bool replace,
                          NumericVector prob = NumericVector::create()
) {
  NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}