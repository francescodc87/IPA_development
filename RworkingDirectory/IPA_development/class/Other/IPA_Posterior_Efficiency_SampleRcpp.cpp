
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
NumericVector sampleRcppExport( NumericVector x,
                                  int size,
                                  bool replace,
                                  NumericVector prob = NumericVector::create()
)
{
  NumericVector result = RcppArmadillo::sample(x, size, replace, prob);
  // printOutVec (result);
  return result;
}