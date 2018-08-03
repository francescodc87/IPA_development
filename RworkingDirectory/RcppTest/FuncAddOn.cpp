#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector addTwo(NumericVector x, NumericVector y) {
  return x + y;
}

/*** R
addTwo(3, 5)
*/
