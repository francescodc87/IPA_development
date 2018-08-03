#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector testHello(CharacterVector x) {
  return x;
}

/*** R
testHello("Hello World!")
*/
