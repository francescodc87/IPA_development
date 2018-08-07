#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// // [[Rcpp::export]]
// arma::sp_mat test(arma::mat x,NumericVector y){
//   arma::sp_mat x_spmat = arma::sp_mat(x);
//   arma::rowvec y_rvec = as<arma::rowvec>(y);
//   x_spmat.submat(1,0,1,3) = y_rvec;
//   return x_spmat;
// }

// // [[Rcpp::export]]
// arma::sp_rowvec test(arma::sp_mat x,int i){
//   arma::sp_rowvec a = x.row(i);
//   return a;
// }

// // [[Rcpp::export]]
// arma::sp_mat test(arma::sp_mat x,arma::rowvec y){
//   x.row(2) = y;
//   return x;
// }


// // [[Rcpp::export]]
// arma::sp_mat computeIntRatio(arma::mat x, arma::rowvec massRatios){
//   arma::sp_mat Iso = arma::sp_mat(x);
//   for(int i = 0; i < Iso.n_rows; i++){
//     Iso.row(i).operator*=(massRatios.at(i));
//   }
//   return Iso;
// }
// // [[Rcpp::export]]
// arma::sp_mat spMatSubView(arma::mat m, arma::rowvec rows){
//   arma::sp_mat x = arma::sp_mat(m);
//   arma::sp_mat result = arma::sp_mat(rows.n_elem, x.n_cols);
//   for(int i = 0; i < rows.n_elem; i++){
//     int rowNum = rows.at(i);
//     arma::sp_mat x_sub = x.row(rowNum);
//     result.row(i) = x_sub;
//   }
//   return result;
// }\

// // [[Rcpp::export]]
// arma::uvec test(arma::mat m, int a){
//   arma::sp_mat x = arma::sp_mat(m);
//   arma::uvec result = arma::find(x == a);
//   return result;
// }

// // [[Rcpp::export]]
// int test(arma::mat m){
//   arma::sp_mat x = arma::sp_mat(m);
//   return x.n_elem;
// }

// // [[Rcpp::export]]
// arma::sp_mat test(arma::mat x,arma::mat y){
//   arma::sp_mat x_spmat = arma::sp_mat(x);
//   arma::sp_mat y_spmat = arma::sp_mat(y);
//   return x_spmat % y_spmat;
// }

// [[Rcpp::export]]
arma::sp_mat test(arma::mat x,int y){
  arma::sp_mat x_spmat = arma::sp_mat(x);
  return x_spmat.
}
