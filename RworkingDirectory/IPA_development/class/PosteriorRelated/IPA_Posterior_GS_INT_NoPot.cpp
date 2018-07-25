
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#include "/usr/include/valgrind/callgrind.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]



NumericVector getRowsAndDoColsum(NumericMatrix mat, NumericVector rowIdices)
{
  NumericVector sum(mat.ncol());
  for(auto n : rowIdices ){
    sum += mat.row(n-1);
  }
  return sum;
}


double_t vecSum(NumericVector vec){
  double sum = 0.0;
  for (auto n : vec){
    sum += n;
  }
  return sum;
}


NumericVector initVec(int num)
{
  NumericVector result(num);
  for(int i = 1; i <= num; i++){
    result[i-1] = i;
  }
  return result;
}

NumericVector sampleRcppUnExport( NumericVector x,
                                  int size,
                                  bool replace,
                                  NumericVector prob = NumericVector::create()
)
{
  NumericVector result = RcppArmadillo::sample(x, size, replace, prob);
  // printOutVec (result);
  return result;
}


arma::mat computePost(int itNum, 
                      int massNum, 
                      int compNum, 
                      int burn, 
                      arma::mat M){
  int valid = itNum - burn;
  arma::mat  posterior_mat = arma::zeros<arma::mat>(massNum,compNum);
  for (int i = 0; i < massNum; i++){
    arma::rowvec massSampCount_rvec = arma::zeros<arma::rowvec>(compNum);
    arma::rowvec col_rvec = M.col(i).t();
    arma::rowvec massSamp_rvec = col_rvec.subvec(burn, itNum - 1);
    for (int o = 0; o < valid; o++){
      int sampCompIdx = massSamp_rvec.at(o) - 1;
      massSampCount_rvec.at(sampCompIdx) +=  1;
    }
    arma::rowvec posteriorPerMass = massSampCount_rvec / valid;
    posterior_mat.row(i) = posteriorPerMass;
  }
  
  return posterior_mat;
  
}

NumericVector listEleIntoVec(List ls, int idx){
  return ls[idx];
}

arma::uvec setDiff(arma::uvec x, arma::uvec y){
  arma::uvec& result = x;
  x = arma::unique(x);
  y = arma::unique(y);
  
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uvec q1 = arma::find(x == y[j]);
    if (!q1.empty()) {
      result.shed_row(q1(0));
    }
  }
  
  return result;
}

arma::mat computeIntRatio(arma::mat Int, arma::rowvec massRatios){
  
  int rowNum = Int.n_rows;
  arma::sp_mat Int_spmat = arma::sp_mat(Int);
  arma::sp_mat::row_iterator start = Int_spmat.begin_row(0);  // start of column 1
  arma::sp_mat::row_iterator end = Int_spmat.end_row(rowNum - 1);    //   end of column 3
  
  for(arma::sp_mat::row_iterator it = start; it != end; ++it)
  { 
    int idx = it.row();
    double x = massRatios.at(idx);
    (*it) = (*it) * x;
  }
  return arma::mat(Int_spmat);
}

// void printOutUvec(arma::uvec x){
//   Rcout << "UVector" << std::endl;
//   x.t().raw_print(std::cout);
// }
// 
// void printOutRvec(arma::rowvec x){
//   Rcout << "RVector" << std::endl;
//   x.raw_print(std::cout);
// }
// void printOutNvec(NumericVector x){
//   for(int i = 0; i < x.length(); i++){
//     Rcout << "NumericVector, " << x.at(i) << std::endl;
//   }
// }


// [[Rcpp::export]]
NumericMatrix GibbsSampling_Int_NoPot(  List removal,
                                        NumericVector sampcomp,
                                        NumericVector potBio,
                                        NumericMatrix add,
                                        NumericMatrix iso,
                                        NumericMatrix bio,
                                        NumericVector massInt,
                                        NumericMatrix prior,
                                        double delAdd, 
                                        double delIso,
                                        double delBio,
                                        double ratioToll,
                                        int itNum,
                                        int burn,
                                        bool v
){
  CALLGRIND_START_INSTRUMENTATION;
  CALLGRIND_TOGGLE_COLLECT;
  bool replace = false;
  int massNum = prior.rows();
  int compNum = prior.cols();
  NumericVector seqMass = initVec(massNum);
  NumericVector seqComp = initVec(compNum);
  arma::mat add_mat = as<arma::mat>(add);
  arma::mat iso_mat = as<arma::mat>(iso);
  arma::mat bio_mat = as<arma::mat>(bio);
  arma::uvec seqMass_uvec = as<arma::uvec>(seqMass);
  arma::rowvec potBio_rvec = as<arma::rowvec>(potBio);
  arma::rowvec sampcomp_rvec = as<arma::rowvec>(sampcomp);
  arma::rowvec massInt_rvec = as<arma::rowvec>(massInt);
  arma::mat allsampcomp_mat = arma::zeros<arma::mat>(itNum,massNum);
  for(int i = 0; i < itNum; i++){
    NumericVector ordine(massNum);
    ordine = sampleRcppUnExport(seqMass, massNum, replace);
    
    for(auto n : ordine ){
      int thism = (n - 1);
      
      // assignmented compounds of masses retained
      NumericVector remMass = listEleIntoVec(removal, thism);
      arma::uvec remMass_uvec = as<arma::uvec>(remMass);
      arma::uvec retainMass_uvec = setDiff(seqMass_uvec,remMass_uvec);
      arma::uvec retainComp_uvec;
      
      if(!retainMass_uvec.is_empty()){
        retainMass_uvec = retainMass_uvec - 1;
        arma::rowvec retainComp_rvec = sampcomp_rvec.elem(retainMass_uvec).t();
        NumericVector retainComp = wrap(retainComp_rvec);
        retainComp = retainComp - 1;
        retainComp_uvec = as<arma::uvec>(retainComp);
      }
      
      // counting adductrelations
      arma::rowvec pAdd_rvec = sum(add_mat.rows(retainComp_uvec),0);
      
      // counting isotoperelations
      // in R:
      // tmp <- matrix(Iso[sampcomp[-ind.rem[[thism]]],],ncol = Nc)*(Int[thism]/Int[-ind.rem[[thism]]])
      //   ind.ones <- which((tmp>=ratio.toll) & (tmp <=(1/ratio.toll))) 
      //   tmp[ind.ones]<-1
      // tmp[tmp!=1] <-0
      // p.iso<-colSums(tmp)
      arma::mat IsoSub_mat = iso_mat.rows(retainComp_uvec);
      arma::rowvec massIntRatio_rvec = massInt_rvec.elem(retainMass_uvec).t();
      IsoSub_mat = computeIntRatio(IsoSub_mat, massIntRatio_rvec);
      arma::uvec ids_one = find(IsoSub_mat >= ratioToll && IsoSub_mat <= 1/ratioToll);   
      arma::uvec ids_all = find_finite(IsoSub_mat);
      arma::uvec ids_zero = setDiff(ids_all,ids_one);
      IsoSub_mat.elem(ids_one).fill(1);       
      IsoSub_mat.elem(ids_zero).fill(0);
      arma::rowvec pIso_rvec = sum(IsoSub_mat,0);
      
      
      // counting biotransformations relations
      int remComp_bio = sampcomp_rvec.at(thism) - 1;
      arma::rowvec remBio_rvec = bio_mat.row(remComp_bio);
      arma::rowvec pBio_rvec = potBio_rvec - remBio_rvec;
      
      // adding penalities, ##only into add relation and iso relation
      
      // normalising with deltas
      pAdd_rvec = (pAdd_rvec + delAdd) / sum(pAdd_rvec + delAdd);
      pIso_rvec = (pIso_rvec + delIso) / sum(pIso_rvec + delIso);
      pBio_rvec = (pBio_rvec + delBio) / sum(pBio_rvec + delBio);
      
      // merging scores, dot product
      arma::rowvec prior_rvec  = prior.row(thism);
      arma::rowvec post_rvec = prior_rvec % pAdd_rvec % pIso_rvec % pBio_rvec;
      
      // eliminate negative value
      // post[post < 0] = 0;
      
      // normalise posterior probability
      arma::rowvec posterior_rvec = post_rvec / sum(post_rvec);
      
      // get the origin sampling result of mass 'thism'
      int oldSamp = sampcomp_rvec.at(thism);
      
      // use normalised probability to re-sample mass 'thism'
      NumericVector posterior = wrap(posterior_rvec);
      NumericVector newSampVec = sampleRcppUnExport(seqComp, 1, true, posterior);
      int newSamp = newSampVec[0];
      
      // update sampcomp_rvec
      sampcomp_rvec.at(thism) = newSamp;
      
      // if re-sample works, then
      if(oldSamp!=newSamp){
        potBio_rvec = potBio_rvec - bio_mat.row(oldSamp-1) + bio_mat.row(newSamp-1);
      }
    }
    
    // when finishing each iteration 
    // get the whole new sampcomp vector, where each element(mess) get resampled using its posterior probability vector
    allsampcomp_mat.submat(i,0,i,massNum-1) = sampcomp_rvec;
    
    if(v){
      Rcout << "Computing Posterior in Rcpp, " << (i * 100) / itNum <<"%" << std::endl;
    }
  }
  
  //calculate posterior probability using allsampcomp which is a sampling distribution matrix (row num: no.its; col num: m, which is the number of masses)
  arma::mat result_mat = computePost(itNum, massNum, compNum, burn, allsampcomp_mat);
  NumericMatrix result = wrap(result_mat);
  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_DUMP_STATS;
  CALLGRIND_STOP_INSTRUMENTATION;
  // decide which kind of value to return, temporarily ignored
  // if(allsamp){
  //   List outList;
  //   outList["post"] = posterior;
  //   outList["allsampcomp"] = sampcomp;
  //   return outList;
  // }else{
  return result;  
  // return allsampcomp;
  // }
}

