
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
// #include "/usr/include/valgrind/callgrind.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


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
      int sampCompIdx = massSamp_rvec.at(o);
      massSampCount_rvec.at(sampCompIdx) +=  1;
    }
    arma::rowvec posteriorPerMass = massSampCount_rvec / valid;
    posterior_mat.row(i) = posteriorPerMass;
  }
  return posterior_mat;
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

arma::mat computeIntRatio(arma::mat Iso, arma::rowvec massRatios){
  for(int i = 0; i < Iso.n_rows; i++){
    Iso.row(i).operator*=(massRatios.at(i));
  }
  return Iso;
}

NumericVector spMatColSum(arma::sp_mat x, arma::uvec y){
  int rowNum = y.n_elem;
  int colNum = x.n_cols;
  NumericVector result(colNum);
  for (int i = 0; i < rowNum; i++){
    int rowIdx = y.at(i);
    arma::sp_rowvec rvec = x.row(rowIdx);
    NumericVector row = wrap(rvec);
    result = result + row;
  }
  return result;
}

// void printOutUvec(arma::uvec x){
//   Rcout << "UVector" << std::endl;
//   x.t().raw_print(std::cout);
// }
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
  // CALLGRIND_START_INSTRUMENTATION;
  // CALLGRIND_TOGGLE_COLLECT;
  int massNum = prior.rows();
  int fomuNum = prior.cols();
  NumericVector seqMass = initVec(massNum);
  NumericVector seqFomu = initVec(fomuNum);
  arma::mat add_mat = as<arma::mat>(add);
  arma::mat iso_mat = as<arma::mat>(iso);
  arma::mat bio_mat = as<arma::mat>(bio);
  arma::mat prior_mat = as<arma::mat>(prior);
  arma::uvec seqMass_uvec = as<arma::uvec>(seqMass);
  arma::rowvec massInt_rvec = as<arma::rowvec>(massInt);
  arma::mat allSampFomu_mat = arma::zeros<arma::mat>(itNum,massNum);
  
  NumericVector sampFomu;
  for(int i = 0; i < massNum; i++){
    arma::rowvec priorRow_rvec = prior_mat.row(i);
    NumericVector priorRow = wrap(priorRow_rvec);
    NumericVector priorRowSample = sampleRcppUnExport(seqFomu, 1, true, priorRow);
    sampFomu.push_back(priorRowSample.at(0));
  }
  sampFomu = sampFomu - 1;
  arma::uvec sampFomu_uvec = as<arma::uvec>(sampFomu);
  arma::rowvec potBio_rvec = sum(bio_mat.rows(sampFomu_uvec), 0);
  arma::rowvec sampFomu_rvec = as<arma::rowvec>(sampFomu);
  NumericVector ordine(massNum);
  for(int i = 0; i < itNum; i++){
    ordine = sampleRcppUnExport(seqMass, massNum, false);
    ordine = ordine - 1;
    for(auto n : ordine ){
      int thism = n;
      // assignmented compounds of masses retained
      NumericVector remMass = removal.at(thism);
      arma::uvec remMass_uvec = as<arma::uvec>(remMass);
      arma::uvec retainMass_uvec = setDiff(seqMass_uvec,remMass_uvec);
      arma::uvec retainFomu_uvec;
      if(!retainMass_uvec.is_empty()){
        retainMass_uvec = retainMass_uvec - 1;
        arma::rowvec retainFomu_rvec = sampFomu_rvec.elem(retainMass_uvec).t();
        NumericVector retainFomu = wrap(retainFomu_rvec);
        retainFomu_uvec = as<arma::uvec>(retainFomu);
      }
      
      // counting adductrelations
      arma::rowvec pAdd_rvec = sum(add_mat.rows(retainFomu_uvec),0);
      
      // counting isotoperelations
      // in R:
      // tmp <- matrix(Iso[sampcomp[-ind.rem[[thism]]],],ncol = Nc)*(Int[thism]/Int[-ind.rem[[thism]]])
      //   ind.ones <- which((tmp>=ratio.toll) & (tmp <=(1/ratio.toll))) 
      //   tmp[ind.ones]<-1
      // tmp[tmp!=1] <-0
      // p.iso<-colSums(tmp)
      arma::mat IsoSub_mat = iso_mat.rows(retainFomu_uvec);
      arma::rowvec massIntRatio_rvec = massInt_rvec.elem(retainMass_uvec).t();
      IsoSub_mat = computeIntRatio(IsoSub_mat, massIntRatio_rvec);
      arma::uvec ids_one = find(IsoSub_mat >= ratioToll && IsoSub_mat <= 1/ratioToll);
      arma::mat IsoSubNew_mat = arma::zeros<arma::mat>(retainFomu_uvec.n_elem, fomuNum);
      IsoSubNew_mat.elem(ids_one).fill(1);       
      arma::rowvec pIso_rvec = sum(IsoSubNew_mat,0);
      
      // counting biotransformations relations
      int remFomu_bio = sampFomu_rvec.at(thism);
      arma::rowvec remBio_rvec = bio_mat.row(remFomu_bio);
      arma::rowvec pBio_rvec = potBio_rvec - remBio_rvec;
      
      // adding penalities, ##only into add relation and iso relation
      
      // normalising with deltas
      pAdd_rvec = (pAdd_rvec + delAdd) / sum(pAdd_rvec + delAdd);
      pIso_rvec = (pIso_rvec + delIso) / sum(pIso_rvec + delIso);
      pBio_rvec = (pBio_rvec + delBio) / sum(pBio_rvec + delBio);
      
      // merging scores, dot product
      arma::rowvec prior_rvec  = prior_mat.row(thism);
      arma::rowvec post_rvec = prior_rvec % pAdd_rvec % pIso_rvec % pBio_rvec;
      
      // eliminate negative value
      // post[post < 0] = 0;
      
      // normalise posterior probability
      arma::rowvec posterior_rvec = post_rvec / sum(post_rvec);
      
      // get the origin sampling result of mass 'thism'
      int oldSamp = sampFomu_rvec.at(thism);
      
      // use normalised probability to re-sample mass 'thism'
      NumericVector posterior = wrap(posterior_rvec);
      NumericVector newSampVec = sampleRcppUnExport(seqFomu, 1, true, posterior);
      int newSamp = newSampVec[0] - 1;
      
      // update sampcomp_rvec
      sampFomu_rvec.at(thism) = newSamp;
      
      // if re-sample works, then
      if(oldSamp != newSamp){
        potBio_rvec = potBio_rvec - bio_mat.row(oldSamp) + bio_mat.row(newSamp);
      }
    }
    
    // when finishing each iteration 
    // get the whole new sampcomp vector, where each element(mess) get resampled using its posterior probability vector
    allSampFomu_mat.submat(i,0,i,massNum-1) = sampFomu_rvec;
    
    if(v){
      Rcout << "Computing Posterior in Rcpp, " << (i * 100) / itNum <<"%" << std::endl;
    }
  }
  
  //calculate posterior probability using allsampcomp which is a sampling distribution matrix (row num: no.its; col num: m, which is the number of masses)
  arma::mat result_mat = computePost(itNum, massNum, fomuNum, burn, allSampFomu_mat);
  NumericMatrix result = wrap(result_mat);
  // CALLGRIND_TOGGLE_COLLECT;
  // CALLGRIND_DUMP_STATS;
  // CALLGRIND_STOP_INSTRUMENTATION;
  
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


