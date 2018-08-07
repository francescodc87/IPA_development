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

arma::rowvec spMatGetDenseRowVec(arma::sp_mat x, int rowNum){
  arma::sp_mat::const_row_iterator xStart = x.begin_row(rowNum);
  arma::sp_mat::const_row_iterator xEnd   = x.end_row(rowNum);
  arma::rowvec result = arma::zeros<arma::rowvec>(x.n_cols);
  for(arma::sp_mat::const_row_iterator it = xStart; it != xEnd; ++it)
  { 
    double p = (*it);
    if(p != 0){
      int colNum = it.col();
      result.at(colNum) = p;
    }
  }
  return result;
}

int binSpMatGetSingleColSum(arma::sp_mat bin, int colNum, arma::uvec subRow){
  arma::sp_mat::const_col_iterator binStart = bin.begin_col(colNum);
  arma::sp_mat::const_col_iterator binEnd   = bin.end_col(colNum);
  NumericVector nonZeroRow;
  for(arma::sp_mat::const_col_iterator it = binStart; it != binEnd; ++it){
    int value = (*it);
    if(value != 0){
      int rowNum = it.row();
      arma::uvec findInSubRow_uvec = arma::find(subRow == rowNum);
      if(!findInSubRow_uvec.is_empty()){
        nonZeroRow.push_back(rowNum);
      }
    }
  }
  if(nonZeroRow.length() == 0){
    return 0;
  }else{
    return nonZeroRow.length();
  }
}

int isoSpMatGetSingleColSum(arma::sp_mat iso, int colNum, arma::uvec subRow, arma::rowvec intRatio, double ratioToll){
  arma::sp_mat::const_col_iterator isoStart = iso.begin_col(colNum);
  arma::sp_mat::const_col_iterator isoEnd   = iso.end_col(colNum);
  NumericVector qualifiedRow;
  for(arma::sp_mat::const_col_iterator it = isoStart; it != isoEnd; ++it){
    double value = (*it);
    if(value != 0){
      int rowNum = it.row();
      double massIntRatio = intRatio.at(rowNum);
      double tmp = value * massIntRatio;
      if(( tmp>= ratioToll) or (tmp <= 1 / ratioToll)){
        arma::uvec findInSubRow_uvec = arma::find(subRow == rowNum);
        if(!findInSubRow_uvec.is_empty()){
          qualifiedRow.push_back(rowNum);
        }
      }
    }
  }
  if(qualifiedRow.length() == 0){
    return 0;
  }else{
    return qualifiedRow.length();
  }
}

// void printOutUvec(arma::uvec x){
//   Rcout << "UVector" << std::endl;
//   x.t().raw_print(std::cout);
// }
void printOutRvec(arma::rowvec x){
  // Rcout << "RVector" << std::endl;
  x.raw_print(std::cout);
}
// void printOutNvec(NumericVector x){
//   for(int i = 0; i < x.length(); i++){
//     Rcout << "NumericVector, " << x.at(i) << std::endl;
//   }
// }

// [[Rcpp::export]]
NumericMatrix GibbsSampling_Efficiency( List removal,
                                        arma::sp_mat add_spMat,
                                        arma::sp_mat iso_spMat,
                                        arma::sp_mat bio_spMat,
                                        NumericVector massInt,
                                        arma::sp_mat prior_spMat,
                                        double delAdd,
                                        double delIso,
                                        double delBio,
                                        double ratioToll,
                                        int itNum,
                                        int burn,
                                        bool v
){int massNum = prior_spMat.n_rows;
  int fomuNum = prior_spMat.n_cols;
  NumericVector seqMass = initVec(massNum);
  NumericVector seqFomu = initVec(fomuNum);
  arma::uvec seqMass_uvec = as<arma::uvec>(seqMass);
  arma::rowvec massInt_rvec = as<arma::rowvec>(massInt);
  arma::mat allSampFomu_mat = arma::zeros<arma::mat>(itNum,massNum);
  NumericVector sampFomu;
  for(int i = 0; i < massNum; i++){
    arma::rowvec priorRow_rvec = spMatGetDenseRowVec(prior_spMat, i);
    NumericVector priorRow = wrap(priorRow_rvec);
    NumericVector priorRowSample = sampleRcppUnExport(seqFomu, 1, true, priorRow);
    sampFomu.push_back(priorRowSample.at(0));
  }
  sampFomu = sampFomu - 1;
  arma::uvec sampFomu_uvec = as<arma::uvec>(sampFomu);
  arma::rowvec potBio_rvec = arma::zeros<arma::rowvec>(fomuNum);
  for(int i = 0; i < fomuNum; i++){
    potBio_rvec.at(i) = binSpMatGetSingleColSum(bio_spMat, i, sampFomu_uvec);
  }
  arma::rowvec sampFomu_rvec = as<arma::rowvec>(sampFomu);
  for(int i = 0; i < itNum; i++){
    NumericVector ordine(massNum);
    ordine = sampleRcppUnExport(seqMass, massNum, false);
    for(auto n : ordine ){
      int thism = (n - 1);
      // assignmented compounds of masses retained
      NumericVector remMass = listEleIntoVec(removal, thism);
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
      arma::rowvec pAdd_rvec = arma::zeros<arma::rowvec>(fomuNum);
      for(int i = 0; i < fomuNum; i++){
        pAdd_rvec.at(i) = binSpMatGetSingleColSum(add_spMat, i, retainFomu_uvec);
      }

      // counting isotoperelations
      // in R:
      // tmp <- matrix(Iso[sampcomp[-ind.rem[[thism]]],],ncol = Nc)*(Int[thism]/Int[-ind.rem[[thism]]])
      //   ind.ones <- which((tmp>=ratio.toll) & (tmp <=(1/ratio.toll)))
      //   tmp[ind.ones]<-1
      // tmp[tmp!=1] <-0
      // p.iso<-colSums(tmp)
      arma::rowvec pIso_rvec = arma::zeros<arma::rowvec>(fomuNum);
      arma::rowvec massIntRatio_rvec = massInt_rvec.elem(retainMass_uvec).t();
      for(int i = 0; i < fomuNum; i++){
        pAdd_rvec.at(i) = isoSpMatGetSingleColSum(iso_spMat, i, retainFomu_uvec, massIntRatio_rvec, ratioToll);
      }

      // counting biotransformations relations
      int remComp_bio = sampFomu_rvec.at(thism);
      arma::rowvec remBio_rvec = spMatGetDenseRowVec(bio_spMat, remComp_bio);
      arma::rowvec pBio_rvec = potBio_rvec - remBio_rvec;
      
      // adding penalities, ##only into add relation and iso relation

      // normalising with deltas
      pAdd_rvec = (pAdd_rvec + delAdd) / sum(pAdd_rvec + delAdd);
      pIso_rvec = (pIso_rvec + delIso) / sum(pIso_rvec + delIso);
      pBio_rvec = (pBio_rvec + delBio) / sum(pBio_rvec + delBio);

      // merging scores, dot product
      arma::rowvec prior_rvec  = spMatGetDenseRowVec(prior_spMat, thism);
      arma::rowvec post_rvec = prior_rvec % pAdd_rvec % pIso_rvec % pBio_rvec;

      // normalise posterior probability
      arma::rowvec posterior_rvec = post_rvec / sum(post_rvec);

      // get the origin sampling result of mass 'thism'
      int oldSamp = sampFomu_rvec.at(thism);

      // use normalised probability to re-sample mass 'thism'
      NumericVector posterior = wrap(posterior_rvec);
      NumericVector newSampVec = sampleRcppUnExport(seqFomu, 1, true, posterior);
      int newSamp = newSampVec[0];

      // update sampcomp_rvec
      sampFomu_rvec.at(thism) = newSamp;

      // if re-sample works, then
      if(oldSamp!=newSamp - 1){
        arma::rowvec potBioRem_rvec = spMatGetDenseRowVec(bio_spMat, oldSamp);
        arma::rowvec potBioAdd_rvec = spMatGetDenseRowVec(bio_spMat, newSamp-1);
        potBio_rvec = potBio_rvec - potBioRem_rvec + potBioAdd_rvec;
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
  // arma::sp_mat result_spMat = arma::sp_mat(result_mat);
  
  // decide which kind of value to return, temporarily ignored
  // if(allsamp){
  //   List outList;
  //   outList["post"] = posterior;
  //   outList["allsampcomp"] = sampcomp;
  //   return outList;
  // }else{
  return wrap(result_mat);
  // return allsampcomp;
  // }
}