
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
// #include <Eigen/Sparse>
// #include "/usr/include/valgrind/callgrind.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// typedef Eigen::SparseMatrix<double, RowMajor> SpMat;

arma::rowvec spVectorGetDense(arma::sp_mat x){
  arma::sp_mat::const_iterator xStart = x.begin();
  arma::sp_mat::const_iterator xEnd   = x.end();
  arma::rowvec result = arma::zeros<arma::rowvec>(x.n_cols);
  for(arma::sp_mat::const_iterator it = xStart; it != xEnd; ++it)
  {
    double p = (*it);
    if(p != 0){
      int colNum = it.col();
      result.at(colNum) = p;
    }
  }
  return result;
}

arma::sp_mat spMatColSum(arma::sp_mat x, arma::rowvec rows){
  arma::sp_mat result = arma::sp_mat(1, x.n_cols);
  for(int i = 0; i < rows.n_elem; i++){
    int rowNum = rows.at(i);
    arma::sp_mat x_sub = x.row(rowNum);
    result = result + x_sub;
  }
  return result;
}

arma::sp_mat spMatSubView(arma::sp_mat x, arma::rowvec rows){
  arma::sp_mat result = arma::sp_mat(rows.n_elem, x.n_cols);
  for(int i = 0; i < rows.n_elem; i++){
    int rowNum = rows.at(i);
    arma::sp_mat x_sub = x.row(rowNum);
    result.row(i) = x_sub;
  }
  return result;
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
  return result;
}


arma::sp_mat computePost( int itNum,
                          int massNum,
                          int compNum,
                          int burn,
                          arma::mat M){
  int valid = itNum - burn;
  arma::sp_mat  posterior_spMat = arma::sp_mat(massNum,compNum);
  for (int i = 0; i < massNum; i++){
    arma::rowvec massSampCount_rvec = arma::zeros<arma::rowvec>(compNum);
    arma::rowvec col_rvec = M.col(i).t();
    arma::rowvec massSamp_rvec = col_rvec.subvec(burn, itNum - 1);
    for (int o = 0; o < valid; o++){
      int sampCompIdx = massSamp_rvec.at(o);
      massSampCount_rvec.at(sampCompIdx) +=  1;
    }
    arma::rowvec posteriorPerMass = massSampCount_rvec / valid;
    posterior_spMat.row(i) = posteriorPerMass;
  }
  return posterior_spMat;
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

List thresholdIntRatio(arma::sp_mat x, arma::rowvec massRatios, double r){
  // for(int i = 0; i < iso.n_rows; i++){
  //   iso.row(i).operator*=(massRatios.at(i));
  // }
  // return iso;
  arma::sp_mat::const_iterator xStart = x.begin();
  arma::sp_mat::const_iterator xEnd   = x.end();
  NumericVector row;
  NumericVector column;
  for(arma::sp_mat::const_iterator it = xStart; it != xEnd; ++it)
  {
    double value = (*it);
    if(value != 0){
      int rowNum = it.row();
      double ratio = value * massRatios.at(rowNum);
      if(ratio >= r && ratio <= 1/r){
        int colNum = it.col();
        row.push_back(rowNum);
        column.push_back(colNum);
      }
    }
  }
  return List::create(Rcpp::Named("rowIdx") = row,
                      Rcpp::Named("colIdx") = column);
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
arma::sp_mat GibbsSampling_Int_NoPot(   List removal,
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
                                        bool log = false
){
  // CALLGRIND_START_INSTRUMENTATION;
  // CALLGRIND_TOGGLE_COLLECT;
  int massNum = prior_spMat.n_rows;
  int fomuNum = prior_spMat.n_cols;
  NumericVector seqMass = initVec(massNum);
  NumericVector seqFomu = initVec(fomuNum);
  // arma::mat add_mat = as<arma::mat>(add);
  // arma::mat iso_mat = as<arma::mat>(iso);
  // arma::mat bio_mat = as<arma::mat>(bio);
  // arma::mat prior_mat = as<arma::mat>(prior);
  arma::uvec seqMass_uvec = as<arma::uvec>(seqMass);
  arma::rowvec massInt_rvec = as<arma::rowvec>(massInt);
  arma::mat allSampFomu_mat = arma::zeros<arma::mat>(itNum,massNum);
  NumericVector sampFomu;
  for(int i = 0; i < massNum; i++){
    arma::sp_mat priorRow_sp = prior_spMat.row(i);
    arma::rowvec priorRow_rvec = spVectorGetDense(priorRow_sp);
    NumericVector priorRowSample = sampleRcppUnExport(seqFomu, 1, true, wrap(priorRow_rvec));
    sampFomu.push_back(priorRowSample.at(0));
  }
  sampFomu = sampFomu - 1;
  arma::rowvec sampFomu_rvec = as<arma::rowvec>(sampFomu);
  arma::sp_mat potBio_sp = spMatColSum(bio_spMat, sampFomu_rvec);
  NumericVector ordine(massNum);
  int v = 0;
  for(int i = 0; i < itNum; i++){
    ordine = sampleRcppUnExport(seqMass, massNum, false);
    ordine = ordine - 1;
    for(auto n : ordine ){
      int thism = n;
      // assignmented compounds of masses retained
      NumericVector remMass = removal.at(thism);
      arma::uvec remMass_uvec = as<arma::uvec>(remMass);
      arma::uvec retainMass_uvec = setDiff(seqMass_uvec,remMass_uvec);
      arma::rowvec retainFomu_rvec;
      if(!retainMass_uvec.is_empty()){
        retainMass_uvec = retainMass_uvec - 1;
        retainFomu_rvec = sampFomu_rvec.elem(retainMass_uvec).t();
      }
      
      // counting adductrelations
      arma::sp_mat pAdd_sp = spMatColSum(add_spMat, retainFomu_rvec);

      // counting isotoperelations
      // in R:
      // tmp <- matrix(Iso[sampcomp[-ind.rem[[thism]]],],ncol = Nc)*(Int[thism]/Int[-ind.rem[[thism]]])
      //   ind.ones <- which((tmp>=ratio.toll) & (tmp <=(1/ratio.toll))) 
      //   tmp[ind.ones]<-1
      // tmp[tmp!=1] <-0
      // p.iso<-colSums(tmp)
      arma::sp_mat isoSub_spMat = spMatSubView(iso_spMat,retainFomu_rvec);
      arma::rowvec massIntRatio_rvec = massInt_rvec.elem(retainMass_uvec).t();
      List l = thresholdIntRatio(isoSub_spMat, massIntRatio_rvec, ratioToll);
      arma::rowvec rowIdx_rvec = l["rowIdx"];
      arma::rowvec colIdx_rvec = l["colIdx"];
      arma::rowvec colIdxUni_rvec = arma::unique(colIdx_rvec);
      arma::rowvec pIso_rvec = arma::zeros<arma::rowvec>(fomuNum);
      for (int k = 0; k < colIdxUni_rvec.n_elem; k++){
        int colIdx = colIdxUni_rvec.at(k);
        arma::uvec colFind_uvec = arma::find(colIdx_rvec == colIdx);
        int num = colFind_uvec.n_elem;
        pIso_rvec.at(colIdx) = num;
      }

      // counting biotransformations relations
      int remFomu_bio = sampFomu_rvec.at(thism);
      arma::sp_mat remBio_sp = bio_spMat.row(remFomu_bio);
      arma::sp_mat pBio_sp = potBio_sp - remBio_sp;

      // adding penalities, ##only into add relation and iso relation
      
      // normalising with deltas
      arma::rowvec pAdd_rvec = spVectorGetDense(pAdd_sp);
      arma::rowvec pBio_rvec = spVectorGetDense(pBio_sp);
      pAdd_rvec = (pAdd_rvec + delAdd) / sum(pAdd_rvec + delAdd);
      pIso_rvec = (pIso_rvec + delIso) / sum(pIso_rvec + delIso);
      pBio_rvec = (pBio_rvec + delBio) / sum(pBio_rvec + delBio);
      arma::rowvec prior_rvec  = spVectorGetDense(prior_spMat.row(thism));
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
        potBio_sp = potBio_sp - bio_spMat.row(oldSamp) + bio_spMat.row(newSamp);
      }
    }
    
    // when finishing each iteration 
    // get the whole new sampcomp vector, where each element(mess) get resampled using its posterior probability vector
    allSampFomu_mat.submat(i,0,i,massNum-1) = sampFomu_rvec;
    
    // print some logs
    if(log){
      int old_v = v;
      v = (i * 100) / itNum;
      if(old_v !=v){
        Rcout << "computing posterior probability: " << v << "%" <<std::endl;
      }
    } 
  }
  
  //calculate posterior probability using allsampcomp which is a sampling distribution matrix (row num: no.its; col num: m, which is the number of masses)
  arma::sp_mat result = computePost(itNum, massNum, fomuNum, burn, allSampFomu_mat);
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


