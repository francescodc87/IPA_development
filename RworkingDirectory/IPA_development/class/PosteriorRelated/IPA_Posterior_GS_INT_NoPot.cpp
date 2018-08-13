
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
// #include <Eigen/Sparse>
// #include "/usr/include/valgrind/callgrind.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// typedef Eigen::SparseMatrix<double, RowMajor> SpMat;

arma::rowvec spVectorGetDenseCol(arma::sp_mat x){
  arma::sp_mat::const_iterator xStart = x.begin();
  arma::sp_mat::const_iterator xEnd   = x.end();
  arma::rowvec result = arma::zeros<arma::rowvec>(x.n_rows);
  for(arma::sp_mat::const_iterator it = xStart; it != xEnd; ++it)
  {
    double p = (*it);
    if(p != 0){
      int rowNum = it.row();
      result.at(rowNum) = p;
    }
  }
  return result;
}

arma::sp_mat spMatRowSum(arma::sp_mat x, arma::rowvec cols){
  arma::sp_mat result = arma::sp_mat(x.n_rows, 1);
  for(int i = 0; i < cols.n_elem; i++){
    int colNum = cols.at(i);
    result += x.col(colNum);
  }
  return result;
}

arma::sp_mat spMatSubColView(arma::sp_mat x, arma::rowvec cols){
  arma::sp_mat result = arma::sp_mat(x.n_rows, cols.n_elem);
  for(int i = 0; i < cols.n_elem; i++){
    int colNum = cols.at(i);
    result.col(i) = x.col(colNum);
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

arma::rowvec thresholdIntRatio(arma::sp_mat x, arma::rowvec cols, arma::rowvec massRatios, double r){
  NumericVector row;
  NumericVector rowSum;
  // Rcout << "inside thresholdIntRatio" << std::endl;
  for(int i = 0; i < cols.n_elem; i++){
    // Rcout << "picked a column" << std::endl;
    int colNum = cols.at(i);
    // Rcout << "colNum: " << colNum << std::endl;
    arma::sp_mat::const_col_iterator xStart = x.begin_col(colNum);
    arma::sp_mat::const_col_iterator xEnd   = x.end_col(colNum);
    for(arma::sp_mat::const_col_iterator it = xStart; it != xEnd; ++it){
      // Rcout << "start sparese matrix column iteration" << std::endl;
      double value = (*it);
      // Rcout << "value: " << value << std::endl;
      if(value != 0){
        // Rcout << "massRatios number: " << massRatios.n_elem << std::endl; 
        // Rcout << "massRatios.at(i): " << massRatios.at(i) << std::endl;
        double ratio = value * massRatios.at(i);
        // Rcout << "ratio: " << ratio << std::endl;
        if(ratio >= r && ratio <= 1/r){
          // Rcout << "ratio accepted" << std::endl;
          int rowNum = it.row();
          // Rcout << "row Num: " << rowNum << std::endl;
          arma::rowvec row_rvec = as<arma::rowvec>(row);
          // Rcout << "start findRow_uvec" << std::endl;
          arma::uvec findRow_uvec = arma::find(row_rvec == rowNum);
          if(findRow_uvec.is_empty()){
            // Rcout << "row not found!" << rowNum << std::endl;
            row.push_back(rowNum);
            rowSum.push_back(1);
          }else{
            // Rcout << "row found!" << rowNum << std::endl;
            int idx = findRow_uvec.at(0);
            rowSum.at(idx) += 1;
          }
        }
      }
    }
  }
  // Rcout << "start sum" << std::endl;
  arma::rowvec result = arma::zeros<arma::rowvec>(x.n_rows);
  for(int i = 0; i < row.length(); i++){
    int idx = row.at(i);
    result.at(idx) = rowSum.at(i);
  }
  // Rcout << "end sum" << std::endl;
  return result;
}

// void printOutUvec(arma::uvec x){
//   // // Rcout << "UVector" << std::endl;
//   x.t().raw_print(std::cout);
// }
// void printOutRvec(arma::rowvec x){
//   // // Rcout << "RVector" << std::endl;
//   x.raw_print(std::cout);
// }
// void printOutNvec(NumericVector x){
//   for(int i = 0; i < x.length(); i++){
//     // // Rcout << "NumericVector, " << x.at(i) << std::endl;
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
  prior_spMat = prior_spMat.t();
  add_spMat = add_spMat.t();
  iso_spMat = iso_spMat.t();
  bio_spMat = bio_spMat.t();
  // Rcout<< "start probability sample" << std::endl;
  for(int i = 0; i < massNum; i++){
    arma::sp_mat priorCol_sp = prior_spMat.col(i);
    arma::rowvec priorCol_rvec = spVectorGetDenseCol(priorCol_sp);
    NumericVector priorColSample = sampleRcppUnExport(seqFomu, 1, true, wrap(priorCol_rvec));
    sampFomu.push_back(priorColSample.at(0));
  }
  // Rcout<< "comlete probability sample" << std::endl;
  sampFomu = sampFomu - 1;
  arma::rowvec sampFomu_rvec = as<arma::rowvec>(sampFomu);
  arma::sp_mat potBio_sp = spMatRowSum(bio_spMat, sampFomu_rvec);
  NumericVector ordine(massNum);
  int v = 0;
  for(int i = 0; i < itNum; i++){
    // Rcout<< "enter into i iteration" << std::endl;
    ordine = sampleRcppUnExport(seqMass, massNum, false);
    ordine = ordine - 1;
    int j = 0;
    int w = 0;
    for(auto n : ordine ){
      // Rcout<< "enter into thism iteration" << std::endl;
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
      // Rcout<< "comlete retainFomu" << std::endl;
      // counting adductrelations
      arma::sp_mat pAdd_sp = spMatRowSum(add_spMat, retainFomu_rvec);
      // Rcout<< "comlete pAdd_sp" << std::endl;
      // counting isotoperelations
      // in R:
      // tmp <- matrix(Iso[sampcomp[-ind.rem[[thism]]],],ncol = Nc)*(Int[thism]/Int[-ind.rem[[thism]]])
      //   ind.ones <- which((tmp>=ratio.toll) & (tmp <=(1/ratio.toll))) 
      //   tmp[ind.ones]<-1
      // tmp[tmp!=1] <-0
      // p.iso<-colSums(tmp)
      // Rcout<< "start pIso_rvec" << std::endl;
      arma::rowvec massIntRatio_rvec =  massInt_rvec.elem(retainMass_uvec).t();
      // Rcout<< "start thresholdIntRatio" << std::endl;
      // Rcout<< "massInt_rvec.at(thism): " << massInt_rvec.at(thism) << std::endl;
      arma::rowvec pIso_rvec = thresholdIntRatio(iso_spMat, retainFomu_rvec, massInt_rvec.at(thism) / massIntRatio_rvec, ratioToll);
      // Rcout<< "comlete pIso_rvec" << std::endl;

      // counting biotransformations relations
      int remFomu_bio = sampFomu_rvec.at(thism);
      arma::sp_mat remBio_sp = bio_spMat.col(remFomu_bio);
      arma::sp_mat pBio_sp = potBio_sp - remBio_sp;
      // Rcout<< "comlete pBio_sp" << std::endl;
      // adding penalities, ##only into add relation and iso relation
      
      // normalising with deltas
      arma::rowvec pAdd_rvec = spVectorGetDenseCol(pAdd_sp);
      arma::rowvec pBio_rvec = spVectorGetDenseCol(pBio_sp);
      pAdd_rvec = (pAdd_rvec + delAdd) / sum(pAdd_rvec + delAdd);
      pIso_rvec = (pIso_rvec + delIso) / sum(pIso_rvec + delIso);
      pBio_rvec = (pBio_rvec + delBio) / sum(pBio_rvec + delBio);
      arma::rowvec prior_rvec  = spVectorGetDenseCol(prior_spMat.col(thism));
      arma::rowvec post_rvec = prior_rvec % pAdd_rvec % pIso_rvec % pBio_rvec;
      // Rcout<< "comlete post_rvec" << std::endl;
      
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
        potBio_sp = potBio_sp - bio_spMat.col(oldSamp) + bio_spMat.col(newSamp);
      }
      if(log){
        int old_w = w;
        w = (j * 100) / massNum;
        if(old_w != w){
          Rcout << "thism iteration: " << w << "%" <<std::endl;
          if(w % 10 == 0){
            Rcout<< "\014" << std::endl;
          }
        }
      } 
      j += 1;
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
        // if(v % 20 == 0){
        //   Rcout<< "\014" << std::endl;
        // }
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


