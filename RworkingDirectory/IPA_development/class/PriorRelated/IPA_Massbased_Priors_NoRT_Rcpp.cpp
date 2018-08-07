#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]



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


// [[Rcpp::export]]
arma::sp_mat ComputePriorRcpp( NumericVector mass,
                               NumericVector compMass,
                               NumericVector pk,  // prior knowledge vector
                               NumericVector RTs,
                               double ppm,
                               List RTranges,     //List RTranges: A list, every element: RTmin, RTmax, RT_in, RT_out, is not applied here
                               double unknownProb,
                               double limit,
                               bool log = false,
                               bool RTrangeList = false
                               ) {
  // get mass number massNum and compound number compNum 
  int massNum = mass.length();
  int compNum = compMass.length();
  
  // set priorKnowledge default value
  if(pk.length() != compNum){
    stop("Error: please check prior knowledge length");
  }
  
  if(RTrangeList){ // test if use RTranges
    // initial 4 vectors related to RT prior computing
    NumericVector RT_min(compNum); 
    NumericVector RT_max(compNum);
    NumericVector RT_in(compNum);
    NumericVector RT_out(compNum);
    if(RTranges.length()==compNum){ // may happen nullexception
      // construct vector RT_min and RT_max, which are used to decide whther inside/outside RT window
      //+construct vector RT_in and RT_out, which are used in RT prior computing
      for(int i = 0; i < compNum; i++){
        NumericVector elemVec = RTranges[i];
        // element: RTmin, RTmax, RT_in, RT_out
        RT_min[i] = elemVec[0];
        RT_max[i] = elemVec[1];
        RT_in[i] = elemVec[2];
        RT_out[i] = elemVec[3];
      }
    }
    else{
      stop("Error: please check RTrange List length");
    }
  }

  // set ppm default value
  if(ppm == NA_REAL){
    Rcout << "NA detected!!!" <<std::endl;
    ppm = 5;
  }
  
  // initiate the prior probability matrix
  arma::sp_mat  prior_spMat = arma::sp_mat(massNum,compNum+1);
  
  int v = 0;
  // prior probability computing
  for(int i = 0; i < massNum; i++){
    // compute sigma vector(1xC)
    NumericVector ppms(compNum, ppm);
    NumericVector sigmaVec = ppms * mass[i] / ( 2 * 1e6);
    
    // compute likelihood(1xC) di based on mass difference
    NumericVector part1 = (-0.5) / pow(sigmaVec, 2);
    NumericVector masses(compNum, mass[i]);
    NumericVector part2 = pow((masses - compMass), 2);
    NumericVector di = exp((part1 * part2));
    // compute likelihood(1xC) rt based on RT 
    // RT.prior <-rep(1,Nc)
    //   if(length(RT.ranges)==Nc){
    //     ind.RT.out <- which(!is.na(RT.ranges))
    //     ind.RT.out <- ind.RT.out[which(RTs[i]<RT.min[ind.RT.out] | RTs[i]>RT.max[ind.RT.out])]
    //     RT.prior[ind.RT.out] <- RT.prior.penality
    //   }
    // NumericVector rt;
    // if(RTrangeList){ // test if use RTranges, here RT is not considered as a factor
    //   rt = initVec(compNum,1);
    // }else{
    //   rt = initVec(compNum,1);
    // }
    
    // prior probability ∝ mass difference likelihood * prior knowledge likelihood * RT prior likelihood
    // pr[i,1:Nc] <- (exp((-0.5*precision[i])*((compounds.mass-mass[i])^2)))*prior.knowledge*RT.prior
    // arma::rowvec x = di * pk * rt;
    arma::rowvec x = di * pk;
    prior_spMat.submat(i,0,i,compNum-1) = x;

    // if probability value is too small, set it to 0
    arma::rowvec  priorRow_rvec = spMatGetDenseRowVec(prior_spMat, i);
    double rowi_sum = sum(priorRow_rvec);
    if(!NumericVector::is_na(limit)){
      double thold = limit * rowi_sum;                  // threshold
      double replace = 0;                               //replace
      arma::uvec ids = find(priorRow_rvec < thold);   // Find indices
      priorRow_rvec.elem(ids).fill(replace);          // Assign value to condition
    }

    // compute unknown probability
    rowi_sum = sum(priorRow_rvec);
    double u = (unknownProb/(1-unknownProb)) * rowi_sum;
    if(u==0){
      u = 1;
    }
    priorRow_rvec[compNum] = u;
    
    //normalise the prior probability
    rowi_sum = sum(priorRow_rvec);
    priorRow_rvec = priorRow_rvec / rowi_sum;
    prior_spMat.row(i) = priorRow_rvec;

    // print some logs
    if(log){
      int old_v = v;
      v = (i * 100) / massNum;
      if(old_v !=v){
        Rcout << "computing prior probability: " << v << "%" <<std::endl;
      }
    } 
  }
  
  
  return prior_spMat;
}



