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

NumericVector initVec(int num)
{
  NumericVector result(num);
  for(int i = 1; i <= num; i++){
    result[i-1] = i;
  }
  return result;
}

// [[Rcpp::export]]
arma::sp_mat ComputePriorRcppThesis( NumericVector mass,
                                     NumericVector compMass,
                                     NumericVector pk,  // prior knowledge vector
                                     NumericVector RTs,
                                     double ppm,
                                     List RTranges,     //List RTranges: A list, every element: RTmin, RTmax, RT_in, RT_out, is not applied here
                                     int u,
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
  
  // build ppm vector
  NumericVector ppms(compNum, ppm);
  arma::rowvec ppms_rvec = as<arma::rowvec>(ppms);
  
  // alter some data type
  arma::rowvec mass_rvec = as<arma::rowvec>(mass);
  arma::rowvec compMass_rvec = as<arma::rowvec>(compMass);
  arma::rowvec pk_rvec = as<arma::rowvec>(pk);
  
  NumericVector compSeq = initVec(compNum) - 1;
  arma::uvec compSeq_uvec = as<arma::uvec>(compSeq);
  
  int v = 0;
  // prior probability computing
  for(int i = 0; i < massNum; i++){
    // compute sigma vector(1xC)
    arma::rowvec prior_rvec(compNum + 1);
    prior_rvec.fill(1);
    double massValue = mass_rvec.at(i);
    arma::rowvec sigmaVec_rvec = (ppms_rvec * massValue) / ( 2 * 1e6);
    for (int j = 0; j < compNum; j++){
      double compMassValue = compMass_rvec.at(j);
      double sigma = sigmaVec_rvec.at(j);
      if ((compMassValue >= massValue - u * sigma) && (compMassValue <= massValue + u * sigma )){
        prior_rvec.at(j) = 0;
      }
    }
    
    if (arma::sum(prior_rvec) == compNum + 1){
      prior_spMat.at(i, compNum) = 1;
    }else{
      // compute likelihood(1xC) di based on mass difference
      arma::rowvec part1_rvec = (-0.5) / arma::pow(sigmaVec_rvec, 2);
      NumericVector masses(compNum, mass[i]);
      NumericVector part2 = pow((masses - compMass), 2);
      arma::rowvec part2_rvec = as<arma::rowvec>(part2);
      arma::rowvec di_rvec = exp((part1_rvec % part2_rvec));
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
      arma::rowvec x_rvec = di_rvec % pk_rvec;
      arma::rowvec xBackUp_rvec = x_rvec;
      // if probability value is too small, set it to 0
      double rowi_sum = arma::sum(x_rvec);
      if(!NumericVector::is_na(limit)){
        double thold = limit * rowi_sum;                  // threshold
        double replace = 0;                               //replace
        arma::uvec ids = find(x_rvec < thold);   // Find indices
        x_rvec.elem(ids).fill(replace);          // Assign value to condition
      }
      // prior_spMat.submat(i,0,i,compNum-1) = x_rvec;
      
      if(arma::sum(x_rvec) == 0){
        xBackUp_rvec = xBackUp_rvec / arma::sum(xBackUp_rvec);
        prior_spMat.submat(i,0,i,compNum-1) = xBackUp_rvec;
      } else{
        x_rvec = x_rvec / arma::sum(x_rvec);
        prior_spMat.submat(i,0,i,compNum-1) = x_rvec;
      }
    }
    
    
    
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



