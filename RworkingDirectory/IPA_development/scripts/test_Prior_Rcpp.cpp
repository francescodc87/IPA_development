#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

NumericVector initVec(int num, double element)
{
  NumericVector result(num);
  for(int i = 1; i <= num; i++){
    result[i-1] = element;
  }
  return result;
}

NumericMatrix initMatrix(int row, int column, double num){
  NumericMatrix result(row,column);
  std::fill( result.begin(), result.end(), num ) ;
  return result;
}

NumericVector computeRTprior(){
  return R_NilValue;
}

CharacterVector getColumnNameVec (int x){
  CharacterVector result(x+1);
  for (int i = 0; i < x; i++){
    result[i] = toString(i+1);
  }
  result[x] = "unknown";
  return result;
}

// [[Rcpp::export]]
NumericMatrix test_computePriorRcpp(NumericVector mass,
                               NumericVector compMass,
                               NumericVector pk,  // prior knowledge vector
                               NumericVector RTs,
                               NumericVector ppm,
                               NumericVector compId,
                               List RTranges,     //List RTranges: A list, every element: RTmin, RTmax, RT_in, RT_out, is not applied here
                               double unknownProb,
                               double limit,
                               int Vit,
                               bool V = false,
                               bool RTrangeList = false
) {
  
  // get mass number massNum and compound number compNum 
  int massNum = mass.length();
  int compNum = compMass.length();
  
  // set priorKnowledge default value
  if(pk == R_NilValue){
    pk = initVec(compNum,1);
  }else{
    if(pk.length() != compNum){
      stop("Error: please check prior knowledge length");
    }
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
  if(ppm == R_NilValue){
    ppm = initVec(compNum,3);
  }else{
    if(ppm.length() != compNum){
      stop("Error: please check ppm length");
    }else{
      // if ppm has been inserted, examine ppm value
      for(int i = 0; i < ppm.length(); i++){
        if(ppm[i] <= 0){
          stop("Error: please check ppm value");
        }
      }
    }
  }
  
  // initiate the prior probability matrix
  arma::mat  prior_arma = arma::zeros<arma::mat>(massNum,compNum+1);
  
  // prior probability computing
  for(int i = 0; i < massNum; i++){
    // compute sigma vector(1xC)
    NumericVector sigmaVec = ppm * mass[i] / ( 2 * 1e6);
    // compute likelihood(1xC) di based on mass difference
    NumericVector part1 = (-0.5) / pow(sigmaVec,2);
    NumericVector part2 = pow((initVec(compNum,mass[i])-compMass),2);
    NumericVector di = exp((part1 * part2));
    //NumericVector testresult = exp((part1 * part2)) ;
    
    // // compute likelihood(1xC) rt based on RT 
    // // RT.prior <-rep(1,Nc)
    // //   if(length(RT.ranges)==Nc){
    // //     ind.RT.out <- which(!is.na(RT.ranges))
    // //     ind.RT.out <- ind.RT.out[which(RTs[i]<RT.min[ind.RT.out] | RTs[i]>RT.max[ind.RT.out])]
    // //     RT.prior[ind.RT.out] <- RT.prior.penality
    // //   }
    NumericVector rt;
    if(RTrangeList){ // test if use RTranges
      rt = initVec(compNum,1);
    }else{
      rt = initVec(compNum,1);
    }

    // prior probability ∝ mass difference likelihood * prior knowledge likelihood * RT prior likelihood
    // pr[i,1:Nc] <- (exp((-0.5*precision[i])*((compounds.mass-mass[i])^2)))*prior.knowledge*RT.prior
    arma::rowvec x = di * pk * rt;
    prior_arma.submat(i,0,i,compNum-1) = x;

    // if probability value is too small, set it to 0
    arma::rowvec  prior_rowi_arma = prior_arma.row(i);
    double rowi_sum = sum(prior_rowi_arma);
    if(!NumericVector::is_na(limit)){
      double thold = limit * rowi_sum;                  // threshold
      double replace = 0;                               //replace
      arma::uvec ids = find(prior_rowi_arma < thold);   // Find indices
      prior_rowi_arma.elem(ids).fill(replace);          // Assign value to condition
    }

    // compute unknown probability
    rowi_sum = sum(prior_rowi_arma);
    double u = (unknownProb/(1-unknownProb)) * rowi_sum;
    if(u==0){
      u = 1;
    }
    prior_rowi_arma[compNum] = u;

    //normalise the prior probability
    rowi_sum = sum(prior_rowi_arma);
    prior_rowi_arma = prior_rowi_arma / rowi_sum;
    prior_arma.row(i) = prior_rowi_arma;
    
    // print some logs
    if(V){
      if(i % Vit == 0) {
        Rcout << "Computing Prior in Rcpp, " << (i * 100) / (massNum - 1) <<"%" << std::endl;
      }
    } 
    
  }
  
  // wrap datatype back
  NumericMatrix prior = wrap(prior_arma);
  
  // // insert the column name of prior probability matrix
  // if(compId != R_NilValue){
  //   CharacterVector colName = getColumnNameVec(compNum);
  //   colnames(prior) = colName;
  // }
  
  return prior;
}



