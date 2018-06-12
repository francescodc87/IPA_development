
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



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

// void printOutVec(NumericVector x){
//   for(int q = 0; q< x.size(); q++){
//     std::cout << "index  "<< q <<std::endl;
//     std::cout << "element  "<< x[q] <<std::endl;
//   }
// }


NumericVector sampleRcppUnExport( NumericVector x,
                                  int size,
                                  bool replace,
                                  NumericVector prob = NumericVector::create()
)
{
  int vecSize = x.size();
  if (vecSize == 1){
    NumericVector y = initVec(x[0]);
    NumericVector ret = RcppArmadillo::sample(y, size, replace, prob);
    // printOutVec (ret);
    return ret;
  }else{
    NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
    // printOutVec (ret);
    return ret;
  }
}


NumericMatrix computePost(int itNum, 
                          int massNum, 
                          int compNum, 
                          int burn, 
                          NumericMatrix allsampcomp){
  int valid = itNum - burn;
  NumericMatrix posterior(massNum, compNum);
  for (int i = 0; i < massNum; i++){
    NumericVector massSampCountVec(compNum);
    NumericVector col = allsampcomp.column(i);
    NumericVector massSampVec = col[Range(burn, itNum - 1)];
    for (int o = 0; o < valid; o++){
      int sampCompIdx = massSampVec[o] - 1;
      massSampCountVec[sampCompIdx] +=  1;
    }
    NumericVector posteriorPerMass = massSampCountVec / valid;
    posterior.row(i) = posteriorPerMass;
  }
  
  return posterior;

}

NumericVector listEleIntoVec(List ls, int idx){
  return ls[idx];
}


// [[Rcpp::export]]
NumericMatrix gibbsSampling(  List removal,
                              NumericVector sampcomp,
                              NumericVector potAdd,
                              NumericVector potIso,
                              NumericVector potBio,
                              NumericMatrix add,
                              NumericMatrix iso,
                              NumericMatrix bio,
                              NumericMatrix prior,
                              double delAdd, 
                              double delIso,
                              double delBio,
                              int itNum,
                              int burn
                              ){
    bool replace = false;
    int massNum = prior.rows();
    int compNum = prior.cols();
    // std::cout << "massNum!!!! "<< massNum <<std::endl;
    // std::cout << "compNUm!!!! "<< compNum <<std::endl;
    NumericMatrix allsampcomp(itNum, massNum);
    // NumericMatrix allposterior(massNum, compNum);
    for(int i = 0; i < itNum; i++){
      // std::cout << "print sampcomp here"<< std::endl;
      // printOutVec(sampcomp);

      NumericVector ordine(massNum);
      NumericVector massNumVec(1);
      massNumVec[0] = massNum;
      ordine = sampleRcppUnExport (massNumVec, massNum, replace);
      // ordine = sampleUnExport(massNum,massNum,false);
      
      // int x = 0;
      // NumericVector result(size);
      
      // std::cout << "print ordine here "<< std::endl;
      // printOutVec(ordine);
      for(const auto n : ordine )    {
        int thism = (n - 1);
        // std::cout << "ttttttttttthismmmmmmmm "<< thism <<std::endl;
        NumericVector remMess = listEleIntoVec(removal, thism) - 1;
        // std::cout << "rrrrrremMesssss "<< thism <<std::endl;
        NumericVector remComp = sampcomp[remMess];
        // std::cout << "rrrrrremCompppp "<< thism <<std::endl;
        
        // counting adductrelations
        NumericVector remAdd = getRowsAndDoColsum(add, remComp);
        NumericVector pAdd = potAdd - remAdd;
        // std::cout << "padd done "<<std::endl;
        
        // counting isotoperelations
        NumericVector remIso = getRowsAndDoColsum(iso, remComp);
        NumericVector pIso = potIso - remIso;
        // std::cout << "piso done "<<std::endl;
        
        // counting biotransformations relations
        NumericVector remBio = bio.row(sampcomp[thism] - 1);
        NumericVector pBio = potBio - remBio;
        // std::cout << "pbio done "<<std::endl;
  
        // adding penalities, ##only into add relation and iso relation

        // normalising with deltas
        pAdd = (pAdd + delAdd) / vecSum(pAdd + delAdd);
        pIso = (pIso + delIso) / vecSum(pIso + delIso);
        pBio = (pBio + delBio) / vecSum(pBio + delBio);
        // std::cout << "normalised with deltas done "<<std::endl;
        
        // merging scores, dot product
        NumericVector post = prior.row(thism) * pAdd * pIso * pBio;
        // NumericVector post = prior.row(thism);
        // NumericVector post = pBio;
        // std::cout << "posterior unnormalised done "<< thism << std::endl;
        
        // normalise posterior probability
        NumericVector posterior = post / vecSum(post);
        // std::cout << "posterior normalised done "<< thism <<std::endl;
        
        // // get the origin sampling result of mass 'thism'
        int oldsamp = sampcomp[thism];
        // std::cout << "oldsamp done "<< "thism " << thism << "oldsamp: " << oldsamp <<std::endl;

        // use normalised probability po to re-sample mass 'thism'
        NumericVector sampVec = initVec(compNum);
        
        // NumericVector newSampVec = multsampleUnExport(sampVec, true, posterior);
        NumericVector newSampVec = sampleRcppUnExport(sampVec, 1, true, posterior);
        int newSamp = newSampVec[0];
        // std::cout << "newsamp done "<< "thism " << thism << "newsamp: " << newSamp <<std::endl;

        // update sampcomp
        sampcomp[thism] = newSamp;
        // std::cout << "sampcomp update done "<< thism <<std::endl;


        // if re-sample works, then
        if(oldsamp!=newSamp){
          potAdd = potAdd - add.row(oldsamp-1) + add.row(newSamp-1);
          potIso = potIso - iso.row(oldsamp-1) + iso.row(newSamp-1);
          potBio = potBio - bio.row(oldsamp-1) + bio.row(newSamp-1);
          // std::cout << "colsum update done "<< thism <<std::endl;
        }
        // std::cout << "no colsum update "<< thism <<std::endl;
    }
      
    // when finishing each iteration 
    // get the whole new sampcomp vector, where each element(mess) get resampled using its po (posterior probability vector)
    allsampcomp.row(i) = sampcomp;
  }
    // std::cout << "allsampcomp update done "<<std::endl;
    
  // calculate posterior probability using allsampcomp which is a sampling distribution matrix (row num: no.its; col num: m, which is the number of masses)
  NumericMatrix posterior = computePost(itNum, massNum, compNum, burn, allsampcomp);
  
  // decide which kind of value to return
  // if(allsamp){
  //   List outList;
  //   outList["post"] = posterior;
  //   outList["allsampcomp"] = sampcomp;
  //   return outList;
  // }else{
    return posterior;  
  // }
}

// // [[Rcpp::export]]
// NumericVector testVecProd(NumericVector X){
//   bool replace = FALSE;
//   const NumericVector post = NumericVector::create(0.8,0.1,0.1);
//   NumericVector newsamp = sampleRcppUnExport(X, 1, replace, post);
//   for (auto n : post){
//     std::cout << n;
//   }
//   return newsamp;
//   
// }
// 
// // [[Rcpp::export]]
// NumericVector testMatrixRow(NumericMatrix X, int Y){
//   NumericVector result = X.row(Y);
//   
//   return result;
//   
// }
// 
// // [[Rcpp::export]]
// NumericVector text(NumericMatrix mat, NumericVector rowIdices)
// {
//   NumericVector sum(mat.ncol());
//   for(auto n : rowIdices ){
//     sum += mat.row(n-1);
//   }
//   return sum;
// }
// //[[Rcpp::export]]
// NumericVector selectCols(NumericMatrix mat, int idx){
//   NumericVector result = mat.column(idx-1);
//   return result;
// }


// // [[Rcpp::export]]
// NumericVector testDivision(NumericVector a, int b)
// { 
//   return a / b ;
// }


// NumericVector sampleUnExport( int x,
//                               int size,
//                               bool replace
// )
// { 
//     // Obtain environment containing function
//     Rcpp::Environment package_env("package:base"); 
//     
//     // Make function callable from C++
//     Rcpp::Function sample = package_env["sample"]; 
//     // NumericVector ret = sample(Named("x",x), Named("size",size), Named("replace",replace),Named("prob",prob));
//     NumericVector ret = sample(x, size, replace);
//     return ret;
// }


// NumericVector multsampleUnExport( NumericVector x,
//                                   bool replace,
//                                   NumericVector prob
// )
// { 
//   // Obtain environment containing function
//   Rcpp::Environment package_env("package:base"); 
//   
//   // Make function callable from C++
//   Rcpp::Function sample = package_env["sample"]; 
//   // NumericVector ret = sample(Named("x",x), Named("size",size), Named("replace",replace),Named("prob",prob));
//   NumericVector ret = sample(x, 1, replace, prob);
//   return ret;
// }
