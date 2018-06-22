
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
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
  NumericVector result = RcppArmadillo::sample(x, size, replace, prob);
  // printOutVec (result);
  return result;
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

// // [[Rcpp::export]]
// arma::rowvec testColSumWithSelectedRows(NumericMatrix a, NumericVector b, NumericVector c){
//   arma::mat a_arma = as<arma::mat>(a); 
//   arma::uvec b_arma = as<arma::uvec>(b); 
//   arma::mat subMat_arma = a_arma.rows(b_arma); 
//   arma::rowvec rvec = sum( subMat_arma, 0 ); 
//   arma::rowvec c_arma = as<arma::rowvec>(c);
//   arma::rowvec result = c_arma - rvec;
//   return result;
// }

// // [[Rcpp::export]]
// arma::rowvec testRowVecMultiply(NumericVector a, NumericVector b){
//   arma::rowvec a_arma = as<arma::rowvec>(a);
//   arma::rowvec b_arma = as<arma::rowvec>(b);
//   arma::rowvec result = a_arma % b_arma;
//   return result;
// }

// // [[Rcpp::export]]
// NumericVector sampleRcppExport( NumericVector x,
//                                   int size,
//                                   bool replace,
//                                   NumericVector prob = NumericVector::create()
// )
// {
//   NumericVector result = RcppArmadillo::sample(x, size, replace, prob);
//   // printOutVec (result);
//   return result;
// }

// [[Rcpp::export]]
NumericMatrix GibbsSampling(        List removal,
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
                                    int burn,
                                    bool v
){
  bool replace = false;
  int massNum = prior.rows();
  int compNum = prior.cols();
  NumericVector seqMass = initVec(massNum);
  NumericVector seqComp = initVec(compNum);
  arma::mat add_arma = as<arma::mat>(add);
  arma::mat iso_arma = as<arma::mat>(iso);
  arma::mat bio_arma = as<arma::mat>(bio);
  arma::rowvec potAdd_arma = as<arma::rowvec>(potAdd);
  arma::rowvec potIso_arma = as<arma::rowvec>(potIso);
  arma::rowvec potBio_arma = as<arma::rowvec>(potBio);
  NumericMatrix allsampcomp(itNum, massNum);
  
  
  for(int i = 0; i < itNum; i++){
    NumericVector ordine(massNum);
    ordine = sampleRcppUnExport(seqMass, massNum, replace);
    
    for(auto n : ordine )    {
      int thism = (n - 1);
      
      // assignmented compounds of masses need to remove
      NumericVector remMass = listEleIntoVec(removal, thism) - 1;
      NumericVector remComp = sampcomp[remMass];
      remComp = remComp - 1;
      arma::uvec remcomp_arma = as<arma::uvec>(remComp);
      // counting adductrelations
      // NumericVector remAdd = getRowsAndDoColsum(add, remComp);
      // NumericVector pAdd = potAdd - remAdd;
      arma::rowvec remAdd_arma = sum(add_arma.rows(remcomp_arma),0);
      arma::rowvec pAdd_arma = potAdd_arma - remAdd_arma;
      
      
      // counting isotoperelations
      // NumericVector remIso = getRowsAndDoColsum(iso, remComp);
      // NumericVector pIso = potIso - remIso;
      arma::rowvec remIso_arma = sum(iso_arma.rows(remcomp_arma),0);
      arma::rowvec pIso_arma = potIso_arma - remIso_arma;
      
      // counting biotransformations relations
      // NumericVector remBio = bio.row(sampcomp[thism] - 1);
      // NumericVector pBio = potBio - remBio;
      NumericVector remComp_bio = sampcomp[thism] - 1;
      arma::uvec remcomp_bio_arma = as<arma::uvec>(remComp_bio);
      arma::rowvec remBio_arma = sum(bio_arma.rows(remcomp_bio_arma),0);
      arma::rowvec pBio_arma = potBio_arma - remBio_arma;
      
      // adding penalities, ##only into add relation and iso relation
      
      // normalising with deltas
      pAdd_arma = (pAdd_arma + delAdd) / sum(pAdd_arma + delAdd);
      pIso_arma = (pIso_arma + delIso) / sum(pIso_arma + delIso);
      pBio_arma = (pBio_arma + delBio) / sum(pBio_arma + delBio);
      
      // merging scores, dot product
      NumericVector priorRow = prior.row(thism);
      arma::rowvec priorRow_arma = as<arma::rowvec>(priorRow);
      arma::rowvec post_arma = priorRow_arma % pAdd_arma % pIso_arma % pBio_arma;
      
      // // eliminate negative value
      // post[post < 0] = 0;
      
      // normalise posterior probability
      arma::rowvec posterior_arma = post_arma / sum(post_arma);
      
      // get the origin sampling result of mass 'thism'
      int oldsamp = sampcomp[thism];
      
      // use normalised probability po to re-sample mass 'thism'
      
      
      // NumericVector newSampVec = multsampleUnExport(sampVec, true, posterior);
      NumericVector posterior = wrap(posterior_arma);
      NumericVector newSampVec = sampleRcppUnExport(seqComp, 1, true, posterior);
      int newSamp = newSampVec[0];
      
      // update sampcomp
      sampcomp[thism] = newSamp;
      
      // if re-sample works, then
      if(oldsamp!=newSamp){
        potAdd_arma = potAdd_arma - add_arma.row(oldsamp-1) + add_arma.row(newSamp-1);
        potIso_arma = potIso_arma - iso_arma.row(oldsamp-1) + iso_arma.row(newSamp-1);
        potBio_arma = potBio_arma - bio_arma.row(oldsamp-1) + bio_arma.row(newSamp-1);
      }
    }
    
    // when finishing each iteration 
    // get the whole new sampcomp vector, where each element(mess) get resampled using its posterior probability vector
    allsampcomp.row(i) = sampcomp;
    if(v){
      Rcout << "Computing Posterior in Rcpp, " << (i * 100) / itNum <<"%" << std::endl;
    }
  }
  
  //calculate posterior probability using allsampcomp which is a sampling distribution matrix (row num: no.its; col num: m, which is the number of masses)
  NumericMatrix result = computePost(itNum, massNum, compNum, burn, allsampcomp);
  
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

