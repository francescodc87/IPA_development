
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
NumericMatrix GibbsSampling(  List removal,
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

    NumericMatrix allsampcomp(itNum, massNum);

    for(int i = 0; i < itNum; i++){
      NumericVector ordine(massNum);
      NumericVector massNumVec(1);
      massNumVec[0] = massNum;
      ordine = sampleRcppUnExport (massNumVec, massNum, replace);

      for(const auto n : ordine )    {
        int thism = (n - 1);
        
        // assignmented compounds of masses need to remove
        NumericVector remMass = listEleIntoVec(removal, thism) - 1;
        NumericVector remComp = sampcomp[remMass];
        
        // counting adductrelations
        NumericVector remAdd = getRowsAndDoColsum(add, remComp);
        NumericVector pAdd = potAdd - remAdd;
        
        // counting isotoperelations
        NumericVector remIso = getRowsAndDoColsum(iso, remComp);
        NumericVector pIso = potIso - remIso;
        
        // counting biotransformations relations
        NumericVector remBio = bio.row(sampcomp[thism] - 1);
        NumericVector pBio = potBio - remBio;
  
        // adding penalities, ##only into add relation and iso relation

        // normalising with deltas
        pAdd = (pAdd + delAdd) / vecSum(pAdd + delAdd);
        pIso = (pIso + delIso) / vecSum(pIso + delIso);
        pBio = (pBio + delBio) / vecSum(pBio + delBio);
        
        // merging scores, dot product
        NumericVector post = prior.row(thism) * pAdd * pIso * pBio;
        
        // // eliminate negative value
        // post[post < 0] = 0;
        
        // normalise posterior probability
        NumericVector posterior = post / vecSum(post);
        
        // // get the origin sampling result of mass 'thism'
        int oldsamp = sampcomp[thism];

        // use normalised probability po to re-sample mass 'thism'
        NumericVector sampVec = initVec(compNum);
        
        // NumericVector newSampVec = multsampleUnExport(sampVec, true, posterior);
        NumericVector newSampVec = sampleRcppUnExport(sampVec, 1, true, posterior);
        int newSamp = newSampVec[0];

        // update sampcomp
        sampcomp[thism] = newSamp;

        // if re-sample works, then
        if(oldsamp!=newSamp){
          potAdd = potAdd - add.row(oldsamp-1) + add.row(newSamp-1);
          potIso = potIso - iso.row(oldsamp-1) + iso.row(newSamp-1);
          potBio = potBio - bio.row(oldsamp-1) + bio.row(newSamp-1);
        }
    }
      
    // when finishing each iteration 
    // get the whole new sampcomp vector, where each element(mess) get resampled using its posterior probability vector
    allsampcomp.row(i) = sampcomp;
      if(v){
        Rcout << "Computing Posterior in Rcpp, " << (i * 100) / itNum <<"%" << std::endl;
      }
  }
    
  // calculate posterior probability using allsampcomp which is a sampling distribution matrix (row num: no.its; col num: m, which is the number of masses)
  NumericMatrix posterior = computePost(itNum, massNum, compNum, burn, allsampcomp);
  
  // decide which kind of value to return, temporarily ignored
  // if(allsamp){
  //   List outList;
  //   outList["post"] = posterior;
  //   outList["allsampcomp"] = sampcomp;
  //   return outList;
  // }else{
    return posterior;  
  // }
}
