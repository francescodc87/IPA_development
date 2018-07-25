#include <RcppArmadillo.h>
using namespace Rcpp;
typedef ListOf<NumericVector> NVList;
// [[Rcpp::depends(RcppArmadillo)]]

List UpdateMassAccuracyInfo(List ma_mean_std, double obs_mean, double obs_std, int obs_num){
  int old_num = ma_mean_std["ma_num"];
  double old_mean = ma_mean_std["ma_mean"];
  double old_std = ma_mean_std["ma_std"];
  
  // use basic statistical knowledge, calculate new mean and std
  int new_num = old_num + obs_num;
  double new_mean = (old_mean * old_num + obs_mean * obs_num) / (old_num + obs_num);
  
  double x = ((old_num - 1) * pow(old_std,2) + (obs_num - 1) * pow(obs_std,2) + old_num * pow(old_mean,2) + obs_num * pow(obs_mean,2) - new_num * pow(new_mean,2)) / (new_num - 1);
  double new_std = sqrt(x);
  ma_mean_std["ma_num"] = new_num;
  ma_mean_std["ma_mean"] = new_mean;
  ma_mean_std["ma_std"] = new_std;
  return ma_mean_std;
}

List UpdateIsoInfo(List IR_mean_std, double obs_mean, double obs_std, int obs_num){
  int old_num = IR_mean_std["iso_num"];
  double old_mean = IR_mean_std["iso_mean"];
  double old_std = IR_mean_std["iso_std"];
  
  // use basic statistical knowledge, calculate new mean and std
  int new_num = old_num + obs_num;
  double new_mean = (old_mean * old_num + obs_mean * obs_num) / (old_num + obs_num);
  
  double x = ((old_num - 1) * pow(old_std,2) + (obs_num - 1) * pow(obs_std,2) + old_num * pow(old_mean,2) + obs_num * pow(obs_mean,2) - new_num * pow(new_mean,2)) / (new_num - 1);
  double new_std = sqrt(x);
  IR_mean_std["iso_num"] = new_num;
  IR_mean_std["iso_mean"] = new_mean;
  IR_mean_std["iso_std"] = new_std;
  return IR_mean_std;
}

List UpdateRtInfo(List RT_mean_std, double obs_mean, double obs_std, int obs_num){
  int old_num = RT_mean_std["rt_num"];
  double old_mean = RT_mean_std["rt_mean"];
  double old_std = RT_mean_std["rt_std"];
  
  // use basic statistical knowledge, calculate new mean and std
  int new_num = old_num + obs_num;
  double new_mean = (old_mean * old_num + obs_mean * obs_num) / (old_num + obs_num);
  
  double x = ((old_num - 1) * pow(old_std,2) + (obs_num - 1) * pow(obs_std,2) + old_num * pow(old_mean,2) + obs_num * pow(obs_mean,2) - new_num * pow(new_mean,2)) / (new_num - 1);
  double new_std = sqrt(x);
  RT_mean_std["rt_num"] = new_num;
  RT_mean_std["rt_mean"] = new_mean;
  RT_mean_std["rt_std"] = new_std;
  return RT_mean_std;
}


double UpdatePriorKnowledge(double old_pk,
                            double scale){
  double new_pk = old_pk  + scale;
  return new_pk;
  
}

double UpdateMassAccuracy(double old_ma,
                          double obs_ma_mean,
                          double scale){
  double new_ma = old_ma + (obs_ma_mean - old_ma) * scale;
  return new_ma;
}

double UpdateIso(double old_iso,
                 double obs_iso_mean,
                 double scale){
  double new_iso = old_iso + (obs_iso_mean - old_iso) * scale;
  return new_iso;
}

void messageWhenDuplicateDetected(){
  Rcout << "-----------------------------------------------------!!!WARNING!!! Overlapped assignments detected!!!-----------------------------------------" << std::endl;
  Rcout << "Here are all information, please look !!!carefully!!! and then choose the correct assignment by type the mass number" << std::endl;
  Rcout << "Once the correct assignment is choosed, other assignment related to this compound will be !!!deleted!!!" << std::endl;
  Rcout << "If you can not decide which one is the correct assignment, type 0, then all assignments related to this compound will be !!!deleted!!!" << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
}

int messageAndChoice(int compNum, arma::rowvec massNum, arma::rowvec mz, arma::rowvec rt){
  Rcout << "-----------------------------------------------------------!!!Folloing are peak information!!!------------------------------------------------" << std::endl;
  Rcout << "Compound number: " << compNum + 1 << std::endl;
  for(int i = 0; i < massNum.n_elem; i++){
    Rcout << "Mass number: " << massNum.at(i) + 1 << ";  m/z: " << mz.at(i) << ";  RT: " << rt.at(i) << std::endl;
  }
  Rcout << "------------------------------------------------!!!Please insert your choice (mass number) in the next line!!!--------------------------------" << std::endl;
  Environment base = Environment("package:base");
  Function readline = base["readline"];
  Function as_numeric = base["as.numeric"];
  int result = as<int>(as_numeric(readline("> ")));
  Rcout << "-------------------------------------------------------------!!!Choice recorded!!!----------------------------------------------------------" << std::endl;
  return result - 1;
}

void messageAfterChoice(){
  Rcout << "--------------------------------------------------------!!!Successfully adjust assignments!!!-------------------------------------------------" << std::endl;
  Rcout << "---------------------------------If you don't see any more print information, then there is no more overlapped assignment;--------------------" << std::endl;
  Rcout << "----------------- ---------------------------Otherwise, !!!please keep solving overlapped assignments!!!--------------------------------------" << std::endl;
}

arma::uvec SetDiff(arma::uvec x, arma::uvec y){
  arma::uvec& result = x;
  
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uvec q1 = arma::find(x == y[j]);
    if (!q1.empty()) {
      result.shed_row(q1(0));
    }
  }
  
  return result;
}


void printOutRvec(arma::rowvec x){
  // Rcout << "UVector" << std::endl;
  for (int i = 0; i < x.n_elem; i++){
    Rcout << "idx: " << i << "; value: " << x.at(i) <<std::endl;
  }
}


template <typename D>
void testIterateFieldList(arma::field<D> x){
  int colNum = x.n_cols;
  int rowNum = x.n_rows;
  
  Rcout << "Row number: " << rowNum << "; Column number: " << colNum << std::endl;
}


arma::field<List> AdjustFieldListStructure(arma::field<List> x, int rowNum, int colNum){
  arma::field<List> result(rowNum,colNum);
  arma::field<List>::iterator start = result.begin();
  arma::field<List>::iterator end = result.end();
  for(arma::field<List>::iterator it = start; it != end; ++it)
  { 
    int idx = it.i;
    Rcout<< "idx: "<< idx << std::endl;
    (*it) = x(idx,0);
  }
  return result;
}

arma::field<NumericVector> AdjustFieldVectorStructure(arma::field<NumericVector> x, int rowNum, int colNum){
  arma::field<NumericVector> result(rowNum,colNum);
  arma::field<NumericVector>::iterator start = result.begin();
  arma::field<NumericVector>::iterator end = result.end();
  for(arma::field<NumericVector>::iterator it = start; it != end; ++it)
  { 
    int idx = it.i ;
    Rcout<< "idx: "<< idx << std::endl;
    (*it) = x(idx,0);
  }
  return result;
}

// // [[Rcpp::export]]
// List UpdateMain(NumericMatrix posterior,                // posterior matrix
//                 NumericVector mass,                     // mass values
//                 NumericVector compMass,                 // compound mass values
//                 NumericVector massInt,                  // intensity values of masses
//                 NumericVector pk,                       // prior knowledge
//                 NumericVector massMZ,                  // mass charge ratio of masses
//                 NumericVector massRT,                   // mass retention time
//                 double ma,                              // mass accuracy
//                 double s_pk,                            // update scale(speed) for pk
//                 double s_ma,                            // update scale for ma
//                 double s_iso,                           // update scale for iso
//                 double threshold                        // threshold used to filter posterior
//                   
// ){
// void printOutUvec(arma::uvec x){
//   // Rcout << "UVector" << std::endl;
//   for (int i = 0; i < x.n_elem; i++){
//     Rcout << "idx: " << i << "; value: " << x.at(i) <<std::endl;
//   }
// }


// // [[Rcpp::export]]
// NVList testListOf(NumericVector x){
//   NVList result;
//   result[0] = x;
//   return result;
// }

// // [[Rcpp::export]]
// arma::uvec testUvec(NumericVector x){
//   arma::uvec result;
//   result = as<arma::uvec>(x);
//   return result;
// }

// [[Rcpp::export]]
List UpdateMainLargeReal(arma::sp_mat post_spMat,                // posterior matrix
                arma::sp_mat iso_spMat,                 // iso matrix
                NumericVector mass,                     // mass values
                NumericVector compMass,                 // compound mass values
                NumericVector massInt,                  // intensity values of masses
                NumericVector pk,                       // prior knowledge
                NumericVector massMZ,                   // mass charge ratio of masses
                NumericVector massRT,                   // mass retention time
                NumericVector obsMA,                    // storage of observed mass accuracies in this update time
                NumericVector recMA,                    // record all observed mass accuracies
                List recMAdistr,                        // record all the distributions of mass accuracies
                List obsIR,            // storage of observed intensity ratios of all detected isotope pairs in this update time
                NumericVector obsIRlComp,               // left compound number of stored observed intensity ratio
                NumericVector obsIRrComp,               // right compound number of stored observed intensity ratio
                List recIRdistr,                // record all distributions of all isotope pairs
                NumericVector recIRdistrLComp,          // left compound number of all recorded distributions of all isotope pairs
                NumericVector recIRdistrRComp,          // right compound number of all recorded distributions of all isotope pairs
                List obsRT,            // storage of observed retention time in this update time
                List recRT,            // record all observed retention time
                List recRTdistr,                // record all the distributions of retention time
                NumericVector obsRTcomp,                // the compound numbers of the observed retention time in this update time
                NumericVector recRTdistrComp,           // the compound numbers of recorded retention time distributions
                double ma,                              // mass accuracy
                double s_pk,                            // update scale(speed) for pk
                double s_ma,                            // update scale for ma
                double s_iso,                           // update scale for iso
                double threshold                        // threshold used to filter posterior
                  
){
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>>>START THE UPDATE STAGE " << std::endl;
  
  NumericVector massSet;                                                  // mass number set after threshold
  NumericVector compSet;                                                  // compound number set after threshold
  
  // posterior probability matrix iterator creation
  arma::sp_mat::const_iterator post_start = post_spMat.begin();
  arma::sp_mat::const_iterator post_end   = post_spMat.end();
  
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>START TO FILTER POSTERIOR MATRIX" << std::endl;
  // filter the posterior probability matrix based on threshold, store qualified mass number, compound number at the same position respectively
  for(arma::sp_mat::const_iterator it = post_start; it != post_end; ++it)
  {
    int rowNum = it.row();
    int colNum = it.col();
    double p = (*it);
    if(p >= threshold){
      // storage step
      massSet.push_back(rowNum);
      compSet.push_back(colNum);
    }
  }
  Rcout << ">>>>>>>>>>>>>>>>>>>>>COMPLETE FILTERING POSTERIOR MATRIX" << std::endl;
  
  // check if overlap assignments happen
  Rcout << ">>>>>>>>>>>>>>>>>>>>>START TO CHECK OVERLAPPED ASSIGNMENT" << std::endl;
  arma::rowvec compSet_rvec = as<arma::rowvec>(compSet);
  arma::rowvec massSet_rvec = as<arma::rowvec>(massSet);
  arma::rowvec compSetUni_rvec = unique(compSet_rvec);
  if(compSetUni_rvec.n_elem != compSet_rvec.n_elem){
    arma::rowvec massMZ_rvec = as<arma::rowvec>(massMZ);
    arma::rowvec massRT_rvec = as<arma::rowvec>(massRT);
    
    // information indicated on the console
    messageWhenDuplicateDetected();
    
    // get duplicated compound number set
    NumericVector compSetUni = wrap(compSetUni_rvec);
    arma::uvec compSet_uvec = as<arma::uvec>(compSet);
    arma::uvec compSetUni_uvec = as<arma::uvec>(compSetUni);
    arma::uvec dupComp_uvec = SetDiff(compSet_uvec, compSetUni_uvec);
    arma::uvec dupCompUni_uvec = unique(dupComp_uvec);                      // that is the duplicated compound numbers
    
    // next, iterate on the duplicate compound number, extract mass numbers, and print helpful information
    for (int i = 0; i < dupCompUni_uvec.n_elem; i++){
      int dupCompNum = dupCompUni_uvec.at(i);                               // first get one of the duplicate compound number
      arma::uvec dupMass_uvec = find(compSet_rvec == dupCompNum);           // search the idx of the duplicate compound number in compSet_rvec
      arma::rowvec dupMass_rvec = massSet_rvec.elem(dupMass_uvec).t();      // get duplicate masses
      
      // provide information about duplication, and record user choice:
      int select_mass = -1;    // get the choice, next delete data accordingly
      if (select_mass == -1){                                                                    // if user cannot decide which mass to choose
        arma::uvec retainComp_uvec = find(compSet_rvec != dupCompNum);                           // delete this duplicated compound
        compSet_rvec = compSet_rvec.elem(retainComp_uvec).t();                                   // filtered compound set
        massSet_rvec = massSet_rvec.elem(retainComp_uvec).t();                                   // filtered mass set
        // messageAfterChoice();                                                                    // print some information after user making the choice each time
      } else{                                                                               // if user chooses a mass as the correct assignment
        arma::uvec allComp_uvec = find_finite(compSet_rvec);                                // all compound indices
        
        // get the intersection, which is the idx of the choosed mass
        arma::uvec selectComp_uvec = find(compSet_rvec == dupCompNum);
        arma::uvec selectMass_uvec = find(massSet_rvec == select_mass);
        arma::uvec retainIdx_uvec = arma::intersect(selectComp_uvec, selectMass_uvec);
        
        // get the indices of the removed masses (compounds) after choice
        arma::uvec remIdx_uvec = SetDiff(selectComp_uvec, retainIdx_uvec);
        
        // then get the diff between the removed idx and the whole indices, which result in the retain indices
        arma::uvec retainComp_uvec = SetDiff(allComp_uvec, remIdx_uvec);
        compSet_rvec = compSet_rvec.elem(retainComp_uvec).t();                                   // filtered compound set
        massSet_rvec = massSet_rvec.elem(retainComp_uvec).t();                                   // filtered mass set
        // messageAfterChoice();                                                                    // print some information after user making the choice each time
      }
    }
  }
  
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>OVERLAPPED ASSIGNMENTS ISSUE FIXED, CONTINUE THE UPDATE STAGE " << std::endl;
  
  // based on assignment pair after thresholding posterior probability matrix: 
  // which are: mass number and compound number,
  // do: update prior knowledge,
  //     record all observed mass accuracies in this update time before updating
  //     and record all observed intensity ratios of observed isotopes in this update time before updating
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START RECORD OBSERVED MASS ACCURACIES AND INTENSITY RATIOS, UPDATE PRIOR KNOWLEDGE " << std::endl;
  arma::rowvec pk_rvec = as<arma::rowvec>(pk);                                                          // prior knowledge vector
  arma::rowvec massInt_rvec = as<arma::rowvec>(massInt);                                                // intensity value vector of masses
  arma::rowvec obsRTComp_rvec = as<arma::rowvec>(obsRTcomp);                                            // compound numbers of observed retention times
  arma::rowvec obsIRlComp_rvec = as<arma::rowvec>(obsIRlComp);                                          // left part of compound numbers of observed intensity values
  arma::rowvec obsIRrComp_rvec = as<arma::rowvec>(obsIRrComp);                                          // right part of compound numbers of observed intensity values
  for(int i = 0; i < compSet_rvec.n_elem; i++){                                                         // iterate on thresholded compound numbers                     
    // after detect and fix overlap assignments, Update Prior Knowledge
    int theCompLeft = compSet_rvec.at(i);                                                               // first time loop, take the value as left part
    int theMassLeft = massSet_rvec.at(i);                                                               // first time loop, take the value as left part
    pk_rvec.at(theCompLeft) = UpdatePriorKnowledge(pk_rvec.at(theCompLeft), s_pk);
    
    // calculate and record observed mass accuracies
    double maValue = ((mass[theMassLeft] - compMass[theCompLeft]) * 1e6 ) / compMass[theCompLeft];
    obsMA.push_back(maValue);
    recMA.push_back(maValue);
    
    // record observed retention times
    double rtValue = massRT(theMassLeft);
    
    
    arma::uvec obsRTCompFind_uvec = find(obsRTComp_rvec == theCompLeft);                // check if we observed other retention time of this compound number before
    // if we have not observed before, push_back retention time value, and push_back the compound number 
    // else we pop the retention time value into the same position we find it, no need to pop_up the compound time cause we have already seen it in obsRTComp_rvec
    if(obsRTCompFind_uvec.is_empty()){  
      NumericVector rt;
      rt.push_back(rtValue);
      obsRT.push_back(rt);
      recRT.push_back(rt);
      obsRTcomp.push_back(theCompLeft);
    } else{
      int obsRTCompFind = obsRTCompFind_uvec.at(0,0);
      NumericVector obsRTsub = obsRT[obsRTCompFind];
      NumericVector recRTsub = recRT[obsRTCompFind];
      obsRTsub.push_back(rtValue);
      recRTsub.push_back(rtValue);
      obsRT[obsRTCompFind] = obsRTsub;
      recRT[obsRTCompFind] = recRTsub;
    }
    
    // next is the operation of recording the observed intensity ratios before updating, we use another iteration of compSet_rvec to create isotope pairs
    for(int j = 0; j < compSet_rvec.n_elem; j++){
      // generate the right part of compounds which will be used as indices to record intensity ratios
      // and generate the right part of masses which will be used as indices to calculate observed intensity ratios
      int theCompRight = compSet_rvec.at(j);
      int theMassRight = massSet_rvec.at(j);
      
      // next is to combine the compound number pair as isotope pair, use mass number pair to compute isotope ratio
      // and complete the storage task before update
      // about the next line judgement code: here we decrease the computational complexity 
      //                                     & we set the condition that only the compound number which is an isotope pair can be applied in the following task
      if((theCompLeft < theCompRight) && (iso_spMat.at(theCompLeft,theCompRight) != 0)){                
        double theIR_LR = massInt_rvec.at(theMassLeft) / massInt_rvec.at(theMassRight);                 // isotope ratio (Left part to right part)
        double theIR_RL = 1 / theIR_LR;                                                                 // isotope ratio (right part to left part)
        // we want to know if we have recorded any isotope ratios of this compound pair before
        // so we take the intersect
        arma::uvec obsIRlCompFind_uvec = find(obsIRlComp_rvec == theCompLeft);
        arma::uvec obsIRrCompFind_uvec = find(obsIRrComp_rvec == theCompRight);
        arma::uvec obsIRisoFind_uvec = arma::intersect(obsIRlCompFind_uvec, obsIRrCompFind_uvec);
        // if we have not recorded the observed isotope ratio of this compound number pair before, we push_back those value and indices
        // else we pop the value into the same position we find it, but we don't pop the indices in this situation cause the indices are already existed
        if(obsIRisoFind_uvec.is_empty()){                        
          obsIR.push_back(theIR_LR);
          obsIR.push_back(theIR_RL);
          obsIRlComp.push_back(theCompLeft);
          obsIRlComp.push_back(theCompRight);
          obsIRrComp.push_back(theCompRight);
          obsIRrComp.push_back(theCompLeft);
        } else{                                                                                 // compound number pair is already recorded
          int obsIRisoFind = obsIRisoFind_uvec.at(0,0);
          NumericVector obsIRlSub = obsIR[obsIRisoFind];
          NumericVector obsIRrSub = obsIR[obsIRisoFind + 1];
          obsIRlSub.push_back(theIR_LR);
          obsIRrSub.push_back(theIR_RL);
          obsIR[obsIRisoFind] = obsIRlSub;
          obsIR[obsIRisoFind + 1] = obsIRrSub;
        }
      }
    }
  }
  
  // now we have recorded observed mass accuracies and intensity ratios, it's time to update the mean and std values respectively
  
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO UPDATE MASS ACCURACY " << std::endl;
  // update mass accuracy
  arma::rowvec obsMA_rvec = as<arma::rowvec>(obsMA);
  if(recMAdistr.length() == 0){                                      // if we have not recorded any distribution of mass accuracy before
    if(obsMA_rvec.n_elem > 1){                                       // if observed more than one mass accuracies in this update time, we can pop at least the first distribution information of mass accuracy 
      double ma_mean = arma::mean(obsMA_rvec);
      double ma_std = arma::stddev(obsMA_rvec);
      double ma_num = obsMA_rvec.n_elem;
      recMAdistr["ma_num"] = ma_num;
      recMAdistr["ma_mean"] = ma_mean;
      recMAdistr["ma_std"] = ma_std;
      ma = UpdateMassAccuracy(ma, recMAdistr["ma_mean"], s_ma);      // update mass accuracy main calculation
      NumericVector tmp;                                           // clear observed mass accuracies
      obsMA = tmp;
    }
  } else{                                                                       // if we have recorded any distribution of mass accuracy before
    if(obsMA_rvec.n_elem > 0){                                                  // if have observed one mass accuracy value in this update time, then we can update the distribution of mass accuracy
      double ma_mean = arma::mean(obsMA_rvec);
      double ma_std = arma::stddev(obsMA_rvec);
      double ma_num = obsMA_rvec.n_elem;
      recMAdistr = UpdateMassAccuracyInfo(recMAdistr, ma_mean, ma_std, ma_num); // update mass accuracy distribution information
      ma = UpdateMassAccuracy(ma, recMAdistr["ma_mean"], s_ma);                 // update mass accuracy main calculation
      NumericVector tmp;                                                       // clear observed mass accuracies
      obsMA = tmp;
    }
  }
  
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO UPDATE ISOTOPE CONNECTIVITIES " << std::endl;
  // update iso observed intensity ratios, distribution information and iso connectivity matrix
  arma::rowvec obsIRlComp_rvec_new = as<arma::rowvec>(obsIRlComp);              // left part of the observed isotope pairs
  arma::rowvec obsIRrComp_rvec_new = as<arma::rowvec>(obsIRrComp);              // right part of the observed isotope pairs
  arma::rowvec recIRdistrLComp_rvec = as<arma::rowvec>(recIRdistrLComp);        // left part of the recorded isotope pairs (of the recorded distribution information)
  arma::rowvec recIRdistrRComp_rvec = as<arma::rowvec>(recIRdistrRComp);        // right part of the recorded isotope pairs (of the recorded distribution information)
  IntegerVector remIRidx;
  // iterate on all observed isotope pairs, 
  // the idea is to find out if the observed isotope pair is already found in the isotope pairs of recorded distribution information about intensity ratios
  for(int i = 0; i < obsIRlComp_rvec_new.n_elem; i = i + 2){              // as you can see here we decrease the computational complexity                
    int obsIRlCompNum = obsIRlComp_rvec_new.at(i);
    int obsIRrCompNum = obsIRrComp_rvec_new.at(i);
    NumericVector theIR_LR = obsIR[i];                                    // the intensity ratio_left part right part
    NumericVector theIR_RL = obsIR[i + 1];                                // the intensity ratio_right part left part
    arma::rowvec theIR_LR_rvec = as<arma::rowvec>(theIR_LR);
    arma::rowvec theIR_RL_rvec = as<arma::rowvec>(theIR_RL);
    arma::uvec recIRdistrLCompFind_uvec = find(recIRdistrLComp_rvec == obsIRlCompNum);
    arma::uvec recIRdistrRCompFind_uvec = find(recIRdistrRComp_rvec == obsIRrCompNum);
    arma::uvec recIRdistrIsoFind_uvec = arma::intersect(recIRdistrLCompFind_uvec, recIRdistrRCompFind_uvec);
    // if we have not recorded the distribution information of this observed compound number pair before, 
    // we can push_back the first distribution information of this observed compound number pair and the indices 
    // only if record more than 2 intensity ratios
    if(recIRdistrIsoFind_uvec.n_elem == 0){
      if(theIR_LR.length() > 1){
        List theIRdistr_LR;
        List theIRdistr_RL;
        // calculate IR numbers, mean and std value for this isotope pair for the first time
        double theIR_LR_mean = arma::mean(theIR_LR_rvec);
        double theIR_LR_std = arma::stddev(theIR_LR_rvec);
        double theIR_RL_mean = arma::mean(theIR_RL_rvec);
        double theIR_RL_std = arma::stddev(theIR_RL_rvec);
        double theIR_num = 2;
        theIRdistr_LR["iso_num"] = theIR_num;
        theIRdistr_LR["iso_mean"] = theIR_LR_mean;
        theIRdistr_LR["iso_std"] = theIR_LR_std;
        theIRdistr_RL["iso_num"] = theIR_num;
        theIRdistr_RL["iso_mean"] = theIR_RL_mean;
        theIRdistr_RL["iso_std"] = theIR_RL_std;
        recIRdistr.push_back(theIRdistr_LR);
        recIRdistr.push_back(theIRdistr_RL);
        recIRdistrLComp.push_back(obsIRlCompNum);
        recIRdistrLComp.push_back(obsIRrCompNum);
        recIRdistrRComp.push_back(obsIRrCompNum);
        recIRdistrRComp.push_back(obsIRlCompNum);
        
        // update iso matrix
        iso_spMat.at(obsIRlCompNum, obsIRrCompNum) = UpdateIso(iso_spMat.at(obsIRlCompNum, obsIRrCompNum), theIR_LR_mean, s_iso);
        iso_spMat.at(obsIRrCompNum, obsIRlCompNum) = UpdateIso(iso_spMat.at(obsIRrCompNum, obsIRlCompNum), theIR_RL_mean, s_iso);
        
        // record the data indices
        remIRidx.push_back(i);
        remIRidx.push_back(i + 1);
      } 
    } else{  
      // if we have recorded the distribution information of this observed compound number pair before,
      // we pop the distribution information into the same position we find it, 
      // but we don't pop the indices in this situation cause the indices are already existed
      int recIRdistrIsoFind = recIRdistrIsoFind_uvec.at(0,0);
      double theIR_LR_mean = theIR_LR_rvec.at(0);
      double theIR_RL_mean = theIR_RL_rvec.at(0);
      double theIR_std = 0;
      double theIR_num = 1;
      
      // update distribution information
      recIRdistr[recIRdistrIsoFind] = UpdateIsoInfo(recIRdistr[recIRdistrIsoFind], theIR_LR_mean, theIR_std, theIR_num);
      recIRdistr[recIRdistrIsoFind + 1] = UpdateIsoInfo(recIRdistr[recIRdistrIsoFind + 1], theIR_RL_mean, theIR_std, theIR_num);
      
      // update iso matrix
      iso_spMat.at(obsIRlCompNum, obsIRrCompNum) = UpdateIso(iso_spMat.at(obsIRlCompNum, obsIRrCompNum), theIR_LR_mean, s_iso);
      iso_spMat.at(obsIRrCompNum, obsIRlCompNum) = UpdateIso(iso_spMat.at(obsIRrCompNum, obsIRlCompNum), theIR_RL_mean, s_iso);
      remIRidx.push_back(i);
      remIRidx.push_back(i + 1);
    }
  }
  
  // after update iso, we delete some observed data from obsIR, obsIRlComp, obsIRrComp
  IntegerVector obsIRidx = seq_len(obsIRlComp.length()) - 1;
  arma::uvec obsIRidx_uvec = as<arma::uvec>(obsIRidx);
  arma::uvec remIRidx_uvec = as<arma::uvec>(remIRidx);
  arma::uvec retainIRidx_uvec = SetDiff(obsIRidx_uvec, remIRidx_uvec);
  IntegerVector retainIRidx = wrap(retainIRidx_uvec.t());
  obsIR = obsIR[retainIRidx];
  obsIRlComp = obsIRlComp[retainIRidx];
  obsIRrComp = obsIRrComp[retainIRidx];
  
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO RECORD RETENTION TIME DISTRIBUTION " << std::endl;
  arma::rowvec obsRTcomp_rvec_new = as<arma::rowvec>(obsRTcomp);                  // all compound numbers of observed retention time before update
  arma::rowvec recRTdistrComp_rvec = as<arma::rowvec>(recRTdistrComp);            // all compound numbers of recorded distribution information of retention time before update
  
  // iterate on all observed compound numbers of retention time, 
  // the idea is to find out if the observed compound is already found in the compound set of recorded distribution information about retention time
  IntegerVector remRTidx;
  for(int i = 0; i < obsRTcomp_rvec_new.n_elem; i++){
    int obsRTcompNum = obsRTcomp_rvec_new.at(i);
    arma::uvec recRTdistrCompFind_uvec = find(recRTdistrComp_rvec == obsRTcompNum);
    NumericVector obsRTsub = obsRT[i];
    // if we have not recorded the distribution information of this observed compoundbefore, 
    // we can push_back the first distribution information of this observed compound and the index (compound number) 
    // only if record more than 2 retention time values of this compound
    if(recRTdistrCompFind_uvec.n_elem == 0){
      if(obsRTsub.length() > 1){
        arma::rowvec obsRT_rvec = as<arma::rowvec>(obsRT[i]);
        double rt_mean = arma::mean(obsRT_rvec);
        double rt_std = arma::stddev(obsRT_rvec);
        double rt_num = 2;
        List rt_mean_std;
        rt_mean_std["rt_num"] = rt_num;
        rt_mean_std["rt_mean"] = rt_mean;
        rt_mean_std["rt_std"] = rt_std;
        recRTdistr.push_back(rt_mean_std); // here we push_back the distribution information
        recRTdistrComp.push_back(obsRTcompNum);      // here we push_back the index (compound number)     
        remRTidx.push_back(i);
      }
    }else{
      // if we have recorded the distribution information of this observed compound before,
      // we pop the distribution information into the same position we find it,
      // but we don't pop the indices in this situation cause the indices are already existed
      int recRTdistrCompFind = recRTdistrCompFind_uvec.at(0,0);
      double rt_mean = obsRTsub.at(0);
      double rt_std = 0;
      double rt_num = 1;
      
      // update retention time distribution information
      recRTdistr[recRTdistrCompFind] = UpdateRtInfo(recRTdistr[recRTdistrCompFind], rt_mean, rt_std, rt_num);
      remRTidx.push_back(i);
    }
  }
  
  // after record retention time distribution information, we remove those observed data from obsRT, obsRTcomp
  IntegerVector obsRTidx = seq_len(obsRTcomp.length()) - 1;
  arma::uvec obsRTidx_uvec = as<arma::uvec>(obsRTidx);
  arma::uvec remRTidx_uvec = as<arma::uvec>(remRTidx);
  arma::uvec retainRTidx_uvec = SetDiff(obsRTidx_uvec, remRTidx_uvec);
  IntegerVector retainRTidx = wrap(retainRTidx_uvec.t());
  obsRT = obsRT[retainRTidx];
  obsRTcomp = obsRTcomp[retainRTidx];
  
  
  
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>END OF UPDATE STAGE" << std::endl;
  
  // next is the return value set, what should we return?
  // pk: prior knowledge vector
  // ma: updated mass accuracy value, observed mass accuracies at this time, all time mass accuracy distributions, all time mass accuracy values 
  // iso: observed intensity ratios at this time + compound pair indices, plus observed intensity ratio numbers at this time 
  //      recorded distributions of intensity ratios of all thime + compound pair indices
  // rt: observed RT in this time, compound number of observed RT in this time, recorded RT of all time, recorded distributions of RT all time, recorded compound number of all thime recorded RT distributions
  
  // normalise prior knowledge
  pk_rvec = pk_rvec / sum(pk_rvec);
  
  // return a list, contain all the updated results
  return List::create(Rcpp::Named("pk") = wrap(pk_rvec),
                      Rcpp::Named("ma") = ma,
                      Rcpp::Named("iso") = iso_spMat,
                      Rcpp::Named("obsMA") = obsMA,
                      Rcpp::Named("recMA") = recMA,
                      Rcpp::Named("recMAdistr") = recMAdistr,
                      Rcpp::Named("recMA") = recMA,
                      Rcpp::Named("obsIR") = obsIR,
                      Rcpp::Named("obsIRlComp") = obsIRlComp,
                      Rcpp::Named("obsIRrComp") = obsIRrComp,
                      Rcpp::Named("recIRdistr") = recIRdistr,
                      Rcpp::Named("recIRdistrLComp") = recIRdistrLComp,
                      Rcpp::Named("recIRdistrRComp") = recIRdistrRComp,
                      Rcpp::Named("obsRT") = obsRT,
                      Rcpp::Named("recRT") = recRT,
                      Rcpp::Named("recRTdistr") = recRTdistr,
                      Rcpp::Named("obsRTcomp") = obsRTcomp,
                      Rcpp::Named("recRTdistrComp") = recRTdistrComp);
  
  // return R_NilValue;
}




// arma::mat UpdateIso(arma::mat iso_old,
//                     arma::rowvec compNum,
//                     arma::rowvec intValue,
//                     arma::rowvec idxRow,
//                     arma::rowvec idxCol,
//                     double scale){
//   
//   // iteration
//   arma::rowvec::const_iterator start_left = compNum.begin();
//   arma::rowvec::const_iterator end_left   = compNum.end();
//   arma::rowvec::const_iterator start_right = compNum.begin();
//   arma::rowvec::const_iterator end_right   = compNum.end();
//   int massIdxLeft = 0;
//   for(arma::rowvec::const_iterator it_left = start_left; it_left != end_left; ++it_left)
//   {
//     int compIdxLeft = (*it_left);
//     int massIdxRight = 0;
//     for (arma::rowvec::const_iterator it_right = start_right; it_right != end_right; ++it_right){
//       int compIdxRight = (*it_right);
//       double intRatio_old = iso_old.at(compIdxLeft, compIdxRight);                        // get the old intensity ratio in iso matrix
//       if((compIdxLeft != (*it_right)) && (intRatio_old != 0) && (!NumericVector::is_na(intRatio_old))){
//         // after get idx c1 (compIdxLeft), c2(compIdxRight), enter into mathematical content of the update
//         double intRatio_obj = intValue.at(massIdxLeft) / intValue.at(massIdxRight);       // compute new intensity ratio
//         double measure = fabs(intRatio_obj - intRatio_old) / (intRatio_obj * scale);
//         if(measure <= 1){
//           iso_old.at(compIdxLeft, compIdxRight) = intRatio_obj;
//         }else{
//           double intRatio_new = intRatio_old + (intRatio_obj - intRatio_old) / measure;
//           iso_old.at(compIdxLeft, compIdxRight) = intRatio_new;
//         }
//       }
//       massIdxRight += 1;
//     }
//     massIdxLeft += 1;
//   }
//   return iso_old;
// }

// // [[Rcpp::export]]                                                                        
// void testIterator(NumericMatrix x){
//   arma::mat x_mat = as<arma::mat>(x);
//   arma::sp_mat x_spmat = arma::sp_mat(x_mat);
//   arma::sp_mat::const_iterator start = x_spmat.begin();
//   arma::sp_mat::const_iterator end   = x_spmat.end();
//   
//   for(arma::sp_mat::const_iterator it = start; it != end; ++it)
//   {
//     Rcout << "location: " << it.row() << "," << it.col() << "  ";
//     Rcout << "value: " << (*it) << std::endl;
//   }
// } 

// // [[Rcpp::export]]   
// arma::rowvec testMatrix(NumericVector a){
//   arma::field<arma::rowvec> x(3, 3);
//   arma::rowvec a_rvec = as<arma::rowvec>(a);
//   x(0,0) = a;
//   return x(0,0);
// }
// 
// // [[Rcpp::export]]
// arma::field<List> testArmaFiel(arma::field<List> x){
//   arma::field<List> result = x;
//   result.at(0,0).push_back(1);
//   result.at(0,0).push_back(2);
//   result.at(0,0).push_back(3);
//   return result;
// }

// // [[Rcpp::export]]
// arma::field<NumericVector> testConvertToField(ListMatrix x, double b){
//   arma::field<NumericVector> result = as
// }

// // [[Rcpp::export]]
// double testArmaStd(arma::rowvec x){
//   return arma::stddev(x);
// }

// // [[Rcpp::export]]
// List testCalMeanStd(double old_mean, double old_std, int old_num, double obs_mean, double obs_std, int obs_num){
// 
//   // use basic statistical knowledge, calculate new mean and std
//   int new_num = old_num + obs_num;
//   double new_mean = (old_mean * old_num + obs_mean * obs_num) / (old_num + obs_num);
// 
//   // double x = (old_num * pow(old_std,2) + obs_num * pow(obs_std,2)) / (old_num + obs_num) + ( old_num * obs_num * pow((old_mean - obs_mean),2) ) / pow((old_num + obs_num),2);
//   double x = ((old_num - 1) * pow(old_std,2) + (obs_num - 1) * pow(obs_std,2) + old_num * pow(old_mean,2) + obs_num * pow(obs_mean,2) - new_num * pow(new_mean,2)) / (new_num - 1);
//   double new_std = sqrt(x);
//   return List::create(Rcpp::Named("num") = new_num,
//                       Rcpp::Named("mean") = new_mean,
//                       Rcpp::Named("std") = new_std);
// }
// 
// // [[Rcpp::export]]
// arma::mat testArmaElem(arma::mat x){
//   return x.elem(arma::find(x>8));
// }

// // [[Rcpp::export]]
// void testIteratorSpMat(arma::mat x){
//   arma::sp_mat x_spmat = arma::sp_mat(x);
//   arma::sp_mat::const_iterator start = x_spmat.begin();
//   arma::sp_mat::const_iterator end   = x_spmat.end();
//   for(arma::sp_mat::const_iterator it = start; it != end; ++it)
//   {
//     int rowNum = it.row();
//     int colNum = it.col();
//     Rcout << "column: "<< colNum << ", row: "<< rowNum << std::endl;
//   }
// }
// 
// // [[Rcpp::export]]
// arma::rowvec testFindDuplicates(arma::rowvec x){
//   NumericVector result;
//   for(int i = 0; i < x.n_elem; i++){
//     if((i < x.n_elem - 1) && (x.at(i) == x.at(i + 1))){
//       result.push_back(x.at(i));
//     }
//   }
//   arma::rowvec result_rvec = as<arma::rowvec>(result);
//   result_rvec = arma::unique(result_rvec);
//   return result_rvec;
// }
// 
// // [[Rcpp::export]]
// arma::uvec testArmaUnique(arma::uvec x){
//   return arma::unique(x);
// }
// 
// // [[Rcpp::export]]
// void testIterateUvec(arma::uvec x){
//   for(int i = 0; i < x.n_elem; i++){
//     Rcout << x.at(i)<< std::endl;
//   }
// }
// 
// 
// 
// // [[Rcpp::export]]
// void testWarningInfoAndConsoleRead(){
//   Rcout << "-----------------------------------------------------!!!WARNING!!! Overlapped assignments detected--------------------------------------------" << std::endl;
//   Rcout << "Here are all information, please look !!!carefully!!! and then choose the correct assignment by type the mass number" << std::endl;
//   Rcout << "Once the correct assignment is choosed, other assignment related to this compound will be !!!deleted!!!" << std::endl;
//   Rcout << "If you can not decide which one is the correct assignment, type 0, then all assignments related to this compound will be !!!deleted!!!" << std::endl;
//   Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
//   Environment base = Environment("package:base");
//   Function readline = base["readline"];
//   Function as_numeric = base["as.numeric"];
//   int input = as<int>(as_numeric(readline("> ")));
//   if(input == 0){
//     Rcout << "test" << std::endl;
//   }
// }
// 

// template <typename D>
// void testIterateFieldList(arma::field<D> x){
//   int colNum = x.n_cols;
//   int rowNum = x.n_rows;
//   
//   Rcout << "Row number: " << rowNum << "; Column number: " << colNum << std::endl;
// }


// arma::field<List> AdjustFieldListStructure(arma::field<List> x, int rowNum, int colNum){
//   arma::field<List> result(rowNum,colNum);
//   arma::field<List>::iterator start = result.begin();
//   arma::field<List>::iterator end = result.end();
//   for(arma::field<List>::iterator it = start; it != end; ++it)
//   { 
//     int idx = it.i;
//     Rcout<< "idx: "<< idx << std::endl;
//     (*it) = x(idx,0);
//   }
//   return result;
// }
// 
// arma::field<NumericVector> AdjustFieldVectorStructure(arma::field<NumericVector> x, int rowNum, int colNum){
//   arma::field<NumericVector> result(rowNum,colNum);
//   arma::field<NumericVector>::iterator start = result.begin();
//   arma::field<NumericVector>::iterator end = result.end();
//   for(arma::field<NumericVector>::iterator it = start; it != end; ++it)
//   { 
//     int idx = it.i ;
//     Rcout<< "idx: "<< idx << std::endl;
//     (*it) = x(idx,0);
//   }
//   return result;
// }
