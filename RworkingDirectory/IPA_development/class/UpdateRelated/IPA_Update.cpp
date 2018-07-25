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

// obsMostIntAdd, otherActuAddComp_uvec, otherActuAddInt
List UpdateRecMostIntAddInfo(List recMostIntAdd, arma::rowvec actuAddComp, arma::rowvec actuAddInt){
  double actuMaxInt = actuAddInt.max();
  // double oldMaxInt = obsMostIntAdd["intValue"];
  int oldAddCompNum = recMostIntAdd["compNum"];
  arma::uvec actuMaxIntIdx = find(actuAddInt == actuMaxInt);
  if (actuMaxIntIdx.n_elem == 1){
    int newAddCompNum = actuAddComp.at(0);
    if(oldAddCompNum == newAddCompNum){
      int recTime = recMostIntAdd["recTime"];
      recMostIntAdd["recTime"] = recTime + 1;
    } else{
      recMostIntAdd["recTime"] = 1;
      recMostIntAdd["compNum"] = newAddCompNum;
    }
  } else{                             //!!! here to code if observe the most intense adduct in multi adducts
    arma::uvec oldInNew_uvec = find(actuMaxIntIdx == oldAddCompNum);
    if(!oldInNew_uvec.is_empty()){    // !!! if the old most intense adduct is observed !!!
      int recTime = recMostIntAdd["recTime"];
      recMostIntAdd["recTime"] = recTime + 1;
    }else{                            // !!! if the old most intense adduct is not observed !!!   ???
      int newAddCompNum = actuMaxIntIdx.at(0);
      recMostIntAdd["recTime"] = 1;
      recMostIntAdd["compNum"] = newAddCompNum;
    }
  }
  return recMostIntAdd;
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
  Rcout<< "old iso: "<< old_iso << std::endl;
  Rcout<< "obs_iso_mean: "<< obs_iso_mean << std::endl;
  double new_iso = old_iso + (obs_iso_mean - old_iso) * scale;
  Rcout<< "new_iso: "<< new_iso << std::endl;
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

// void printOutRvec(arma::rowvec x){
//   // Rcout << "UVector" << std::endl;
//   for (int i = 0; i < x.n_elem; i++){
//     Rcout << "idx: " << i << "; value: " << x.at(i) <<std::endl;
//   }
// }

// void printOutUvec(arma::uvec x){
//   // Rcout << "UVector" << std::endl;
//   for (int i = 0; i < x.n_elem; i++){
//     Rcout << "idx: " << i << "; value: " << x.at(i) <<std::endl;
//   }
// }

// // [[Rcpp::export]]
// double testGetValueFromArmaUvec(NumericVector x){
//   arma::uvec result = as<arma::uvec>(x);
//   return result(0,2);
// }

// // [[Rcpp::export]]
// int testIndexMax(NumericVector x){
//   arma::rowvec x_rvec = as<arma::rowvec>(x);
//   arma::uword result = x_rvec.index_max();
//   return result;
// }

// // [[Rcpp::export]]
// arma::uvec findOldInNew(NumericVector x, int a){
//   arma::uvec x_uvec = as<arma::uvec>(x);
//   arma::uvec result = find(x_uvec == a);
//   return result;
// }

// [[Rcpp::export]]
List UpdateMain(arma::sp_mat post_spMat,                // posterior matrix
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
                List obsIR,                             // storage of observed intensity ratios of all detected isotope pairs in this update time
                NumericVector obsIRlComp,               // left compound number of stored observed intensity ratio
                NumericVector obsIRrComp,               // right compound number of stored observed intensity ratio
                List recIRdistr,                        // record all distributions of all isotope pairs
                NumericVector recIRdistrLComp,          // left compound number of all recorded distributions of all isotope pairs
                NumericVector recIRdistrRComp,          // right compound number of all recorded distributions of all isotope pairs
                NumericVector monoAddComp,              // compound number of all monoisotopic adducts
                List allAdd,                          // compound number of all other adducts except monoisotopic adducts
                List recMostIntAddInfo,                 // observed most intense adduct information of each monoisotopic adduct
                NumericVector recMainAdd,               // observed main adduct information of each monoisotopic adduct
                List obsRT,                             // storage of observed retention time in this update time
                List recRT,                             // record all observed retention time
                List recRTdistr,                        // record all the distributions of retention time
                NumericVector obsRTcomp,                // the compound numbers of the observed retention time in this update time
                NumericVector recRTdistrComp,           // the compound numbers of recorded retention time distributions
                double ma,                              // mass accuracy
                double s_pk,                            // update scale(speed) for pk
                double s_ma,                            // update scale for ma
                double s_iso,                           // update scale for iso
                int t_mainAdd,                          // threshold observed time value for main adduct shuffle
                double threshold                        // threshold used to filter posterior
                  
){
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>>>START THE WHOLE UPDATE STAGE" << std::endl;

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
      int select_mass = messageAndChoice(dupCompNum, dupMass_rvec, massMZ_rvec, massRT_rvec);    // get the choice, next delete data accordingly
      if (select_mass == -1){                                                                    // if user cannot decide which mass to choose
        arma::uvec retainComp_uvec = find(compSet_rvec != dupCompNum);                           // delete this duplicated compound
        compSet_rvec = compSet_rvec.elem(retainComp_uvec).t();                                   // filtered compound set
        massSet_rvec = massSet_rvec.elem(retainComp_uvec).t();                                   // filtered mass set
        messageAfterChoice();                                                                    // print some information after user making the choice each time
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
        messageAfterChoice();                                                                    // print some information after user making the choice each time
      }
    }
  }

  Rcout << ">>>>>>>>>>>>>>>>>>>>>>OVERLAPPED ASSIGNMENTS ISSUE FIXED, CONTINUE THE UPDATE STAGE" << std::endl;

  // based on assignment pair after thresholding posterior probability matrix: 
  // which are: mass number and compound number,
  // do: update prior knowledge,
  //     record all observed mass accuracies in this update time before updating
  //     update observed most intense adduct information and observed main adduct information of monoisotopic adduct
  //     and record all observed intensity ratios of observed isotopes in this update time before updating
  Rcout << ">>>>>>>START RECORD OBSERVED MASS ACCURACIES, RETENTION TIME AND INTENSITY RATIOS;" << std::endl;
  Rcout << "       PREPARE FOR MAIN ADDUCT SHUFFLE, AND UPDATE PRIOR KNOWLEDGE" << std::endl;
  arma::rowvec pk_rvec = as<arma::rowvec>(pk);                                                          // prior knowledge vector
  arma::rowvec massInt_rvec = as<arma::rowvec>(massInt);                                                // intensity value vector of masses
  arma::rowvec obsRTComp_rvec = as<arma::rowvec>(obsRTcomp);                                            // compound numbers of observed retention times
  arma::rowvec obsIRlComp_rvec = as<arma::rowvec>(obsIRlComp);                                          // left part of compound numbers of observed intensity values
  arma::rowvec obsIRrComp_rvec = as<arma::rowvec>(obsIRrComp);                                          // right part of compound numbers of observed intensity values
  arma::rowvec monoAddComp_rvec = as<arma::rowvec>(monoAddComp);                                        // compound number indices of all the monoisotopic adducts
  arma::uvec compSet_uvec = as<arma::uvec>(wrap(compSet_rvec));
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
      int obsRTCompFind = obsRTCompFind_uvec.at(0);
      NumericVector obsRTsub = obsRT[obsRTCompFind];
      NumericVector recRTsub = recRT[obsRTCompFind];
      obsRTsub.push_back(rtValue);
      recRTsub.push_back(rtValue);
      obsRT[obsRTCompFind] = obsRTsub;
      recRT[obsRTCompFind] = recRTsub;
    }
    
    // next: update observed most intense adduct information of each monoisotopic adduct
    // theCompLeft, theMassLeft
    arma::uvec monoAddFind_uvec = find(monoAddComp_rvec == theCompLeft);               // judge if the compound is a monoisotopic adduct
    if(!monoAddFind_uvec.is_empty()){                                                  // if the selected compound is a monoisotopic adduct
      // get all other types of adducts connected with this monoisotopic adduct, which are theoretical indices
      int monoAddFind = monoAddFind_uvec.at(0);
      NumericVector allTheoAddComp = allAdd.at(monoAddFind);
      
      // get all other types of adducts connected with this monoisotopic adduct, which are actual indices after threshold
      arma::uvec allTheoAddComp_uvec = as<arma::uvec>(allTheoAddComp);
      arma::uvec allActuAddComp_uvec = arma::intersect(allTheoAddComp_uvec, compSet_uvec);
      NumericVector allActuAddInt;
      NumericVector allActuAddComp;
      // if we can get any other types of adducts connected with this monoisotopic adduct after threshold
      // iterate and update observed most intense adduct information of each monoisotopic adduct
      if(!allActuAddComp_uvec.is_empty()){     
        List recMostIntAdd = recMostIntAddInfo.at(monoAddFind);
        for(int i = 0; i < allActuAddComp_uvec.n_elem; i++){
          int otherAddCompNum = allActuAddComp_uvec.at(i);
          arma::uvec otherAddCompFind = find(compSet_rvec == otherAddCompNum);
          int otherAddCompidx = otherAddCompFind.at(0);
          int otherAddMassNum = massSet_rvec.at(otherAddCompidx);
          allActuAddComp.push_back(otherAddCompNum);
          allActuAddInt.push_back(massInt_rvec.at(otherAddMassNum));
          
        }
        arma::rowvec allActuAddInt_rvec = as<arma::rowvec>(allActuAddInt);
        arma::rowvec allActuAddComp_rvec = as<arma::rowvec>(allActuAddComp);
        recMostIntAdd = UpdateRecMostIntAddInfo(recMostIntAdd, allActuAddComp_rvec, allActuAddInt_rvec);
        recMostIntAddInfo.at(monoAddFind) = recMostIntAdd;
      }
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
          NumericVector obsIRlSub;
          NumericVector obsIRrSub;
          obsIRlSub.push_back(theIR_LR);
          obsIRrSub.push_back(theIR_RL);
          obsIR.push_back(obsIRlSub);
          obsIR.push_back(obsIRrSub);
          obsIRlComp.push_back(theCompLeft);
          obsIRlComp.push_back(theCompRight);
          obsIRrComp.push_back(theCompRight);
          obsIRrComp.push_back(theCompLeft);
        } else{                                                                                 // compound number pair is already recorded
          int obsIRisoFind = obsIRisoFind_uvec.at(0);
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

  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO UPDATE MASS ACCURACY" << std::endl;
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

  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO UPDATE ISOTOPE CONNECTIVITIES" << std::endl;
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
        double oldIntValueLR = iso_spMat.at(obsIRlCompNum, obsIRrCompNum);
        double oldIntValueRL = iso_spMat.at(obsIRrCompNum, obsIRlCompNum);
        iso_spMat.at(obsIRlCompNum, obsIRrCompNum) = UpdateIso(oldIntValueLR, theIR_LR_mean, s_iso);
        iso_spMat.at(obsIRrCompNum, obsIRlCompNum) = UpdateIso(oldIntValueRL, theIR_RL_mean, s_iso);
        
        // record the data indices
        remIRidx.push_back(i);
        remIRidx.push_back(i + 1);
      } 
    } else{  
      // if we have recorded the distribution information of this observed compound number pair before,
      // we pop the distribution information into the same position we find it, 
      // but we don't pop the indices in this situation cause the indices are already existed
      int recIRdistrIsoFind = recIRdistrIsoFind_uvec.at(0);
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
  
  
  // next we shuffle the main adduct relations based on the observed adduct intensities of monoisotopic adducts
  // monoAddComp, otherAddComp, recMostIntAddInfo, obsMainAddInfo
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO SHUFFLE MAIN ADDUCT RELATIONS" << std::endl;
  // based on the recorded most intense adduct information of each mono-isotopic adduct, update main adduct information
  for(int i = 0; i < recMostIntAddInfo.length(); i++){
    List recMostIntAdd = recMostIntAddInfo.at(i);
    int recMostIntCompNum = recMostIntAdd["compNum"];
    int recMainAddCompNum = recMainAdd.at(i);
    if(recMainAddCompNum != recMostIntCompNum){
      int recTime = recMostIntAdd["recTime"];
      if(recTime >= t_mainAdd){
        recMainAdd.at(i) = recMostIntCompNum;
      }
    }
  }
  
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO RECORD RETENTION TIME DISTRIBUTION" << std::endl;
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
      int recRTdistrCompFind = recRTdistrCompFind_uvec.at(0);
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
  


  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>END OF WHOLE UPDATE STAGE" << std::endl;

  // next is the return value set, what should we return?
  // pk: prior knowledge vector
  // ma: updated mass accuracy value, observed mass accuracies at this time, all time mass accuracy distributions, all time mass accuracy values 
  // iso: observed intensity ratios at this time + compound pair indices, plus observed intensity ratio numbers at this time 
  //      recorded distributions of intensity ratios of all thime + compound pair indices
  // add: 
  // rt: observed RT in this time, compound number of observed RT in this time, recorded RT of all time, recorded distributions of RT all time, recorded compound number of all thime recorded RT distributions

  // normalise prior knowledge
  pk_rvec = pk_rvec / sum(pk_rvec);

  // return a list, contain all the updated results
  List model = List::create(Rcpp::Named("pk") = wrap(pk_rvec),
                            Rcpp::Named("ma") = ma,
                            Rcpp::Named("iso") = iso_spMat);
  
  List obs_MA = List::create(Rcpp::Named("obsMA") = obsMA);
  
  List rec_MA = List::create(Rcpp::Named("recMA") = recMA,
                             Rcpp::Named("recMAdistr") = recMAdistr);
  
  List obs_IR = List::create(Rcpp::Named("obsIR") = obsIR,
                             Rcpp::Named("obsIRlComp") = obsIRlComp,
                             Rcpp::Named("obsIRrComp") = obsIRrComp);
  
  List rec_IR = List::create(Rcpp::Named("recIRdistr") = recIRdistr,
                             Rcpp::Named("recIRdistrLComp") = recIRdistrLComp,
                             Rcpp::Named("recIRdistrRComp") = recIRdistrRComp);
  
  List obs_RT = List::create(Rcpp::Named("obsRT") = obsRT,
                             Rcpp::Named("obsRTcomp") = obsRTcomp);
  
  
  List rec_RT = List::create(Rcpp::Named("recRT") = recRT,
                             Rcpp::Named("recRTdistr") = recRTdistr,
                             Rcpp::Named("recRTdistrComp") = recRTdistrComp);
  
  List rec_ADD = List::create(Rcpp::Named("recMostIntAddInfo") = recMostIntAddInfo,
                              Rcpp::Named("recMainAdd") = recMainAdd);

  List record = List::create(Rcpp::Named("ma") = rec_MA,
                             Rcpp::Named("ir") = rec_IR,
                             Rcpp::Named("add") = rec_ADD,
                             Rcpp::Named("rt") = rec_RT);
  
  List other = List::create(Rcpp::Named("ma") = obs_MA,
                             Rcpp::Named("ir") = obs_IR,
                             Rcpp::Named("rt") = obs_RT);
  
  return List::create(Rcpp::Named("model") = model,
                      Rcpp::Named("record") = record,
                      Rcpp::Named("other") = other);

  // return R_NilValue;
}

