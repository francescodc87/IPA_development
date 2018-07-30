#include <RcppArmadillo.h>
using namespace Rcpp;
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

List UpdateIsoInfo(List IR_mean_std, double obs_mean, double obs_logMean, double obs_logStd, int obs_num){
  int old_num = IR_mean_std["iso_num"];
  double old_mean = IR_mean_std["iso_mean"];
  double old_logMean = IR_mean_std["iso_logMean"];
  double old_logStd = IR_mean_std["iso_logStd"];
  int new_num = old_num + obs_num;
  // use basic statistical knowledge, calculate new mean and std
  double new_mean = pow(pow(old_mean, old_num) * pow(obs_mean, obs_num), 1. / new_num );
  double new_logMean = (old_logMean * old_num + obs_logMean * obs_num) / new_num;
  double new_logVar = ((old_num - 1) * pow(old_logStd,2) + (obs_num - 1) * pow(obs_logStd, 2) + old_num * pow(old_logMean, 2) + obs_num * pow(obs_logMean, 2) - new_num * pow(new_logMean, 2)) / (new_num - 1);
  double new_logStd = sqrt(new_logVar);
  IR_mean_std["iso_num"] = new_num;
  IR_mean_std["iso_mean"] = new_mean;
  IR_mean_std["iso_logMean"] = new_logMean;
  IR_mean_std["iso_std"] = exp(new_logStd);
  IR_mean_std["iso_logStd"] = new_logStd;
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
NumericVector UpdateAddMostIntTime(NumericVector intIdx, arma::rowvec intValue, NumericVector mostIntTime){
  double maxInt = intValue.max();
  arma::uvec findMaxInt_uvec = arma::find(intValue == maxInt);
  for(int i = 0; i < findMaxInt_uvec.n_elem; i++){
    int maxInt_idx = intIdx.at(findMaxInt_uvec.at(i));
    mostIntTime.at(maxInt_idx) = mostIntTime.at(maxInt_idx) + 1;
  }
  return mostIntTime;
}

double UpdateAddPriorKnowledge(double old_pk,
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
  double x = 1 / scale; 
  double new_iso = pow((pow(old_iso, x) * obs_iso_mean), 1. / (x + 1));
  return new_iso;
}

int MsgAndChoWhenOverlapAssignDetected(int num){
  Rcout << " " << std::endl;
  Rcout << "-----------------------------------------------------!!!WARNING!!! " << num << " sets of overlapped assignments detected--------------------------------------------" << std::endl;
  Rcout << "This step is designed to initially decrease the number of overlapped assignments, it is not the step of choose masses in overlapped assignments" << std::endl;
  Rcout << "Type 0 to choose ignoring all overlapped assignments, then all overlapped assignments will be deleted" << std::endl;
  Rcout << "Type 1 to choose software help, the rule of filtering overlapped assignments are: " << std::endl;
  Rcout << "     I. prefer assignments in which the masses show higher intensities" << std::endl;
  Rcout << "     II. prefer assignments in which detecting lower massaccuracies" << std::endl;
  Rcout << "     III. prefer assignments in which more connections in formula levels are found based on the mass retention times" << std::endl;
  Rcout << "Please insert your choice in the next line and press 'Enter'" << std::endl;
  Environment base = Environment("package:base");
  Function readline = base["readline"];
  Function as_numeric = base["as.numeric"];
  int result = as<int>(as_numeric(readline("> ")));
  Rcout << "Your choice is successfully accepted!" << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;
  return result;
}

void MsgWhenIgnoreAllOverlapAssign(){
  Rcout << " " << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << "This is a hint that you have choosed to ignore all overlapped assignments" << std::endl;
  Rcout << "So the whole update stage will operate normally with the delete of all overlapped assignments, and there is no need to select masses in overlapped assignments" << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;
}

void MsgWhenChooseSoftwareHelp(){
  Rcout << " " << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << "This is a hint that you have choosed software help to decrease the number of overlapped assignments" << std::endl;
  Rcout << "Next step is to choose the masses in the remained overlapped assignments manually" << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;
}

int MsgAndChoiceOfMasses(int fomuNum, arma::rowvec massNum, arma::rowvec score){
  Rcout << " " << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << "Here are all the information of masses in an overlapped assignment, please look carefully and then choose the correct assignment by type the mass number" << std::endl;
  Rcout << "Once the input number is inserted, other assignments related to the same chemical fomula will be deleted" << std::endl;
  Rcout << "However, if you can not decide which one is the correct assignment, type 0, then all assignments related to this chemical fomula will be deleted" << std::endl;
  Rcout << "    Fomula number: " << fomuNum + 1 << std::endl;
  for(int i = 0; i < massNum.n_elem; i++){
    // Rcout << "    Mass number: " << massNum.at(i) + 1 << ";  m/z: " << mz.at(massNum.at(i)) << ";  RT: " << rt.at(massNum.at(i)) << std::endl;
    Rcout << "    Mass number: " << massNum.at(i) + 1 << std::endl;
  }
  Rcout << "Please insert your choice (mass number) in the next line and press 'Enter'" << std::endl;
  Environment base = Environment("package:base");
  Function readline = base["readline"];
  Function as_numeric = base["as.numeric"];
  int result = as<int>(as_numeric(readline("> ")));
  Rcout << "Your choice is successfully accepted!" << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;
  return result - 1;
}

void MsgAfterChoiceOfMasses(){
  Rcout << " " << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << "If you don't see any more print information, then there is no more overlapped assignment" << std::endl;
  Rcout << "Otherwise, please keep solving overlapped assignments" << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;
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

arma::rowvec UpdateCompPriorKnowledge(arma::rowvec pkComp_rvec, arma::uvec idx_uvec, double scale){
  pkComp_rvec.elem(idx_uvec).operator+=(scale);
  return pkComp_rvec;
}

int CompareToShuffleMainAdd(int oldMainAddFomu, NumericVector addFomu, NumericVector mostIntTime){
  arma::rowvec mostIntTime_rvec = as<arma::rowvec>(mostIntTime);
  int maxValue = mostIntTime_rvec.max();
  arma::uvec findMax_uvec = arma::find(mostIntTime_rvec == maxValue);
  int result;
  if (findMax_uvec.n_elem == 1){
    result = addFomu.at(findMax_uvec.at(0));
  }else{
    arma::rowvec addFomu_rvec = as<arma::rowvec>(addFomu);
    arma::rowvec mainAddCandidate_rvec = addFomu_rvec.elem(findMax_uvec).t();
    arma::uvec findOldMainAdd_uvec = arma::find(mainAddCandidate_rvec == oldMainAddFomu);
    if (findOldMainAdd_uvec.n_elem != 0){
      result = oldMainAddFomu;
    }else{
      // here to choose the main adduct when main adduct candidates do not contain the old main adduct
      // impossible to happen
      result = mainAddCandidate_rvec.at(0);
    }
  }
  
  return result;
}

// geometric weighted mean
double GeoMean(double a, double b){
  return pow((a * b), 0.5);
}

double GeoStd(NumericVector intRatios){
  // exp(std(ln(A)))
  NumericVector intLn = log(intRatios);
  arma::rowvec intLn_rvec = as<arma::rowvec>(intLn);
  double std = arma::stddev(intLn_rvec);
  return exp(std);
}

double GeoWeiMean(double mean_a, int num_a, double mean_b, int num_b){
  return pow(pow(mean_a, num_a) * pow(mean_b, num_b), 1. / (num_a + num_b));
}

int CountLinkInSpMatRow(arma::sp_mat x, int rowNum, arma::uvec dupFomu, arma::rowvec allFomu, arma::rowvec allMass, arma::rowvec massRT, double rtValue, double t){
  arma::sp_mat::const_row_iterator xStart = x.begin_row(rowNum);
  arma::sp_mat::const_row_iterator xEnd   = x.end_row(rowNum);
  int result = 0;
  for(arma::sp_mat::const_row_iterator it = xStart; it != xEnd; ++it)
  { 
    int colNum = it.col();
    double ir = (*it);
    if(ir != 0){
      arma::uvec fomuFind = arma::find(allFomu == colNum);
      if(!fomuFind.is_empty()){
        arma::uvec dupFind = arma::find(dupFomu == colNum);
        if(dupFind.is_empty()){
          int idx = fomuFind.at(0);
          int massNum = allMass.at(idx);
          int obsRTvalue = massRT.at(massNum);
          if(fabs(rtValue - obsRTvalue) <= t){
            result += 1;
          }
        }
      }
    }
  }
  return result;
}

arma::uvec SelectLinkWinner(NumericVector isoLink, NumericVector addLink){
  NumericVector allLink = isoLink + addLink;
  arma::rowvec allLink_rvec = as<arma::rowvec>(allLink);
  double maxNum = arma::max(allLink_rvec);
  return arma::find(allLink_rvec == maxNum);
}


// void printOutRvec(arma::rowvec x){
//   // Rcout << "UVector" << std::endl;
//   for (int i = 0; i < x.n_elem; i++){
//     Rcout << "idx: " << i << "; value: " << x.at(i) <<std::endl;
//   }
// }
// 
// void printOutUvec(arma::uvec x){
//   // Rcout << "UVector" << std::endl;
//   for (int i = 0; i < x.n_elem; i++){
//     Rcout << "idx: " << i << "; value: " << x.at(i) <<std::endl;
//   }
// }

// [[Rcpp::export]]
List UpdateMain(arma::sp_mat post_spMat,                // posterior matrix
                arma::sp_mat iso_spMat,                 // iso matrix
                arma::sp_mat add_spMat,                 // add matrix
                arma::sp_mat bio_spMat,                 // bio matrix
                NumericVector fomuCompId,               // compound id of fomulas
                NumericVector pkFomu,                   // prior knowledge of all fomulas
                NumericVector pkComp,                   // prior knowledge of all compounds
                double ma,                              // mass accuracy
                NumericVector massMZ,                   // mass-charge-ratio values of masses
                NumericVector fomuMZ,                   // mass-charge-ratio values of chemical fomulas
                NumericVector obsMA,                    // storage of observed mass accuracies in this update time
                NumericVector recMA,                    // record of all observed mass accuracies
                List recMAdistr,                        // record of all the distributions of mass accuracies
                NumericVector massInt,                  // intensity values of masses
                List obsIR,                             // storage of observed intensity ratios of all detected isotope pairs in this update time
                NumericVector obsIRlFomu,               // storage of all the left fomula numbers of all the stored observed intensity ratios in this update time
                NumericVector obsIRrFomu,               // storage of all the right fomula numbers of all the stored observed intensity ratios in this update time
                List recIRdistr,                        // record of all the distributions of all recorded intensity ratios
                NumericVector recIRdistrLFomu,          // record of all the left fomula numbers of all recorded distributions
                NumericVector recIRdistrRFomu,          // record of all the right fomula numbers of all recorded distributions
                NumericVector recMonoFomu,              // record of all the fomula numbers of all monoisotopic adducts
                NumericVector recMonoMainAddFomu,       // record of all the fomula numbers of all monoisotopic adduct's main adducts
                List recMonoFomuWithInSameComp,         // the fomula numbers of all monoisotopic adducts of each monoisotopic compound (share the same compound)
                List recMITwithInSameComp,              // the most intense time values of monoisotopic adducts of each monoisotopic compound (share the same compound)
                NumericVector recCompMainAddFomu,       // fomula numbers of the main adducts of compounds
                List recCompMonoFomu,                   // the monoisotopic adduct fomula numbers of compounds
                List recCompBioLink,                    // the compound numbers of compounds (linked with biochemical reactions)
                NumericVector massRT,                   // mass retention time values of all the masses
                List obsFomuRT,                         // storage of observed retention time values of fomulas in this update time
                NumericVector obsRTfomu,                // the fomula numbers of the observed retention time in this update time
                List recFomuRT,                         // record all observed retention time values of fomulas
                NumericVector recRTfomu,                // the fomula numbers of the recorded retention time values
                List recFomuRTdistr,                    // record all the distributions of retention time values of fomulas
                NumericVector recRTdistrFomu,           // the fomula numbers of the recorded retention time distributions
                double s_pkAdd,                         // update scale(speed) for pkFomu
                double s_pkComp,                        // update scale(speed) for pkComp
                double s_ma,                            // update scale for ma (0< <=1)
                double s_iso,                           // update scale for iso (0< <=1)
                double t_RT,                            // threshold used to count iso and add connection when solving overlapped assignments automatically
                double t_post                           // threshold used to filter posterior

){
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>>>START THE WHOLE UPDATE STAGE" << std::endl;
  Rcout << " " << std::endl;

  NumericVector massSet;                                                  // mass number set after threshold
  NumericVector fomuSet;                                                  // fomula number set after threshold
  // NumericVector fomuCompIdSet;                                                // fomula id set after threshold
  int monoAddNum = recMonoFomu.length();                                  // monoisotopic number
  int compNum = pkComp.length();                                          // compound number

  // posterior probability matrix iterator creation
  arma::sp_mat::const_iterator postStart = post_spMat.begin();
  arma::sp_mat::const_iterator postEnd   = post_spMat.end();

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>START TO FILTER POSTERIOR MATRIX" << std::endl;
  Rcout << " " << std::endl;
  // filter the posterior probability matrix based on threshold, store qualified mass number, fomula number at the same position respectively
  for(arma::sp_mat::const_iterator it = postStart; it != postEnd; ++it)
  {
    int rowNum = it.row();
    int colNum = it.col();
    double p = (*it);
    if(p >= t_post){
      // storage step
      massSet.push_back(rowNum);
      fomuSet.push_back(colNum);
      // fomuCompIdSet.push_back(fomuCompId.at(colNum));
    }
  }

  // check if overlap assignments happen
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>START TO CHECK OVERLAPPED ASSIGNMENT" << std::endl;
  Rcout << " " << std::endl;
  arma::rowvec fomuSet_rvec = as<arma::rowvec>(fomuSet);
  arma::rowvec fomuMZ_rvec = as<arma::rowvec>(fomuMZ);
  arma::rowvec massSet_rvec = as<arma::rowvec>(massSet);
  arma::rowvec massMZ_rvec = as<arma::rowvec>(massMZ);
  arma::rowvec massRT_rvec = as<arma::rowvec>(massRT);
  arma::rowvec massInt_rvec = as<arma::rowvec>(massInt);
  // arma::rowvec fomuCompId_rvec = as<arma::rowvec>(fomuCompId);
  arma::rowvec fomuSetUni_rvec = unique(fomuSet_rvec);
  if(fomuSetUni_rvec.n_elem != fomuSet_rvec.n_elem){
    // get duplicated compound number set
    NumericVector fomuSetUni = wrap(fomuSetUni_rvec);
    arma::uvec fomuSet_uvec = as<arma::uvec>(fomuSet);
    arma::uvec fomuSetUni_uvec = as<arma::uvec>(fomuSetUni);
    arma::uvec dupFomu_uvec = SetDiff(fomuSet_uvec, fomuSetUni_uvec);
    arma::uvec dupFomuUni_uvec = unique(dupFomu_uvec);                      // that is the unique duplicated fomula numbers

    // first inspect the number of overlapped assignments and then make choice of decreasing the number
    int beforeSelectMass = MsgAndChoWhenOverlapAssignDetected(dupFomuUni_uvec.n_elem);
    if(beforeSelectMass == 0){
      MsgWhenIgnoreAllOverlapAssign();
      for(int i =0; i < dupFomuUni_uvec.n_elem; i++){
        int dupFomuNum = dupFomuUni_uvec.at(i);                                     // first get one of the duplicate fomula number
        arma::uvec dupFomuIdx_uvec = find(fomuSet_rvec == dupFomuNum);              // search the idx of the duplicate fomula number in fomula row vector
        arma::uvec allFomuIdx_uvec = find_finite(fomuSet_rvec);                     // all idx in fomula row vector
        arma::uvec retainFomuIdx_uvec = SetDiff(allFomuIdx_uvec, dupFomuIdx_uvec);  // the idx of retained fomula (which is also the idx of the retained masses)
        massSet_rvec = massSet_rvec.elem(retainFomuIdx_uvec).t();                   // the retained mass
        fomuSet_rvec = fomuSet_rvec.elem(retainFomuIdx_uvec).t();                   // the retained formula
      }
    } else{
      MsgWhenChooseSoftwareHelp();

      // first we follow the rules to decrease the numbers of overlapped assignments, which are:
      // I. prefer assignments in which the masses show higher intensities;
      // II. prefer assignments in which detecting lower massaccuracies;
      // III. prefer assignments in which more connections in formula levels are found based on the mass retention times;
      
      // first get all indices of non-duplicate fomula, pop them into selectedIdx
      NumericVector selectedIdx;
      arma::uvec nonDupFomu_uvec = SetDiff(fomuSetUni_uvec, dupFomuUni_uvec);
      for(int i = 0; i < nonDupFomu_uvec.n_elem; i++){
        int nonDupFomuNum = nonDupFomu_uvec.at(i);
        arma::uvec nonDupIdx_uvec = arma::find(fomuSet_rvec == nonDupFomuNum);
        for(int j = 0; j < nonDupIdx_uvec.n_elem; j++){
          selectedIdx.push_back(nonDupIdx_uvec.at(j));
        }
      }
      
      // then select mass of every duplicate fomula, pop the index into selectedIdx
      for(int i =0; i < dupFomuUni_uvec.n_elem; i++){
        int dupFomuNum = dupFomuUni_uvec.at(i);                                   // first get one of the duplicate fomula number
        arma::uvec dupFomuIdx_uvec = find(fomuSet_rvec == dupFomuNum);            // search the idx of the duplicate fomula number in fomula row vector
        arma::rowvec dupSubMassSet_rvec = massSet_rvec.elem(dupFomuIdx_uvec).t(); // get duplicated masses
        
        // now we have the sub duplicated set of masses/fomula from each fomula's all assignments, the next is to perform score criterion and automatic select the mass with the highest score
        // prepare some data which will used in scoring
        NumericVector dupSubMassIntSet;
        NumericVector dupSubMAset;
        NumericVector dupSubMassRTset;
        NumericVector dupSubIsoLink;
        NumericVector dupSubAddLink;
        arma::rowvec dupSubScore_rvec = arma::zeros<arma::rowvec>(dupSubMassSet_rvec.n_elem);  // initial score to 0
        for(int j = 0; j < dupSubMassSet_rvec.n_elem; j++){
          int massNum = dupSubMassSet_rvec.at(j);
          double massIntValue = massInt_rvec.at(massNum);
          double fomuMZvalue = fomuMZ_rvec.at(dupFomuNum);
          double massMZvalue = massMZ_rvec.at(massNum);
          double maValue = ((massMZvalue - fomuMZvalue) * 1e6 ) / fomuMZvalue;
          double RTvalue = massRT_rvec.at(massNum);
          dupSubMassIntSet.push_back(massIntValue);
          dupSubMAset.push_back(fabs(maValue));
          dupSubMassRTset.push_back(RTvalue);
          int isoLink = CountLinkInSpMatRow(iso_spMat, dupFomuNum, dupFomuUni_uvec, fomuSet_rvec, massSet_rvec, massRT_rvec, RTvalue, t_RT);
          int addLink = CountLinkInSpMatRow(add_spMat, dupFomuNum, dupFomuUni_uvec, fomuSet_rvec, massSet_rvec, massRT_rvec, RTvalue, t_RT);
          // int addLink = 0;
          dupSubIsoLink.push_back(isoLink);
          dupSubAddLink.push_back(addLink);
        }

        arma::rowvec dupSubMassIntSet_rvec = as<arma::rowvec>(dupSubMassIntSet);
        arma::rowvec dupSubMAset_rvec = as<arma::rowvec>(dupSubMAset);
        // arma::rowvec dupSubMassRTset_rvec = as<arma::rowvec>(dupSubMassRTset);

        // rule I: prefer assignments in which the masses show higher intensities
        double maxInt = arma::max(dupSubMassIntSet_rvec);
        arma::uvec dupSubMaxIntIdx = arma::find(dupSubMassIntSet_rvec == maxInt);
        dupSubScore_rvec.elem(dupSubMaxIntIdx) += 1;       // score add

        // rule II: prefer assignments in which detecting lower massaccuracies
        double minMA = arma::min(dupSubMAset_rvec);
        arma::uvec dupSubMinMAidx = arma::find(dupSubMAset_rvec == minMA);
        dupSubScore_rvec.elem(dupSubMinMAidx) += 2;       // score add

        // rule III: prefer assignments in which more connections in formula levels are found based on the mass retention times
        arma::uvec dupSubLinkWinner = SelectLinkWinner(dupSubIsoLink,dupSubAddLink);
        dupSubScore_rvec.elem(dupSubLinkWinner) += 1;    // score add
        
        // after thresholded by three rules, perform automatic selection of mass in the subvec of duplicate assignments by the final score
        // dupSubScore_rvec, dupSubFomuSet_rvec, dupSubMassSet_rvec, dupFomuIdx_uvec
        int maxScore = arma::max(dupSubScore_rvec);
        arma::uvec maxScoreFind_uvec = arma::find(dupSubScore_rvec == maxScore);
        if(maxScoreFind_uvec.n_elem == 1){
          selectedIdx.push_back(dupFomuIdx_uvec.at(maxScoreFind_uvec.at(0)));
        } else{
          arma::rowvec manualDupMass_rvec = dupSubMassSet_rvec.elem(maxScoreFind_uvec).t();
          arma::rowvec manualDupScore_rvec = dupSubScore_rvec.elem(maxScoreFind_uvec).t();
          arma::uvec manualDupFomuIdx_uvec = dupFomuIdx_uvec.elem(maxScoreFind_uvec);
          int select_mass = MsgAndChoiceOfMasses(dupFomuNum, manualDupMass_rvec, manualDupScore_rvec);    
          if(select_mass != -1){                                                                              
            arma::uvec selectMassFind_uvec = find(manualDupMass_rvec == select_mass);
            selectedIdx.push_back(dupFomuIdx_uvec.at(selectMassFind_uvec.at(0)));
            MsgAfterChoiceOfMasses();   // print some information after user making the choice each time                                                                 
          } 
        }
      }
      arma::uvec selectedIdx_uvec = as<arma::uvec>(selectedIdx);
      fomuSet_rvec = fomuSet_rvec.elem(selectedIdx_uvec).t();                       // filtered fomula set
      massSet_rvec = massSet_rvec.elem(selectedIdx_uvec).t();                       // filtered mass set
    }
  }

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>OVERLAPPED ASSIGNMENTS ISSUE FIXED, CONTINUE THE UPDATE STAGE" << std::endl;
  Rcout << " " << std::endl;
  // based on assignment pair after thresholding posterior probability matrix:
  // which are: mass numbers and fomula numbers,
  // do: update prior knowledge (both for fomula and compound),
  //     record all observed mass accuracies in this update time before updating mass accuracy
  //     and record all observed intensity ratios of observed isotope pairs in this update time before updating iso matrix
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START RECORD OBSERVED MASS ACCURACIES, RETENTION TIME AND INTENSITY RATIOS;" << std::endl;
  Rcout << "                       UPDATE MOST INTENSE TIME TO PREPARE MAIN ADDUCT SHUFFLE, AND UPDATE PRIOR KNOWLEDGE" << std::endl;
  Rcout << " " << std::endl;
  arma::rowvec pkFomu_rvec = as<arma::rowvec>(pkFomu);                                                  // prior knowledge of all the fomulas
  arma::rowvec pkComp_rvec = as<arma::rowvec>(pkComp);                                                  // prior knowledge of all the compounds
  arma::rowvec obsRTfomu_rvec = as<arma::rowvec>(obsRTfomu);                                            // fomula numbers of observed retention times
  arma::rowvec obsIRlFomu_rvec = as<arma::rowvec>(obsIRlFomu);                                          // left part of fomula numbers of observed intensity values
  arma::rowvec obsIRrFomu_rvec = as<arma::rowvec>(obsIRrFomu);                                          // right part of fomula numbers of observed intensity values
  arma::uvec fomuSet_uvec = as<arma::uvec>(wrap(fomuSet_rvec));

  // update compound prior knowledge
  arma::uvec compMainAdd_uvec = as<arma::uvec>(recCompMainAddFomu);                                           // fomula numbers of all the compound main adducts
  arma::uvec compMainAddFind_uvec = arma::intersect(compMainAdd_uvec, fomuSet_uvec);                       // we want to see how many compound main adducts survive after thresholding the posterior matrix
  if(!compMainAddFind_uvec.is_empty()){
    // if any compound main adducts survive, update compound prior knowledge
    for(int i = 0; i < compMainAddFind_uvec.n_elem; i++){
      int compMainAddFomuNum = compMainAddFind_uvec.at(i);
      pkComp_rvec = UpdateCompPriorKnowledge(pkComp_rvec, arma::find(compMainAdd_uvec == compMainAddFomuNum), s_pkComp);
    }
  }

  for(int i = 0; i < fomuSet_rvec.n_elem; i++){                                                            // iterate on thresholded fomula numbers
    // after detect and fix overlap assignments, Update fomula Prior Knowledge
    int theLeftFomuNum = fomuSet_rvec.at(i);                                                               // first time loop, take the value as left part
    int theLeftMassNum = massSet_rvec.at(i);                                                               // first time loop, take the value as left part
    pkFomu_rvec.at(theLeftFomuNum) = UpdateAddPriorKnowledge(pkFomu_rvec.at(theLeftFomuNum), s_pkAdd);     // update fomula prior knowledge

    // calculate and record observed mass accuracies (negative values allowed)
    double maValue = ((massMZ[theLeftMassNum] - fomuMZ[theLeftFomuNum]) * 1e6 ) / fomuMZ[theLeftFomuNum];
    obsMA.push_back(maValue);
    recMA.push_back(maValue);

    // record observed retention times
    double rtValue = massRT(theLeftMassNum);

    // check if we observed other retention time of this fomula number before
    // if we have not observed before, push_back retention time value, and push_back the fomula number
    // else we pop the retention time value into the same position we find it, no need to pop_up the compound time cause we have already seen it in obsRTfomu_rvec
    arma::uvec obsRTFomuFind_uvec = find(obsRTfomu_rvec == theLeftFomuNum);
    if(obsRTFomuFind_uvec.is_empty()){
      NumericVector rt;
      rt.push_back(rtValue);
      obsFomuRT.push_back(rt);
      recFomuRT.push_back(rt);
      recRTfomu.push_back(theLeftFomuNum);
      obsRTfomu.push_back(theLeftFomuNum);
    } else{
      int obsRTFomuNum = obsRTFomuFind_uvec.at(0);
      NumericVector obsRTsub = obsFomuRT[obsRTFomuNum];
      NumericVector recRTsub = recFomuRT[obsRTFomuNum];
      obsRTsub.push_back(rtValue);
      recRTsub.push_back(rtValue);
      obsFomuRT[obsRTFomuNum] = obsRTsub;
      recFomuRT[obsRTFomuNum] = recRTsub;
    }

    // record the observed intensity ratios before updating iso matrix, we use another iteration of fomuSet_rvec to create isotope pairs
    for(int j = 0; j < fomuSet_rvec.n_elem; j++){
      // generate the right part of fomula numbers which will be used as indices to record intensity ratios
      // and generate the right part of masses which will be used as indices to calculate observed intensity ratios
      int theFomuRight = fomuSet_rvec.at(j);
      int theMassRight = massSet_rvec.at(j);

      // next is to combine the fomula number pair as isotope pair, use mass number pair to compute isotope ratio
      // and complete the storage task before update
      // about the next line judgement code: here we decrease the computational complexity in the second iteration, actually I iterate half, as you can see from the condition in the next if() statement
      //                                     & the condition also considerds that only the fomula which is an isotope pair can be applied in the following task
      if((theLeftFomuNum < theFomuRight) && (iso_spMat.at(theLeftFomuNum,theFomuRight) != 0)){
        double theIRvalue_LR = massInt_rvec.at(theLeftMassNum) / massInt_rvec.at(theMassRight);              // isotope ratio (Left part to right part)
        double theIRvalue_RL = 1 / theIRvalue_LR;                                                                 // isotope ratio (right part to left part)
        // we want to know if we have recorded any isotope ratios of this fomula pair before
        // so we take the intersect
        arma::uvec obsIRlFomuFind_uvec = find(obsIRlFomu_rvec == theLeftFomuNum);
        arma::uvec obsIRrFomuFind_uvec = find(obsIRrFomu_rvec == theFomuRight);
        arma::uvec obsIRisoFind_uvec = arma::intersect(obsIRlFomuFind_uvec, obsIRrFomuFind_uvec);
        // if we have not recorded the observed isotope ratio of this fomula number pair before, we push_back those value and indices
        // else we pop the value into the same position we find it, but we don't pop the indices in this situation cause the indices are already existed
        if(obsIRisoFind_uvec.is_empty()){
          NumericVector obsIRlSub;
          NumericVector obsIRrSub;
          obsIRlSub.push_back(theIRvalue_LR);
          obsIRrSub.push_back(theIRvalue_RL);
          obsIR.push_back(obsIRlSub);
          obsIR.push_back(obsIRrSub);
          obsIRlFomu.push_back(theLeftFomuNum);
          obsIRlFomu.push_back(theFomuRight);
          obsIRrFomu.push_back(theFomuRight);
          obsIRrFomu.push_back(theLeftFomuNum);
        } else{    // already spoken situation that fomula number pair is already recorded
          int obsIRisoFind = obsIRisoFind_uvec.at(0);
          NumericVector obsIRlSub = obsIR[obsIRisoFind];
          NumericVector obsIRrSub = obsIR[obsIRisoFind + 1];
          obsIRlSub.push_back(theIRvalue_LR);
          obsIRrSub.push_back(theIRvalue_RL);
          obsIR[obsIRisoFind] = obsIRlSub;
          obsIR[obsIRisoFind + 1] = obsIRrSub;
        }
      }
    }
  }

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>UPDATE OBSERVED MOST INTENSE TIME VALUES OF MONOISOTOPIC ADDUCTS" << std::endl;
  Rcout << " " << std::endl;
  // update recorded most intense time values of each monoisotopic adduct
  for(int i = 0; i < recMonoFomu.length(); i++){
    // get all other adducts connected with this monoisotopic adduct (linked by the same compound)
    NumericVector monoAddWithInSameComp = recMonoFomuWithInSameComp.at(i);                            // get all adducts except isotopes of this monoisotopic adduct
    arma::uvec monoAddWithInSameComp_uvec = as<arma::uvec>(monoAddWithInSameComp);

    // next: decide which one is the most intense, and update the most intense time in recMITwithInSameComp
    if(!monoAddWithInSameComp_uvec.is_empty()){
      NumericVector monoAddInt;
      NumericVector monoAddIdx;
      arma::uvec monoAddAllIdx = find_finite(monoAddWithInSameComp_uvec.t());
      for(int j = 0; j < monoAddWithInSameComp_uvec.n_elem; j++){                            // iterate on all "other" monoisotopic adducts
        int MonoAddFomuNum = monoAddWithInSameComp_uvec.at(j);
        arma::uvec MonoAddIdxInFomu_uvec = find(fomuSet_rvec == MonoAddFomuNum);
        if(!MonoAddIdxInFomu_uvec.is_empty()){
          int MonoAddIdxInFomu = MonoAddIdxInFomu_uvec.at(0);
          int MonoAddMassNum = massSet_rvec.at(MonoAddIdxInFomu);
          monoAddInt.push_back(massInt_rvec.at(MonoAddMassNum));
          monoAddIdx.push_back(monoAddAllIdx.at(j));
        }
      }
      // based on actuMonoAdd_addIdx and actuMonoAdd_addInt, update recMITwithInSameComp
      if(monoAddIdx.length() != 0){
        arma::rowvec monoAdd_addInt_rvec = as<arma::rowvec>(monoAddInt);
        recMITwithInSameComp.at(i) = UpdateAddMostIntTime(monoAddIdx, monoAdd_addInt_rvec, recMITwithInSameComp.at(i));
      }
    }
  }

  // we have recorded observed mass accuracies and intensity ratios, it's time to update the mean and std values respectively, and update the mass accuracy based on the distributions respectively
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO UPDATE MASS ACCURACY" << std::endl;
  Rcout << " " << std::endl;
  // update mass accuracy, we record all mass accuracies, but when it talks to build distributions and update mass accuracy we only use the observed mass accuracies in this time
  arma::rowvec obsMA_rvec = as<arma::rowvec>(obsMA);
  if(recMAdistr.length() == 0){                                      // if we have not recorded any distribution of mass accuracy before
    if(obsMA_rvec.n_elem > 1){                                       // if observed more than one mass accuracies in this update time, we can pop at least the first distribution information of mass accuracy
      double ma_mean = arma::mean(arma::abs(obsMA_rvec));
      double ma_std = arma::stddev(arma::abs(obsMA_rvec));
      double ma_num = obsMA_rvec.n_elem;
      recMAdistr["ma_num"] = ma_num;
      recMAdistr["ma_mean"] = ma_mean;
      recMAdistr["ma_std"] = ma_std;
      ma = UpdateMassAccuracy(ma, ma_mean, s_ma);                    // update mass accuracy main calculation
      NumericVector tmp;                                             // clear observed mass accuracies if the observed mass accuracies have been used to build the recorded distribution
      obsMA = tmp;
    }
  } else{                                                                       // if we have recorded any distribution of mass accuracy before
    if(obsMA_rvec.n_elem > 1){                                                  // if have observed one mass accuracy value in this update time, then we can update the distribution of mass accuracy
      double ma_mean = arma::mean(arma::abs(obsMA_rvec));
      double ma_std = arma::stddev(arma::abs(obsMA_rvec));
      double ma_num = obsMA_rvec.n_elem;
      recMAdistr = UpdateMassAccuracyInfo(recMAdistr, ma_mean, ma_std, ma_num); // update mass accuracy distribution information
      ma = UpdateMassAccuracy(ma, ma_mean, s_ma);                               // update mass accuracy main calculation
      NumericVector tmp;                                                        // clear observed mass accuracies if the observed mass accuracies have been used to build the recorded distribution
      obsMA = tmp;
    }
  }

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO UPDATE ISOTOPE CONNECTIVITIES" << std::endl;
  Rcout << " " << std::endl;
  // update iso observed intensity ratios (if not enough number to build distributions), distribution information and iso connectivity matrix
  arma::rowvec obsNewIRlFomu_rvec = as<arma::rowvec>(obsIRlFomu);              // left part of the observed isotope pairs
  arma::rowvec obsNewIRrFomu_rvec = as<arma::rowvec>(obsIRrFomu);              // right part of the observed isotope pairs
  arma::rowvec recIRdistrLFomu_rvec = as<arma::rowvec>(recIRdistrLFomu);       // left part of the recorded isotope pairs (of the recorded distribution information)
  arma::rowvec recIRdistrRFomu_rvec = as<arma::rowvec>(recIRdistrRFomu);       // right part of the recorded isotope pairs (of the recorded distribution information)
  IntegerVector remIRidx;
  // iterate on all observed isotope pairs,
  // the idea is to find out if the observed isotope pair is already found in the isotope pairs of recorded distribution information about intensity ratios
  for(int i = 0; i < obsNewIRlFomu_rvec.n_elem; i = i + 2){                    // as you can see here we decrease the computational complexity
    int obsIRlFomuNum = obsNewIRlFomu_rvec.at(i);
    int obsIRrFomuNum = obsNewIRrFomu_rvec.at(i);
    NumericVector theIRvalue_LR = obsIR[i];                                    // the intensity ratio_left part right part
    NumericVector theIRvalue_RL = obsIR[i + 1];                                // the intensity ratio_right part left part
    arma::uvec recIRdistrLFomuFind_uvec = find(recIRdistrLFomu_rvec == obsIRlFomuNum);
    arma::uvec recIRdistrRFomuFind_uvec = find(recIRdistrRFomu_rvec == obsIRrFomuNum);
    arma::uvec recIRdistrIsoFind_uvec = arma::intersect(recIRdistrLFomuFind_uvec, recIRdistrRFomuFind_uvec);
    // if we have not recorded the distribution information of this observed fomula number pair before,
    // we can push_back the first distribution information of this observed fomula number pair and the indices
    // only if record more than 2 intensity ratios
    if(recIRdistrIsoFind_uvec.n_elem == 0){
      if(theIRvalue_LR.length() > 1){
        List theIRvalueDistr_LR;
        List theIRvalueDistr_RL;
        NumericVector theIRlogValue_LR = log(theIRvalue_LR);
        NumericVector theIRlogValue_RL = log(theIRvalue_RL);
        arma::rowvec theIRlogValue_LR_rvec = as<arma::rowvec>(theIRlogValue_LR);
        arma::rowvec theIRlogValue_RL_rvec = as<arma::rowvec>(theIRlogValue_RL);
        // calculate IR numbers, mean and std value for this isotope pair for the first time
        double theIRvalue_LR_mean = GeoMean(theIRvalue_LR.at(0), theIRvalue_LR.at(1));
        double theIRvalue_LR_logMean = mean(theIRlogValue_LR);
        double theIRvalue_LR_std = GeoStd(theIRvalue_LR);
        double theIRvalue_LR_logStd = arma::stddev(theIRlogValue_LR_rvec);
        double theIRvalue_RL_mean = GeoMean(theIRvalue_RL.at(0), theIRvalue_RL.at(1));
        double theIRvalue_RL_logMean = mean(theIRlogValue_RL);
        double theIRvalue_RL_std = GeoStd(theIRvalue_RL);
        double theIRvalue_RL_logStd = arma::stddev(theIRlogValue_RL_rvec);
        double theIRvalue_num = 2;
        theIRvalueDistr_LR["iso_num"] = theIRvalue_num;
        theIRvalueDistr_LR["iso_mean"] = theIRvalue_LR_mean;
        theIRvalueDistr_LR["iso_logMean"] = theIRvalue_LR_logMean;
        theIRvalueDistr_LR["iso_std"] = theIRvalue_LR_std;
        theIRvalueDistr_LR["iso_logStd"] = theIRvalue_LR_logStd;
        theIRvalueDistr_RL["iso_num"] = theIRvalue_num;
        theIRvalueDistr_RL["iso_mean"] = theIRvalue_RL_mean;
        theIRvalueDistr_RL["iso_logMean"] = theIRvalue_RL_logMean;
        theIRvalueDistr_RL["iso_std"] = theIRvalue_RL_std;
        theIRvalueDistr_RL["iso_logStd"] = theIRvalue_RL_logStd;
        recIRdistr.push_back(theIRvalueDistr_LR);
        recIRdistr.push_back(theIRvalueDistr_RL);
        recIRdistrLFomu.push_back(obsIRlFomuNum);
        recIRdistrLFomu.push_back(obsIRrFomuNum);
        recIRdistrRFomu.push_back(obsIRrFomuNum);
        recIRdistrRFomu.push_back(obsIRlFomuNum);

        // update iso matrix
        double oldIntValueLR = iso_spMat.at(obsIRlFomuNum, obsIRrFomuNum);
        double oldIntValueRL = iso_spMat.at(obsIRrFomuNum, obsIRlFomuNum);
        iso_spMat.at(obsIRlFomuNum, obsIRrFomuNum) = UpdateIso(oldIntValueLR, theIRvalue_LR_mean, s_iso);
        iso_spMat.at(obsIRrFomuNum, obsIRlFomuNum) = UpdateIso(oldIntValueRL, theIRvalue_RL_mean, s_iso);

        // record the data indices
        remIRidx.push_back(i);
        remIRidx.push_back(i + 1);
      }
    } else{
      // if we have recorded the distribution information of this observed fomula number pair before,
      // we pop the distribution information into the same position we find it,
      // but we don't pop the indices in this situation cause the indices are already existed
      int recIRdistrIsoIdx = recIRdistrIsoFind_uvec.at(0);
      double theIRvalue_LR_mean = theIRvalue_LR.at(0);
      double theIRvalue_RL_mean = theIRvalue_RL.at(0);

      // update distribution information
      recIRdistr[recIRdistrIsoIdx] = UpdateIsoInfo(recIRdistr[recIRdistrIsoIdx], theIRvalue_LR_mean, log(theIRvalue_LR_mean), 0, 1);
      recIRdistr[recIRdistrIsoIdx + 1] = UpdateIsoInfo(recIRdistr[recIRdistrIsoIdx + 1], theIRvalue_RL_mean, log(theIRvalue_RL_mean), 0, 1);

      // get new intensity ratio mean value
      List theIRvalueDistr_LR = recIRdistr[recIRdistrIsoIdx];
      List theIRvalueDistr_RL = recIRdistr[recIRdistrIsoIdx + 1];
      theIRvalue_LR_mean = theIRvalueDistr_LR["iso_mean"];
      theIRvalue_RL_mean = theIRvalueDistr_RL["iso_mean"];

      // update iso matrix
      iso_spMat.at(obsIRlFomuNum, obsIRrFomuNum) = UpdateIso(iso_spMat.at(obsIRlFomuNum, obsIRrFomuNum), theIRvalue_LR_mean, s_iso);
      iso_spMat.at(obsIRrFomuNum, obsIRlFomuNum) = UpdateIso(iso_spMat.at(obsIRrFomuNum, obsIRlFomuNum), theIRvalue_RL_mean, s_iso);
      remIRidx.push_back(i);
      remIRidx.push_back(i + 1);
    }
  }

  // after update iso, we delete some observed data which is used in building the distribution from obsIR, obsIRlFomu, obsIRrFomu
  IntegerVector obsIRidx = seq_len(obsIRlFomu.length()) - 1;
  arma::uvec obsIRidx_uvec = as<arma::uvec>(obsIRidx);
  arma::uvec remIRidx_uvec = as<arma::uvec>(remIRidx);
  arma::uvec retainIRidx_uvec = SetDiff(obsIRidx_uvec, remIRidx_uvec);
  IntegerVector retainIRidx = wrap(retainIRidx_uvec.t());
  obsIR = obsIR[retainIRidx];
  obsIRlFomu = obsIRlFomu[retainIRidx];
  obsIRrFomu = obsIRrFomu[retainIRidx];

  // next we shuffle the main adduct relations based on the observed adduct intensities of monoisotopic adducts
  // which include main adduct relations of monoisotopic adducts and main adduct relations of compounds
  // recCompMainAddFomu, recCompMonoFomu
  // recMonoMainAddFomu, recMonoFomuWithInSameComp, recMITwithInSameComp
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO SHUFFLE MAIN ADDUCT RELATIONS" << std::endl;
  Rcout << " " << std::endl;
  NumericVector monoNewMainAddFomu(monoAddNum);
  NumericVector compNewMainAddFomu(compNum);
  int monoFlag = 0;
  for(int i = 0; i < recCompMainAddFomu.length(); i++){
    NumericVector monoFomuVec = recCompMonoFomu.at(i);
    // based on the recorded most intense adduct time of each mono-isotopic adduct, update main adduct information
    int compNewMainAdd = CompareToShuffleMainAdd(recMonoMainAddFomu.at(monoFlag), recMonoFomuWithInSameComp.at(monoFlag), recMITwithInSameComp.at(monoFlag));
    for(int j = monoFlag; j < monoFlag + monoFomuVec.length(); j++){
      monoNewMainAddFomu.at(j) = compNewMainAdd;
    }
    compNewMainAddFomu.at(i) = compNewMainAdd;
    monoFlag = monoFlag + monoFomuVec.length();
  }

  // after shuffle the main adducts, we can update adduct connectivity and bio connectivity now
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO UPDATE ADDUCT CONNECTIVITIES" << std::endl;
  Rcout << " " << std::endl;
  // recMonoFomu, recMonoMainAddFomu, monoNewMainAddFomu, add_spMat
  for (int i = 0; i < recMonoFomu.length(); i++){
    int oldMainAddFomu = recMonoMainAddFomu.at(i);
    int newMainAddFomu = monoNewMainAddFomu.at(i);
    int monoAdd = recMonoFomu.at(i);
    if(oldMainAddFomu != newMainAddFomu && monoAdd != newMainAddFomu) {
      add_spMat.at(monoAdd, oldMainAddFomu) = 0;
      add_spMat.at(oldMainAddFomu, monoAdd) = 0;
      add_spMat.at(monoAdd, newMainAddFomu) = 1;
      add_spMat.at(newMainAddFomu, monoAdd) = 1;
    }
  }

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO UPDATE BIOCHEMICAL CONNECTIVITIES" << std::endl;
  Rcout << " " << std::endl;
  // recCompBioLink, recCompMainAddFomu, compNewMainAddFomu, bio_spMat
  for (int i = 0; i < recCompBioLink.length(); i++){
    NumericVector bioConnectVec = recCompBioLink.at(i);
    if (bioConnectVec.length() != 0){
      for (int j = 0; j < bioConnectVec.length(); j++){
        int compNum = bioConnectVec.at(j);
        int oldCompLMainAdd = recCompMainAddFomu.at(i);
        int oldCompRMainAdd = recCompMainAddFomu.at(compNum);
        int newCompLMainAdd = compNewMainAddFomu.at(i);
        int newCompRMainAdd = compNewMainAddFomu.at(compNum);
        if (!(oldCompLMainAdd == oldCompRMainAdd && newCompLMainAdd == newCompRMainAdd)){
          bio_spMat.at(oldCompLMainAdd, oldCompRMainAdd) = 0;
          bio_spMat.at(newCompLMainAdd, newCompRMainAdd) = 1;
        }
      }
    }
  }

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO RECORD RETENTION TIME DISTRIBUTION" << std::endl;
  Rcout << " " << std::endl;
  arma::rowvec obsNewRTfomu_rvec = as<arma::rowvec>(obsRTfomu);                  // all fomula numbers of observed retention time before update
  arma::rowvec recRTdistrFomu_rvec = as<arma::rowvec>(recRTdistrFomu);           // all fomula numbers of recorded distribution information of retention time before update

  // iterate on all observed fomula numbers of retention time,
  // the idea is to find out if the observed fomula is already found in the fomula set of recorded distribution information about retention time
  IntegerVector remRTidx;
  for(int i = 0; i < obsNewRTfomu_rvec.n_elem; i++){
    int obsRTfomuNum = obsNewRTfomu_rvec.at(i);
    arma::uvec recRTdistrFomuFind_uvec = find(recRTdistrFomu_rvec == obsRTfomuNum);
    NumericVector obsFomuRTsub = obsFomuRT[i];
    // if we have not recorded the distribution information of this observed fomula before,
    // we can push_back the first distribution information of this observed fomula and the index (fomula number)
    // only if record more than 2 retention time values of this fomula
    if(recRTdistrFomuFind_uvec.n_elem == 0){
      if(obsFomuRTsub.length() > 1){
        arma::rowvec obsFomuRT_rvec = as<arma::rowvec>(obsFomuRTsub);
        double rt_mean = arma::mean(obsFomuRT_rvec);
        double rt_std = arma::stddev(obsFomuRT_rvec);
        double rt_num = 2;
        List rt_distr;
        rt_distr["rt_num"] = rt_num;
        rt_distr["rt_mean"] = rt_mean;
        rt_distr["rt_std"] = rt_std;
        recFomuRTdistr.push_back(rt_distr);          // here we push_back the distribution information
        recRTdistrFomu.push_back(obsRTfomuNum);      // here we push_back the index (compound number)
        remRTidx.push_back(i);
      }
    }else{
      // if we have recorded the distribution information of this observed fomula before,
      // we pop the distribution information into the same position we find it,
      // but we don't pop the indices in this situation cause the indices are already existed
      int recRTdistrIdx = recRTdistrFomuFind_uvec.at(0);
      double rt_mean = obsFomuRTsub.at(0);
      double rt_std = 0;
      double rt_num = 1;

      // update retention time distribution information
      recFomuRTdistr[recRTdistrIdx] = UpdateRtInfo(recFomuRTdistr[recRTdistrIdx], rt_mean, rt_std, rt_num);
      remRTidx.push_back(i);
    }
  }

  // after record retention time distribution information, we remove those observed data from obsFomuRT, obsRTfomu
  IntegerVector obsRTidx = seq_len(obsRTfomu.length()) - 1;
  arma::uvec obsRTidx_uvec = as<arma::uvec>(obsRTidx);
  arma::uvec remRTidx_uvec = as<arma::uvec>(remRTidx);
  arma::uvec retainRTidx_uvec = SetDiff(obsRTidx_uvec, remRTidx_uvec);
  IntegerVector retainRTidx = wrap(retainRTidx_uvec.t());
  obsFomuRT = obsFomuRT[retainRTidx];
  obsRTfomu = obsRTfomu[retainRTidx];

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>END OF WHOLE UPDATE STAGE" << std::endl;
  Rcout << " " << std::endl;

  // next is the return value set, what should we return?
  // pkFomu: prior knowledge vector (both fomula and compound)
  // ma: updated mass accuracy value, observed mass accuracies at this time, all time mass accuracy distributions, all time mass accuracy values
  // iso: observed intensity ratios at this time + fomula pair indices, plus observed intensity ratio numbers at this time
  //      recorded distributions of intensity ratios of all thime + fomula pair indices
  //      updated iso matrix
  // add: fomula number of all monoisotopic adducts;
  //      fomula number of monoisotopic adduct's main adduct
  //      the adduct numbers from monoisotopic & adduct connection
  //      the most intense time of adduct being observed as the most adduct in monoisotopic & adduct connection
  //      updated add matrix
  // bio: list of compounds which has bio connection with each compound
  //      updated main adduct fomula number of each compound
  //      updated bio matrix
  // rt: observed RT in this time, fomula number of observed RT in this time, recorded RT of all time, recorded distributions of RT all time, recorded fomula number of all thime recorded RT distributions

  // normalise prior knowledge
  pkFomu_rvec = pkFomu_rvec / sum(pkFomu_rvec);
  pkComp_rvec = pkComp_rvec / sum(pkComp_rvec);

  // return a list, contain all the updated results
  List obsMAlist = List::create(Rcpp::Named("obsMA") = obsMA);

  List recMAlist = List::create(Rcpp::Named("recMA") = recMA,
                                Rcpp::Named("recMAdistr") = recMAdistr);

  List obsIRlist = List::create(Rcpp::Named("obsIR") = obsIR,
                                Rcpp::Named("obsIRlFomu") = obsIRlFomu,
                                Rcpp::Named("obsIRrFomu") = obsIRrFomu);

  List recIRlist = List::create(Rcpp::Named("recIRdistr") = recIRdistr,
                                Rcpp::Named("recIRdistrLFomu") = recIRdistrLFomu,
                                Rcpp::Named("recIRdistrRFomu") = recIRdistrRFomu);

  List recADDlist = List::create(Rcpp::Named("recMonoFomu") = recMonoFomu,
                                 Rcpp::Named("recMonoMainAddFomu") = monoNewMainAddFomu,
                                 Rcpp::Named("recMonoFomuWithInSameComp") = recMonoFomuWithInSameComp,
                                 Rcpp::Named("recMITwithInSameComp") = recMITwithInSameComp);

  List recBIOlist = List::create(Rcpp::Named("recCompBioLink") = recCompBioLink,
                                 Rcpp::Named("recCompMonoFomu") = recCompMonoFomu,
                                 Rcpp::Named("recCompMainAddFomu") = compNewMainAddFomu);

  List obsRTlist = List::create(Rcpp::Named("obsFomuRT") = obsFomuRT,
                                Rcpp::Named("obsRTfomu") = obsRTfomu);


  List recRTlist = List::create(Rcpp::Named("recFomuRT") = recFomuRT,
                                Rcpp::Named("recRTfomu") = recRTfomu,
                                Rcpp::Named("recFomuRTdistr") = recFomuRTdistr,
                                Rcpp::Named("recRTdistrFomu") = recRTdistrFomu);

  List model = List::create(Rcpp::Named("pkFomu") = wrap(pkFomu_rvec),
                            Rcpp::Named("pkComp") = wrap(pkComp_rvec),
                            Rcpp::Named("ma") = ma,
                            Rcpp::Named("iso") = iso_spMat,
                            Rcpp::Named("add") = add_spMat,
                            Rcpp::Named("bio") = bio_spMat);

  List record = List::create(Rcpp::Named("ma") = recMAlist,
                             Rcpp::Named("iso") = recIRlist,
                             Rcpp::Named("add") = recADDlist,
                             Rcpp::Named("bio") = recBIOlist,
                             Rcpp::Named("rt") = recRTlist);

  List other = List::create(Rcpp::Named("ma") = obsMAlist,
                            Rcpp::Named("iso") = obsIRlist,
                            Rcpp::Named("rt") = obsRTlist);

  return List::create(Rcpp::Named("model") = model,
                      Rcpp::Named("record") = record,
                      Rcpp::Named("other") = other);

  // return R_NilValue;
}

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

// // [[Rcpp::export]]
// arma::rowvec testAbs(arma::rowvec x){
//   return arma::abs(x);
// }
// // [[Rcpp::export]]
// arma::uvec testUvecFind(NumericVector x, int y){
//   arma::uvec x_uvec = as<arma::uvec>(x);
//   arma::uvec result = arma::find(x_uvec == y);
//   return result;
// }

// // [[Rcpp::export]]
// NumericVector UpdateCompPriorKnowledge(arma::rowvec pkComp_rvec, NumericVector idx, double scale){
//   arma::uvec idx_uvec = as<arma::uvec>(idx);
//   pkComp_rvec.elem(idx_uvec).operator+=(scale);
//   return wrap(pkComp_rvec);
// }


// NumericVector getIdxOfSubVec(arma::rowvec x, arma::uvec x_sub){
//   int n = x_sub.n_elem;
//   NumericVector result;
//   for(int i = 0; i < n; i++){
//     int sub_elem = x_sub.at(i);
//     arma::uvec findSub_uvec = arma::find(x == sub_elem);
//     result.push_back(findSub_uvec.at(0));
//   }
//   
//   return result;
// }

// // [[Rcpp::export]]
// arma::sp_mat testAlterSpMat(arma::mat x){
//   arma::sp_mat x_spMat = arma::sp_mat(x);
//   x_spMat.at(0,0) = 0;
//   x_spMat.at(0,1) = 99;
//   return x_spMat;
// }

