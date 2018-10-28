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

double UpdateIR(double old_iso,
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
  Rcout << "Type 1 to choose manually, caution if choose to select mass manually you can not ask for software help anymore" << std::endl;
  Rcout << "Type 2 to choose software help, the rule of filtering overlapped assignments are: " << std::endl;
  Rcout << "     I. prefer assignments in which the masses show higher intensities" << std::endl;
  Rcout << "     II. prefer assignments in which detecting lower massaccuracies" << std::endl;
  Rcout << "     III. prefer assignments in which more isotope connections are found based on the mass retention times" << std::endl;
  Rcout << "     IV. prefer assignments in which more adduct connections are found based on the mass retention times" << std::endl;
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

int MsgAndChoiceOfMasses(int fomuNum, arma::rowvec massNum, arma::rowvec mz, arma::rowvec rt){
  Rcout << " " << std::endl;
  Rcout << "----------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  Rcout << "Here are all the information of masses in an overlapped assignment, please look carefully and then choose the correct assignment by type the mass number" << std::endl;
  Rcout << "Once the input number is inserted, other assignments related to the same chemical fomula will be deleted" << std::endl;
  Rcout << "However, if you can not decide which one is the correct assignment, type 0, then all assignments related to this chemical fomula will be deleted" << std::endl;
  Rcout << "    Fomula number: " << fomuNum << std::endl;
  for(int i = 0; i < massNum.n_elem; i++){
    Rcout << "    Mass number: " << massNum.at(i) + 1 << ";  m/z: " << mz.at(i) << ";  RT: " << rt.at(i) << std::endl;
    // Rcout << "    Mass number: " << massNum.at(i) + 1 << std::endl;
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

// arma::uvec SelectLinkWinner(NumericVector isoLink, NumericVector addLink){
//   Rcout << "isoLink" << std::endl;
//   print(isoLink);
//   Rcout << "addLink" << std::endl;
//   print(addLink);
//   NumericVector allLink = isoLink + addLink;
//   arma::rowvec allLink_rvec = as<arma::rowvec>(allLink);
//   double maxNum = arma::max(allLink_rvec);
//   return arma::find(allLink_rvec == maxNum);
// }


void printOutRvec(arma::rowvec x){
  // Rcout << "UVector" << std::endl;
  for (int i = 0; i < x.n_elem; i++){
    Rcout << "idx: " << i << "; value: " << x.at(i) <<std::endl;
  }
}

void printOutUvec(arma::uvec x){
  Rcout << "UVector" << std::endl;
  for (int i = 0; i < x.n_elem; i++){
    Rcout << "idx: " << i << "; value: " << x.at(i) <<std::endl;
  }
}

// [[Rcpp::export]]
List UpdateMain(NumericVector colNameIdx,
                arma::sp_mat post_spMat,                // posterior matrix
                arma::sp_mat iso_spMat,                 // iso matrix
                arma::sp_mat add_spMat,                 // add matrix
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
                NumericVector obsIRmonoFomu,            // storage of all the mono fomula numbers of all the stored observed intensity ratios in this update time
                NumericVector obsIRisoFomu,             // storage of all the iso fomula numbers of all the stored observed intensity ratios in this update time
                List recIRdistr,                        // record of all the distributions of all recorded intensity ratios
                NumericVector recIRdistrMonoFomu,       // record of all the mono fomula numbers of all recorded distributions
                NumericVector recIRdistrIsoFomu,        // record of all the iso fomula numbers of all recorded distributions
                NumericVector recMonoFomu,              // record of all the fomula numbers of all monoisotopic adducts
                NumericVector recAllAddFomu,            // record of all the adducts
                NumericVector recAllAddComp,            // record of all the compounds of according adducts
                NumericVector recComp,                  // compound id
                NumericVector recCompMainAddFomu,       // fomula numbers of the main adducts of compounds
                List recCompMonoFomu,                   // the monoisotopic adduct fomula numbers of compounds
                List recMITwithInSameComp,              // the most intense time values of monoisotopic adducts of each compound (share the same compound)
                NumericVector massRT,                   // mass retention time values of all the masses
                NumericVector obsCompRT,                // storage of observed retention time values of fomulas in this update time
                NumericVector obsRTcomp,                // the fomula numbers of the observed retention time in this update time
                List recCompRT,                         // record all observed retention time values of fomulas
                NumericVector recRTcomp,                // the fomula numbers of the recorded retention time values
                List recCompRTdistr,                    // record all the distributions of retention time values of fomulas
                NumericVector recRTdistrComp,           // the fomula numbers of the recorded retention time distributions
                double s_pkAdd,                         // update scale(speed) for pkFomu
                double s_pkComp,                        // update scale(speed) for pkComp
                double s_ma,                            // update scale for ma (0< <=1)
                double s_iso,                           // update scale for iso (0< <=1)
                double w_int,                           // weight for most intense mass when solving overlapped assignment
                double w_ma,                            // weight for min mass accuracy mass when solving overlapped assignment
                double w_addLink,                       // weight for adduct connection of mass when solving overlapped assignment
                double w_isoLink,                       // weight for isotope connection of of mass when solving overlapped assignment
                double t_RT,                            // threshold used to count iso and add connection when solving overlapped assignments automatically
                double t_post,                          // threshold used to filter posterior
                bool isTest = false                     // is it a test?

){
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>>>START THE WHOLE UPDATE STAGE" << std::endl;
  Rcout << " " << std::endl;

  NumericVector massSet;                                                  // mass number set after threshold
  NumericVector fomuSet;                                                  // fomula number set after threshold
  NumericVector fomuIdSet;                                                // formula id in database set after threshold   
  // NumericVector fomuCompIdSet;                                         // fomula id set after threshold
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
      fomuIdSet.push_back(colNameIdx.at(colNum + 1));
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
  arma::rowvec fomuIdSet_rvec = as<arma::rowvec>(fomuIdSet);
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
    int firstChoice = MsgAndChoWhenOverlapAssignDetected(dupFomuUni_uvec.n_elem);
    if(firstChoice == 0){
      MsgWhenIgnoreAllOverlapAssign();
      for(int i = 0; i < dupFomuUni_uvec.n_elem; i++){
        int dupFomuNum = dupFomuUni_uvec.at(i);                                     // first get one of the duplicate fomula number
        arma::uvec dupFomuIdx_uvec = arma::find(fomuSet_rvec == dupFomuNum);        // search the idx of the duplicate fomula number in fomula row vector
        arma::uvec allFomuIdx_uvec = find_finite(fomuSet_rvec);                     // all idx in fomula row vector
        arma::uvec retainFomuIdx_uvec = SetDiff(allFomuIdx_uvec, dupFomuIdx_uvec);  // the idx of retained fomula (which is also the idx of the retained masses)
        massSet_rvec = massSet_rvec.elem(retainFomuIdx_uvec).t();                   // the retained mass
        fomuSet_rvec = fomuSet_rvec.elem(retainFomuIdx_uvec).t();                   // the retained formula number
        fomuIdSet_rvec = fomuIdSet_rvec.elem(retainFomuIdx_uvec).t();               // the retained formula id
        // Rcout<< "0;massSet_rvec" << std::endl;
        // printOutRvec(massSet_rvec);
        // Rcout<< "0;fomuSet_rvec" << std::endl;
        // printOutRvec(fomuSet_rvec);
      }
    } else if(firstChoice == 1){
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
      for(int i = 0; i < dupFomuUni_uvec.n_elem; i++){
        int dupFomuNum = dupFomuUni_uvec.at(i);                                     // first get one of the duplicate fomula number
        arma::uvec dupFomuIdx_uvec = arma::find(fomuSet_rvec == dupFomuNum);        // search the idx of the duplicate fomula number in fomula row vector
        arma::rowvec dupMasses_rvec = massSet_rvec.elem(dupFomuIdx_uvec).t();       // get the duplicate masses
        arma::uvec dupMasses_uvec = as<arma::uvec>(wrap(dupMasses_rvec));
        // Rcout<< "dupMasses_uvec" << std::endl;
        // printOutUvec(dupMasses_uvec);
        arma::rowvec dupMassesMZ_rvec = massMZ_rvec.elem(dupMasses_uvec).t();       // get the mass accuracy of the duplicate masses
        arma::rowvec dupMassesRT_rvec = massRT_rvec.elem(dupMasses_uvec).t();       // get the retention time of the duplicate masses
        // Rcout<< "dupMassesMZ_rvec" << std::endl;
        // printOutRvec(dupMassesMZ_rvec);
        // Rcout<< "dupMassesRT_rvec" << std::endl;
        // printOutRvec(dupMassesRT_rvec);
        int select_mass = MsgAndChoiceOfMasses(colNameIdx.at(dupFomuNum), dupMasses_rvec, dupMassesMZ_rvec, dupMassesRT_rvec);
        if (select_mass != -1){
          arma::uvec selectMassIdx_uvec = arma::find(dupMasses_rvec == select_mass);
          selectedIdx.push_back(dupFomuIdx_uvec.at(selectMassIdx_uvec.at(0)));        // push back the index 
          MsgAfterChoiceOfMasses();   // print some information after user making the choice each time
        }
      }
      arma::uvec selectedIdx_uvec = as<arma::uvec>(selectedIdx);
      massSet_rvec = massSet_rvec.elem(selectedIdx_uvec).t();                   // the retained mass
      fomuSet_rvec = fomuSet_rvec.elem(selectedIdx_uvec).t();                   // the retained formula num
      fomuIdSet_rvec = fomuIdSet_rvec.elem(selectedIdx_uvec).t();               // the retained formula id
      // Rcout<< "1;massSet_rvec" << std::endl;
      // printOutRvec(massSet_rvec);
      // Rcout<< "1;fomuSet_rvec" << std::endl;
      // printOutRvec(fomuSet_rvec);
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
        int dupFomuNum = dupFomuUni_uvec.at(i);                                         // first get one of the duplicate fomula number
        arma::uvec dupFomuIdx_uvec = arma::find(fomuSet_rvec == dupFomuNum);            // search the idx of the duplicate fomula number in fomula row vector
        arma::rowvec dupSubMassSet_rvec = massSet_rvec.elem(dupFomuIdx_uvec).t();       // get duplicated masses
        
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
          dupSubIsoLink.push_back(isoLink);
          dupSubAddLink.push_back(addLink);
        }

        arma::rowvec dupSubMassIntSet_rvec = as<arma::rowvec>(dupSubMassIntSet);
        arma::rowvec dupSubMAset_rvec = as<arma::rowvec>(dupSubMAset);
        // Rcout<< "dupSubMassSet_rvec" <<std::endl;
        // printOutRvec(dupSubMassSet_rvec);
        // rule I: prefer assignments in which the masses show higher intensities
        double maxInt = arma::max(dupSubMassIntSet_rvec);
        arma::uvec dupSubMaxIntIdx = arma::find(dupSubMassIntSet_rvec == maxInt);
        dupSubScore_rvec.elem(dupSubMaxIntIdx) += w_int;       // score add
        // Rcout<< "dupSubMaxIntIdx" <<std::endl;
        // printOutUvec(dupSubMaxIntIdx);
        
        // rule II: prefer assignments in which detecting lower massaccuracies
        double minMA = arma::min(dupSubMAset_rvec);
        arma::uvec dupSubMinMAidx = arma::find(dupSubMAset_rvec == minMA);
        dupSubScore_rvec.elem(dupSubMinMAidx) += w_ma;       // score add
        // Rcout<< "dupSubMinMAidx" <<std::endl;
        // printOutUvec(dupSubMinMAidx);

        // rule III: prefer assignments in which more isotope connections are found based on the mass retention times
        arma::rowvec dupSubIsoLink_rvec = as<arma::rowvec>(dupSubIsoLink);
        double maxIsoLink = arma::max(dupSubIsoLink_rvec);
        arma::uvec dupSubMaxIsoLink = arma::find(dupSubIsoLink_rvec == maxIsoLink);
        dupSubScore_rvec.elem(dupSubMaxIsoLink) += w_isoLink;    // score add
        // Rcout<< "dupSubMaxIsoLink" <<std::endl;
        // printOutUvec(dupSubMaxIsoLink);
        
        // rule IV: prefer assignments in which more adduct connections are found based on the mass retention times
        arma::rowvec dupSubAddLink_rvec = as<arma::rowvec>(dupSubAddLink);
        double maxAddLink = arma::max(dupSubAddLink_rvec);
        arma::uvec dupSubMaxAddLink = arma::find(dupSubAddLink_rvec == maxAddLink);
        dupSubScore_rvec.elem(dupSubMaxAddLink) += w_addLink;       // score add
        // Rcout<< "dupSubMaxAddLink" <<std::endl;
        // printOutUvec(dupSubMaxAddLink);
        
        // after thresholded by three rules, perform automatic selection of mass in the subvec of duplicate assignments by the final score
        // dupSubScore_rvec, dupSubFomuSet_rvec, dupSubMassSet_rvec, dupFomuIdx_uvec
        int maxScore = arma::max(dupSubScore_rvec);
        arma::uvec maxScoreFind_uvec = arma::find(dupSubScore_rvec == maxScore);
        if(maxScoreFind_uvec.n_elem == 1){
          selectedIdx.push_back(dupFomuIdx_uvec.at(maxScoreFind_uvec.at(0)));
        } else{
          int select_mass = -1;
          arma::rowvec manualDupMass_rvec = dupSubMassSet_rvec.elem(maxScoreFind_uvec).t();
          arma::rowvec manualDupScore_rvec = dupSubScore_rvec.elem(maxScoreFind_uvec).t();
          arma::uvec manualDupFomuIdx_uvec = dupFomuIdx_uvec.elem(maxScoreFind_uvec);
          arma::uvec manualDupMass_uvec = as<arma::uvec>(wrap(manualDupMass_rvec));
          arma::rowvec manualDupMassMZ_rvec = massMZ_rvec.elem(manualDupMass_uvec).t();       // get the mass accuracy of the duplicate masses
          arma::rowvec manualDupMassRT_rvec = massRT_rvec.elem(manualDupMass_uvec).t();       // get the retention time of the duplicate masses
          if(!isTest){
            select_mass = MsgAndChoiceOfMasses(colNameIdx.at(dupFomuNum), manualDupMass_rvec, manualDupMassMZ_rvec, manualDupMassRT_rvec);
          }
          if(select_mass != -1){                                                                              
            arma::uvec selectMassFind_uvec = arma::find(manualDupMass_rvec == select_mass);
            selectedIdx.push_back(dupFomuIdx_uvec.at(selectMassFind_uvec.at(0)));
            MsgAfterChoiceOfMasses();   // print some information after user making the choice each time                                                                 
          } 
        }
      }
      arma::uvec selectedIdx_uvec = as<arma::uvec>(selectedIdx);
      fomuSet_rvec = fomuSet_rvec.elem(selectedIdx_uvec).t();                       // filtered fomula num set
      fomuIdSet_rvec = fomuIdSet_rvec.elem(selectedIdx_uvec).t();                   // the retained formula id
      massSet_rvec = massSet_rvec.elem(selectedIdx_uvec).t();                       // filtered mass set
      // Rcout<< "2;massSet_rvec" << std::endl;
      // printOutRvec(massSet_rvec);
      // Rcout<< "2;fomuSet_rvec" << std::endl;
      // printOutRvec(fomuSet_rvec);
    }
  }

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>OVERLAPPED ASSIGNMENTS ISSUE FIXED, CONTINUE THE UPDATE STAGE" << std::endl;
  Rcout << " " << std::endl;
  // based on assignment pair after thresholding posterior probability matrix:
  // which are: mass numbers and fomula numbers,
  // do: update prior knowledge,
  //     record all observed mass accuracies in this update time before updating mass accuracy
  //     and record all observed intensity ratios of observed isotope pairs in this update time before updating iso
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START RECORD OBSERVED MASS ACCURACIES, RETENTION TIME AND INTENSITY RATIOS;" << std::endl;
  Rcout << "                       UPDATE MOST INTENSE TIME TO PREPARE MAIN ADDUCT SHUFFLE, AND UPDATE FORMULA PRIOR KNOWLEDGE" << std::endl;
  Rcout << " " << std::endl;
  arma::rowvec pkFomu_rvec = as<arma::rowvec>(pkFomu);                                                  // prior knowledge of all the fomulas
  arma::rowvec obsIRmonoFomu_rvec = as<arma::rowvec>(obsIRmonoFomu);                                    // mono fomula numbers of observed intensity values
  arma::rowvec obsIRisoFomu_rvec = as<arma::rowvec>(obsIRisoFomu);                                      // iso fomula numbers of observed intensity values
  arma::uvec fomuIdSet_uvec = as<arma::uvec>(wrap(fomuIdSet_rvec));
  arma::rowvec recMonoFomu_rvec = as<arma::rowvec>(recMonoFomu);  
  arma::rowvec recAllAddFomu_rvec = as<arma::rowvec>(recAllAddFomu);
  NumericVector obsRTvalues;
  NumericVector obsRTcompIds;
  
  // update compound prior knowledge: 1st part
  arma::uvec compMainAdd_uvec = as<arma::uvec>(recCompMainAddFomu);                                     // fomula numbers of all the compound main adducts
  arma::uvec compMainAddIntsc_uvec = arma::intersect(compMainAdd_uvec, fomuIdSet_uvec);                  // we want to see how many compound main adducts survive after thresholding the posterior matrix
  arma::rowvec compPriorMark_rvec = arma::zeros<arma::rowvec>(compNum);
  if(!compMainAddIntsc_uvec.is_empty()){
    // if any compound main adducts survive, set a mark
    for(int i = 0; i < compMainAddIntsc_uvec.n_elem; i++){
      int compMainAddFomuId = compMainAddIntsc_uvec.at(i);
      arma::uvec findcompMain_uvec = arma::find(compMainAdd_uvec == compMainAddFomuId);
      compPriorMark_rvec.at(findcompMain_uvec.at(0)) += 1;
    }
  }
  for(int i = 0; i < fomuSet_rvec.n_elem; i++){                                                            // iterate on thresholded fomula numbers
    // after detect and fix overlap assignments, Update fomula Prior Knowledge
    int theLeftFomuNum = fomuSet_rvec.at(i);                                                                
    int theLeftMassNum = massSet_rvec.at(i);                                                                 
    int theLeftFomuId = fomuIdSet_rvec.at(i);
    pkFomu_rvec.at(theLeftFomuNum) = UpdateAddPriorKnowledge(pkFomu_rvec.at(theLeftFomuNum), s_pkAdd);     // update fomula prior knowledge
    
    // calculate and record observed mass accuracies (negative values allowed)
    double maValue = ((massMZ[theLeftMassNum] - fomuMZ[theLeftFomuNum]) * 1e6 ) / fomuMZ[theLeftFomuNum];
    if (maValue >= 50){
      Rcout << "theLeftMassNum: " << theLeftMassNum << "; theLeftFomuNum: " << theLeftFomuNum << std::endl;
    }
    obsMA.push_back(maValue);
    recMA.push_back(maValue);
    
    // record observed retention times
    arma::uvec findAddId_uvec = arma::find(recAllAddFomu_rvec == theLeftFomuId);
    if(!findAddId_uvec.is_empty()){
      obsRTcompIds.push_back(recAllAddComp.at(findAddId_uvec.at(0)));
      obsRTvalues.push_back(massRT(theLeftMassNum));
    }
    
    // record the observed intensity ratios before updating iso matrix, we use another iteration of fomuSet_rvec to create isotope pairs
    for(int j = 0; j < fomuSet_rvec.n_elem; j++){
      // generate the right part of fomula numbers which will be used as indices to record intensity ratios
      // and generate the right part of masses which will be used as indices to calculate observed intensity ratios
      int theRightFomuNum = fomuSet_rvec.at(j);
      int theRightMassNum = massSet_rvec.at(j);
      int theRightFomuId = fomuIdSet_rvec.at(j);

      // next is to combine the fomula number pair as isotope pair, use mass number pair to compute isotope ratio
      // and complete the storage task before update
      // about the next line judgement code: here we decrease the computational complexity in the second iteration, actually I iterate half, as you can see from the condition in the next if() statement
      //                                     & the condition also considerds that only the fomula which is an isotope pair can be applied in the following task
      if((iso_spMat.at(theLeftFomuNum, theRightFomuNum) != 0)){
        arma::uvec findmonoAdd = arma::find(recMonoFomu_rvec == theLeftFomuId);
        if (findmonoAdd.is_empty()){
          double theIRvalue = massInt_rvec.at(theRightMassNum) / massInt_rvec.at(theLeftMassNum);              // isotope ratio
          // we want to know if we have observed any isotope ratios of this fomula pair before
          // so we take the intersect
          arma::uvec obsIRmonoFomuFind_uvec = arma::find(obsIRmonoFomu_rvec == theRightFomuId);
          arma::uvec obsIRisoFomuFind_uvec = arma::find(obsIRisoFomu_rvec == theLeftFomuId);
          arma::uvec obsIRisoFind_uvec = arma::intersect(obsIRmonoFomuFind_uvec, obsIRisoFomuFind_uvec);
          // if we have not recorded the observed isotope ratio of this fomula number pair before, we push_back those value and indices
          // else we pop the value into the same position we find it, but we don't pop the indices in this situation cause the indices are already existed
          if(obsIRisoFind_uvec.is_empty()){
            NumericVector obsIRSub;
            obsIRSub.push_back(theIRvalue);
            obsIR.push_back(theIRvalue);
            obsIRmonoFomu.push_back(theRightFomuId);
            obsIRisoFomu.push_back(theLeftFomuId);
          } else{            // already spoken situation that fomula number pair is already recorded
            int obsIRisoFind = obsIRisoFind_uvec.at(0);
            NumericVector obsIRSub = obsIR[obsIRisoFind];
            obsIRSub.push_back(theIRvalue);
            obsIR[obsIRisoFind] = obsIRSub;
          }
        } else{
          double theIRvalue = massInt_rvec.at(theLeftMassNum) / massInt_rvec.at(theRightMassNum);              // isotope ratio
          // we want to know if we have recorded any isotope ratios of this fomula pair before
          // so we take the intersect
          arma::uvec obsIRmonoFomuFind_uvec = arma::find(obsIRmonoFomu_rvec == theLeftFomuId);
          arma::uvec obsIRisoFomuFind_uvec = arma::find(obsIRisoFomu_rvec == theRightFomuId);
          arma::uvec obsIRisoFind_uvec = arma::intersect(obsIRmonoFomuFind_uvec, obsIRisoFomuFind_uvec);
          // if we have not recorded the observed isotope ratio of this fomula number pair before, we push_back those value and indices
          // else we pop the value into the same position we find it, but we don't pop the indices in this situation cause the indices are already existed
          if(obsIRisoFind_uvec.is_empty()){
            NumericVector obsIRSub;
            obsIRSub.push_back(theIRvalue);
            obsIR.push_back(obsIRSub);
            obsIRmonoFomu.push_back(theLeftFomuId);
            obsIRisoFomu.push_back(theRightFomuId);
          } else{           // already spoken situation that fomula number pair is already recorded
            int obsIRisoFind = obsIRisoFind_uvec.at(0);
            NumericVector obsIRSub = obsIR[obsIRisoFind];
            obsIRSub.push_back(theIRvalue);
            obsIR[obsIRisoFind] = obsIRSub;
          }
        }
      }
    }
  }

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>UPDATE OBSERVED MOST INTENSE TIME VALUES OF MONOISOTOPIC ADDUCTS" << std::endl;
  Rcout << " " << std::endl;
  // update recorded most intense time values of each monoisotopic adduct
  for(int i = 0; i < recComp.length(); i++){
    // get all other adducts connected with this monoisotopic adduct (linked by the same compound)
    NumericVector recMonoWithInSameComp = recCompMonoFomu.at(i);                            // get all adducts except isotopes with in this compound
    arma::uvec recMonoWithInSameComp_uvec = as<arma::uvec>(recMonoWithInSameComp);

    // next: decide which one is the most intense, and update the most intense time in recMITwithInSameComp
    NumericVector monoAddInt;
    NumericVector monoAddIdx;
    arma::uvec monoAddAllIdx = find_finite(recMonoWithInSameComp_uvec.t());
    for(int j = 0; j < recMonoWithInSameComp_uvec.n_elem; j++){                            // iterate on all monoisotopic adducts
      int monoAddFomuId = recMonoWithInSameComp_uvec.at(j);
      arma::uvec findMonoAddId_uvec = arma::find(fomuIdSet_rvec == monoAddFomuId);
      if(!findMonoAddId_uvec.is_empty()){
        int monoAddIdxInFomu = findMonoAddId_uvec.at(0);
        int monoAddMassNum = massSet_rvec.at(monoAddIdxInFomu);
        monoAddInt.push_back(massInt_rvec.at(monoAddMassNum));
        monoAddIdx.push_back(monoAddAllIdx.at(j));
      } 
    }
    // based on monoAddInt and monoAddIdx, update recMITwithInSameComp
    if(monoAddIdx.length() != 0){
      arma::rowvec monoAddInt_rvec = as<arma::rowvec>(monoAddInt);
      recMITwithInSameComp.at(i) = UpdateAddMostIntTime(monoAddIdx, monoAddInt_rvec, recMITwithInSameComp.at(i));
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
    if(obsMA_rvec.n_elem >= 1){                                                  // if have observed one mass accuracy value in this update time, then we can update the distribution of mass accuracy
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
  arma::rowvec obsNewIRmonoFomu_rvec = as<arma::rowvec>(obsIRmonoFomu);            // mono part of the observed isotope pairs
  arma::rowvec obsNewIRisoFomu_rvec = as<arma::rowvec>(obsIRisoFomu);              // iso part of the observed isotope pairs
  arma::rowvec recIRdistrMonoFomu_rvec = as<arma::rowvec>(recIRdistrMonoFomu);     // mono part of the recorded isotope pairs (of the recorded distribution information)
  arma::rowvec recIRdistrIsoFomu_rvec = as<arma::rowvec>(recIRdistrIsoFomu);       // iso part of the recorded isotope pairs (of the recorded distribution information)
  NumericVector newIRmonoFomu;
  NumericVector newIRisoFomu;
  NumericVector newIRirValue;
  IntegerVector remIRidx;
  // iterate on all observed isotope pairs,
  // the idea is to find out if the observed isotope pair is already found in the isotope pairs of recorded distribution information about intensity ratios
  for(int i = 0; i < obsNewIRmonoFomu_rvec.n_elem; i++){
    int obsIRmonoFomuId = obsNewIRmonoFomu_rvec.at(i);
    int obsIRisoFomuId = obsNewIRisoFomu_rvec.at(i);
    NumericVector theIRvalue = obsIR[i];                                    // the intensity ratio
    arma::uvec recIRdistrMonoFomuFind_uvec = arma::find(recIRdistrMonoFomu_rvec == obsIRmonoFomuId);
    arma::uvec recIRdistrIsoFomuFind_uvec = arma::find(recIRdistrIsoFomu_rvec == obsIRisoFomuId);
    arma::uvec recIRdistrIsoFind_uvec = arma::intersect(recIRdistrMonoFomuFind_uvec, recIRdistrIsoFomuFind_uvec);
    // if we have not recorded the distribution information of this observed mono-iso pair before,
    // we can push_back the first distribution information of this observed mono-iso pair and the indices
    // only if observed more than 2 intensity ratios
    if(recIRdistrIsoFind_uvec.n_elem == 0){ 
      if(theIRvalue.length() > 1){
        List theIRvalueDistr;
        NumericVector theIRlogValue = log(theIRvalue);
        arma::rowvec theIRlogValue_rvec = as<arma::rowvec>(theIRlogValue);
        // calculate IR numbers, mean and std value for this isotope pair for the first time
        double theIRvalue_mean = GeoMean(theIRvalue.at(0), theIRvalue.at(1));
        double theIRvalue_logMean = mean(theIRlogValue);
        double theIRvalue_std = GeoStd(theIRvalue);
        double theIRvalue_logStd = arma::stddev(theIRlogValue_rvec);
        double theIRvalue_num = 2;
        theIRvalueDistr["iso_num"] = theIRvalue_num;
        theIRvalueDistr["iso_mean"] = theIRvalue_mean;
        theIRvalueDistr["iso_logMean"] = theIRvalue_logMean;
        theIRvalueDistr["iso_std"] = theIRvalue_std;
        theIRvalueDistr["iso_logStd"] = theIRvalue_logStd;
        recIRdistr.push_back(theIRvalueDistr);
        recIRdistrMonoFomu.push_back(obsIRmonoFomuId);
        recIRdistrIsoFomu.push_back(obsIRisoFomuId);

        // record the mono-iso pairs in which the intensity ratio needs update
        newIRmonoFomu.push_back(obsIRmonoFomuId);
        newIRisoFomu.push_back(obsIRisoFomuId);
        arma::uvec findMonoFomuId_uvec = arma::find(fomuIdSet_rvec == obsIRmonoFomuId);
        arma::uvec findIsoFomuId_uvec = arma::find(fomuIdSet_rvec == obsIRisoFomuId);
        double oldIRvalue = iso_spMat.at(fomuSet_rvec.at(findMonoFomuId_uvec.at(0)), fomuSet_rvec.at(findIsoFomuId_uvec.at(0)));
        newIRirValue.push_back(UpdateIR(oldIRvalue, theIRvalue_mean, s_iso));
        remIRidx.push_back(i);
      }
    } else{
      // if we have recorded the distribution information of this observed mono-iso pair before,
      // we pop the distribution information into the same position we find it,
      // but we don't pop the indices in this situation cause the indices are already existed
      int recIRdistrIsoIdx = recIRdistrIsoFind_uvec.at(0);
      double theIRvalue_mean = theIRvalue.at(0);

      // update distribution information
      recIRdistr[recIRdistrIsoIdx] = UpdateIsoInfo(recIRdistr[recIRdistrIsoIdx], theIRvalue_mean, log(theIRvalue_mean), 0, 1);

      // get new intensity ratio mean value
      List theIRvalueDistr = recIRdistr[recIRdistrIsoIdx];
      theIRvalue_mean = theIRvalueDistr["iso_mean"];

      // record the mono-iso pairs in which the intensity ratio needs update
      newIRmonoFomu.push_back(obsIRmonoFomuId);
      newIRisoFomu.push_back(obsIRisoFomuId);
      arma::uvec findMonoFomuId_uvec = arma::find(fomuIdSet_rvec == obsIRmonoFomuId);
      arma::uvec findIsoFomuId_uvec = arma::find(fomuIdSet_rvec == obsIRisoFomuId);
      double oldIRvalue = iso_spMat.at(fomuSet_rvec.at(findMonoFomuId_uvec.at(0)), fomuSet_rvec.at(findIsoFomuId_uvec.at(0)));
      newIRirValue.push_back(UpdateIR(oldIRvalue, theIRvalue_mean, s_iso));
      remIRidx.push_back(i);
    }
  }

  // after update iso, we record the mono-iso pairs in which the observed data needs to be wiped. and record the mono-iso pairs in which the observed data needs to be updated
  IntegerVector obsIRidx = seq_len(obsNewIRmonoFomu_rvec.n_elem) - 1;
  arma::uvec obsIRidx_uvec = as<arma::uvec>(obsIRidx);
  arma::uvec remIRidx_uvec = as<arma::uvec>(remIRidx);
  arma::uvec retainIRidx_uvec = SetDiff(obsIRidx_uvec, remIRidx_uvec);
  IntegerVector retainIRidx = wrap(retainIRidx_uvec.t());

  // next we shuffle the main adduct relations based on the observed most intense time of monoisotopic adducts
  // which contanis main adduct relations between compounds and mono adducts
  // recCompMainAddFomu, recCompMonoFomu recMITwithInSameComp
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO SHUFFLE MAIN ADDUCT RELATIONS" << std::endl;
  Rcout << " " << std::endl;
  NumericVector newMainAddCompIdx;
  NumericVector newMainAddId;
  for(int i = 0; i < recCompMainAddFomu.length(); i++){
    // based on the recorded most intense adduct time of each mono-isotopic adduct linked with the same compound, update main adduct information
    int compOldMainAdd = recCompMainAddFomu.at(i);
    int compNewMainAdd = CompareToShuffleMainAdd(compOldMainAdd, recCompMonoFomu.at(i), recMITwithInSameComp.at(i));
    if (compOldMainAdd != compNewMainAdd){
      newMainAddCompIdx.push_back(i);
      newMainAddId.push_back(compNewMainAdd);
      compPriorMark_rvec.at(i) += 1;
    }
  }
  
  // after shuffle the main adducts, we can start update compound prior knowledge now
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO UPDATE COMPOUND PRIOR KNOWLEDGE" << std::endl;
  Rcout << " " << std::endl;
  // compPriorMark_rvec
  arma::rowvec pkComp_rvec = as<arma::rowvec>(pkComp);                                                  // prior knowledge of all the compounds
  arma::uvec compNeedPriorUpdate_uvec = arma::find(compPriorMark_rvec != 0);
  pkComp_rvec = UpdateCompPriorKnowledge(pkComp_rvec, compNeedPriorUpdate_uvec, s_pkComp);
  
  // the last one: update retention time information of compounds
  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>START TO RECORD RETENTION TIME DISTRIBUTION" << std::endl;
  Rcout << " " << std::endl;
  // obsRTcompIds obsRTvalues
  // obsCompRT, obsRTcomp, recCompRT, recRTcomp, recCompRTdistr, recRTdistrComp
  // first integrate the observed retention time information
  arma::rowvec obsRTcompIds_rvec = as<arma::rowvec>(obsRTcompIds);
  arma::rowvec obsRTcompIdUnique_rvec = unique(obsRTcompIds_rvec);
  arma::rowvec obsRTvalues_rvec = as<arma::rowvec>(obsRTvalues);
  arma::rowvec obsRTcomp_rvec = as<arma::rowvec>(obsRTcomp);                                      
  arma::rowvec recRTcomp_rvec = as<arma::rowvec>(recRTcomp);
  arma::rowvec recRTdistrComp_rvec = as<arma::rowvec>(recRTdistrComp);
  NumericVector remObsRTidx;
  for(int i = 0; i < obsRTcompIdUnique_rvec.n_elem; i++){
    int obsRTcompId = obsRTcompIdUnique_rvec.at(i);
    arma::uvec findComp_uvec = arma::find(obsRTcompIds_rvec == obsRTcompId);
    arma::rowvec obsRTset_rvec = obsRTvalues_rvec.elem(findComp_uvec).t();
    double obsRTvalue = arma::mean(obsRTset_rvec);
    // check if we recorded retention time info of this compound before
    arma::uvec recRTcompFind_uvec = arma::find(recRTcomp_rvec == obsRTcompId);
    if(recRTcompFind_uvec.is_empty()){    // if not recorded before, we check if observed before
      arma::uvec obsRTcompFind_uvec = arma::find(obsRTcomp_rvec == obsRTcompId);
      if(obsRTcompFind_uvec.is_empty()){  // if not observed before as well
        obsCompRT.push_back(obsRTvalue);
        obsRTcomp.push_back(obsRTcompId);
      } else{                             // if observed before, directly wrap them into distribution
        int obsRTidx = obsRTcompFind_uvec.at(0);
        double obsRToldValue = obsCompRT.at(obsRTidx);
        NumericVector rt;
        rt.push_back(obsRToldValue);
        rt.push_back(obsRTvalue);
        recCompRT.push_back(rt);
        recRTcomp.push_back(obsRTcompId);
        arma::rowvec rt_rvec = as<arma::rowvec>(rt);
        double rt_mean = arma::mean(rt_rvec);
        double rt_std = arma::stddev(rt_rvec);
        double rt_num = 2;
        List rt_distr;
        rt_distr["rt_num"] = rt_num;
        rt_distr["rt_mean"] = rt_mean;
        rt_distr["rt_std"] = rt_std;
        recCompRTdistr.push_back(rt_distr);          // here we push_back the distribution information
        recRTdistrComp.push_back(obsRTcompId);       // here we push_back the index (compound id)
        remObsRTidx.push_back(obsRTidx);
      }
    } else{                                          // if recorded before, which means we do not need to store it into observed format any more, just store it into recorded format and update the distribution
      int recRTidx = recRTcompFind_uvec.at(0);
      NumericVector recRToldValues = recCompRT[recRTidx];
      recRToldValues.push_back(obsRTvalue);
      recCompRT[recRTidx] = recRToldValues;
      double rt_mean = obsRTvalue;
      double rt_std = 0;
      double rt_num = 1;
      // find the idx of this compound in recCompRTdistr(should be the same, just in case)
      arma::uvec findRecRTdistr_uvec = arma::find(recRTdistrComp_rvec == obsRTcompId);
      int recRTdistrIdx = findRecRTdistr_uvec.at(0);
      // update retention time distribution information
      recCompRTdistr[recRTdistrIdx] = UpdateRtInfo(recCompRTdistr[recRTdistrIdx], rt_mean, rt_std, rt_num);
    }
  }

  // after record retention time distribution information, we find those observed data from obsCompRT, obsRTcomp need to remove
  IntegerVector allObsRTidx = seq_len(obsRTcomp.length()) - 1;
  arma::uvec allObsRTidx_uvec = as<arma::uvec>(allObsRTidx);
  arma::uvec remObsRTidx_uvec = as<arma::uvec>(remObsRTidx);
  arma::uvec retainObsRTidx_uvec = SetDiff(allObsRTidx_uvec, remObsRTidx_uvec);
  IntegerVector retainObsRTidx = wrap(retainObsRTidx_uvec.t());

  Rcout << " " << std::endl;
  Rcout << ">>>>>>>>>>>>>>>>>>>>>>>>END OF WHOLE UPDATE STAGE" << std::endl;
  Rcout << " " << std::endl;

  // next is the return value set, what should we return?
  //
  // pk: pkFomu_rvec, pkComp_rvec
  //
  // ma: ma, obsMA, recMAdistr, recMA
  //
  // iso: obsIR, obsIRmonoFomu, obsIRisoFomu, removeIRidx, retainIRidx 
  //      newIRmonoFomu, newIRisoFomu, newIRirValue
  //      recIRdistr, recIRdistrMonoFomu, recIRdistrIsoFomu
  //
  // add: recComp, recCompMonoFomu, recMITwithInSameComp, newMainAddCompIdx, newMainAddId 
  //
  // rt: obsCompRT, obsRTcomp, retainObsRTidx
  //     recCompRT, recRTcomp
  //     recCompRTdistr, recRTdistrComp

  // // normalise prior knowledge
  // pkFomu_rvec = pkFomu_rvec / sum(pkFomu_rvec);
  // pkComp_rvec = pkComp_rvec / sum(pkComp_rvec);

  // return a list, contain all the updated results
  List pkList = List::create(Rcpp::Named("pkFomu") = wrap(pkFomu_rvec),
                             Rcpp::Named("pkComp") = wrap(pkComp_rvec));

  List maList = List::create(Rcpp::Named("ma") = ma,
                             Rcpp::Named("obsMA") = obsMA,
                             Rcpp::Named("recMA") = recMA,
                             Rcpp::Named("recMAdistr") = recMAdistr);
  
  List isoList = List::create(Rcpp::Named("obsIR") = obsIR,
                             Rcpp::Named("obsIRmonoFomu") = obsIRmonoFomu,
                             Rcpp::Named("obsIRisoFomu") = obsIRisoFomu,
                             Rcpp::Named("removeIRidx") = remIRidx + 1,
                             Rcpp::Named("retainIRidx") = retainIRidx + 1,
                             Rcpp::Named("newIRmonoFomu") = newIRmonoFomu,
                             Rcpp::Named("newIRisoFomu") = newIRisoFomu,
                             Rcpp::Named("newIRirValue") = newIRirValue,
                             Rcpp::Named("recIRdistr") = recIRdistr,
                             Rcpp::Named("recIRdistrMonoFomu") = recIRdistrMonoFomu,
                             Rcpp::Named("recIRdistrIsoFomu") = recIRdistrIsoFomu);
  
  List addList = List::create(Rcpp::Named("recComp") = recComp,
                              Rcpp::Named("recCompMonoFomu") = recCompMonoFomu,
                              Rcpp::Named("recMITwithInSameComp") = recMITwithInSameComp,
                              Rcpp::Named("newMainAddCompIdx") = newMainAddCompIdx + 1,
                              Rcpp::Named("newMainAddId") = newMainAddId);
  
  List rtList = List::create(Rcpp::Named("obsCompRT") = obsCompRT,
                              Rcpp::Named("obsRTcomp") = obsRTcomp,
                              Rcpp::Named("retainObsRTidx") = retainObsRTidx + 1,
                              Rcpp::Named("recCompRT") = recCompRT,
                              Rcpp::Named("recRTcomp") = recRTcomp,
                              Rcpp::Named("recCompRTdistr") = recCompRTdistr,
                              Rcpp::Named("recRTdistrComp") = recRTdistrComp);

  return List::create(Rcpp::Named("pk") = pkList,
                      Rcpp::Named("ma") = maList,
                      Rcpp::Named("iso") = isoList,
                      Rcpp::Named("add") = addList,
                      Rcpp::Named("rt") = rtList);

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

