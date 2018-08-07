### using correlations
"CheckingCorr" <- function(corrVec,i, corrThld){
  # no need to record i, because correlation(i,i) = 1, always < corr.thr, so i always included
  out <- which(corrVec<corrThld)
  out
}