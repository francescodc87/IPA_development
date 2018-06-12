### using correlations
"checking.corr" <- function(Corr.vec,i, corr.thr){
  # no need to record i, because correlation(i,i) = 1, always < corr.thr, so i always included
  out <- which(Corr.vec<corr.thr)
  out
}