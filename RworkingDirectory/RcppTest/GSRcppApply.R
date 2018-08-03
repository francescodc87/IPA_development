# function to sample from a defined uniform probability distribution
"RcppAsample" <- function(P){
  o <- Rcppsample(1:length(P), size = 1, replace = TRUE, P)
  o
}