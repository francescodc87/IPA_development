# function to sample from a defined uniform probability distribution
"multsampleRcpp" <- function(P){
  o <- sampleRcpp(1:length(P), size = 1, replace = TRUE, P)
  o
}