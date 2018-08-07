# function to compute the posterior probabilities for the allsampcom
# post <- t(apply(allsampcomp,2,ComputePost, burn=burn, it=it, Nc=Nc))
# Nc: row number of adducts connections matrix
"ComputePost" <- function(allsampcomp, it, burn, Nc){
  # table function: assignment distribution count for every column data(every mass) of all sampcomp, the following division turn them into probabilities
  # tmp example:    1    12
  #               1753   247         this means that mass is assigned to compound 1 1753 times, to compound 12 247 times
  tmp <- table(allsampcomp[(burn+1):it])/(it-burn)
  out <- rep(0,Nc)
  # insert tmp value into out vector, as.numeric(names(tmp) assures the value being inserted into correct column number (compound number)
  # using the example above, out should be : 1753 0 0 0 0 0 0 0 0 0 0 247 0 0 0 0 0 0
  out[as.numeric(names(tmp))] <- tmp
  out
}