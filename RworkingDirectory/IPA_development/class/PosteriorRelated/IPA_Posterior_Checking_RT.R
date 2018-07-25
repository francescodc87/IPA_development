# selecting things not to remove in each iteration according to RT
"CheckingRT" <- function(RT,i, RTwin){
  out <- c(i,which(abs(RT-RT[i]) > RTwin))
  out
}