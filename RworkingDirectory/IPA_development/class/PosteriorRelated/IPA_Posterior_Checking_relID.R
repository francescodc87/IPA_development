### checking rel.id
"CheckingRelId" <- function(relId,i){
  out <- c(i,which(relId!=relId[i]))
  out
}