"ComputePosteriorRcpp_Add_Iso_Bio_Int_NoPot" <- function( P, Add, Iso, Bio, Int,
                                                          RT=NULL, relId=NULL,
                                                          corrMat=NULL,  RTwin=3,
                                                          corrThld=.80,it=1100, 
                                                          burn=100, delAdd=.5, delIso=.5,
                                                          delBio=1, ratioToll=0.8, log = F){
  if(is.null(RT) & is.null(relId) & is.null(corrMat)){
    cat("\n Missing RT, corrMat or relId")
    stop()
  }
  
  if(!is.null(RT)){
    counter <- 0
    remIdx <- lapply(RT, function(x){
      counter <<- counter + 1                                      # a flag
      CheckingRT(RT,counter, RTwin)                              # get a vector of which the elements are [counter, indexes] in RT which abs(RT[index] - RT[counter]) > RTwin
    })
  }else if(!is.null(relId)){
    counter <- 0
    remIdx <- lapply(relId, function(x){
      counter <<- counter + 1                                      # a flag
      CheckingRelId(relId,counter)                              # get a vector of which the elements are [counter, indexes] in relId which relId[index] != relId[counter]
    })
  }else{
    counter <- 0
    remIdx <- lapply(RT, function(x){
      counter <<- counter + 1                                      # a flag 
      CheckingCorr(corrMat[,counter],counter, corrThld)       # get a vector of which the elements are [indexes] in corrMat[,counter] which corrMat[,counter][index] < corr.thr
    })
  }

  post <- GibbsSampling_Int_NoPot(remIdx,
                                  Add,
                                  Iso,
                                  Bio,
                                  Int,
                                  P,
                                  delAdd,
                                  delIso,
                                  delBio,
                                  ratioToll,
                                  it,
                                  burn,
                                  log
  )
  return (post)
}