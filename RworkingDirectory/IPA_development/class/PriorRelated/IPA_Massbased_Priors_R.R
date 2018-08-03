# IPA code - R vesion - ONLY MASSES - COMPUTING ONLY THE PROBABILITIES ONLY BASED ON MASS
### mass:  the measured mass
### compMass: a vector containing the C theoretical masses 
### ppm: accuracy, this parameter is used to compute the precision parameter (gamma in the paper). It should be a multiple of the
### averege accuracy of the instrument


# list of optional parameters
### precision: mass noise precision (gamma in the paper), default 1  ### gotta change this with accuracy 
### noisetype: 0 for'absolute' or 1 for 'relative' noise (default 0)
### limit: the minimum probability accepted (default = 1e-06)

# out: list of outputs
### pr: a vectorcontaining C prior probabilities


"ComputePriorR" <- function(mass, compMass, ppm, limit =1e-02, unknownProb=0.05,
                                  pk=rep(1,length(compMass)),compId=NULL,
                                  RTs=NULL, RTranges=NA, RT_prior_penality=.5,
                                  v=FALSE, log=100){
  
  
  "SolveRTrange" <- function(RTrange,min_or_max){
    RTrange <- unlist(strsplit(RTrange, split=","))
      if(min_or_max=="min"){
        RTrange <- min(RTrange)
      }else{
        RTrange <- max(RTrange)
      }
      return(as.numeric(RTrange))
  }
  
  # loading packages..no need in the future
  library(Matrix)
  compMass <- as.numeric(compMass)
  
  # checking the parameters
  if(ppm <=0){
    cat('\n not allowed parameters values')
    stop()
  }
  if(length(compMass)!=length(pk)){
    cat('\n compMass and prior knowledge are not compatible')
    stop()
  }
  
  # defing the number of masses and the number of the compounds
  Nc <- length(compMass)
  M <- length(mass)
  
  # evaluating precision
  deltaMs <- ppm*mass*(1e-06)
  sigma <- deltaMs/2
  precision <- 1/(sigma^2)
  
  rm(sigma, deltaMs)
  
  # evaluation of prior probabilities (likelihood based only on mass), Nc + 1 is taking unknown into consideration
  # initialize some variables
  pr <-Matrix(0,M,(Nc+1))
  
  
  unknownProbs <- NULL
  
  # considering RTs and RT ranges for the 
  if(length(RTranges)==Nc){
    RTmin <- apply(matrix(RTranges,ncol=1),1,SolveRTrange,min_or_max="min")
    RTmax <- apply(matrix(RTranges,ncol=1),1,SolveRTrange,min_or_max="max")
  }
  
  
  
  for(i in 1:M){
    # computing RT probability
    RTprior <-rep(1,Nc)
    if(length(RTranges)==Nc){
    RTidx <- which(!is.na(RTranges))
    RTpenIdx <- RTidx[which(RTs[i]<RTmin[RTidx] | RTs[i]>RTmax[RTidx])]
    RTprior[RTpenIdx] <- RT_prior_penality
    }
    
    # main computing of prior probability
    pr[i,1:Nc] <- (exp((-0.5*precision[i])*((compMass-mass[i])^2)))*pk*RTprior
    # if probability value is too small, set it to 0
    if(!is.na(limit)){
      zerosIdx <- which(pr[i,]<(limit*sum(pr[i,])))
      if(length(zerosIdx)>0){
        pr[i,zerosIdx] <- 0  
      }
    }
    
    # computing "the unknown" probability
    unknownProbs  <- (unknownProb/(1-unknownProb))*sum(pr[i,])
    if(unknownProbs==0){unknownProbs <- 1}
    pr[i,Nc+1] <- unknownProbs
    pr[i,] <- pr[i,]/sum(pr[i,])
    
    
    if(v){
      if(i %% log==0) {
        # Print on the screen some message
        cat("Computing Prior in R, ",paste0(round((i*100)/M,0), "%", "\n"))
      }
    } 
    
  }
  
  if(!is.null(compId)){
    colnames(pr) <- c(compId,"unknown")
  }

    return(pr)
}
