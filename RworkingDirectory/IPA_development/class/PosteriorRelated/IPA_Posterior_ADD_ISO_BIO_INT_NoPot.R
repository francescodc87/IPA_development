"ComputePosteriorR_Add_Iso_Bio_Int_NoPot" <- function( P, Add, Iso, Bio, Int,
                                                       RT=NULL, relId=NULL,
                                                       corrMat=NULL,  RTwin=3,
                                                       corrThld=.80,it=1100, 
                                                       burn=100, delAdd=.5, delIso=.5,
                                                       delBio=1, allSamp=F,
                                                       unknownPen=NULL, ratioToll=0.8, log = F){
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
  
  
  Nc <- nrow(Add)
  M <- nrow(P)
  
  sampcomp <- apply(P,1,multsample)
  pot.bio <-apply(Bio[sampcomp,],2,sum)
  
  allsampcomp <- matrix(0,it, M) # matrix of all the samples
  v <- 0
  for (i in 1:it){
    w <- 0
    # cat("1")
    ordine <- sample(M)  # randomising the order used to cheack all assignments
    j <- 0
    for (thism in ordine){
      # cat("2")
      #counting adducts
      p.add <- colSums(matrix(Add[sampcomp[-remIdx[[thism]]],], ncol=Nc))
      # cat("3")
      ###counting isotopes
      tmp <- matrix(Iso[sampcomp[-remIdx[[thism]]],],ncol = Nc)*(Int[thism]/Int[-remIdx[[thism]]])
      ind.ones <- which((tmp>=ratioToll) & (tmp <=(1/ratioToll))) 
      # cat("ind.ones")
      # print(ind.ones)
      tmp[ind.ones]<-1
      tmp[tmp!=1] <-0
      p.iso<-colSums(tmp)
      # cat("4")
      ##counting biotransformations
      p.bio <- pot.bio - colSums(matrix(Bio[sampcomp[thism],], ncol=Nc))
      # cat("5")
      ### adding penalities
      if(!is.null(unknownPen)){
        p.add[length(p.add)] <- unknownPen
        p.iso[length(p.iso)] <- unknownPen
        
      } 
      
      ## normalising with deltas
      p.add <- (p.add + delAdd)/sum(p.add+delAdd)
      p.iso <- (p.iso + delIso)/sum(p.iso+delIso)
      p.bio <- (p.bio + delBio)/sum(p.bio+delBio)

      # cat("6")
      # cat("p.add")
      # print(p.add)
      # cat("p.iso")
      # print(p.iso)
      # cat("p.bio")
      # print(p.bio)
      ## merging scores
      po <- p.add*p.iso*p.bio*P[thism,]
      # cat("po")
      # print(po)
      po <- po/sum(po)
      oldval <- sampcomp[thism]
      sampcomp[thism] <- multsample(po)
      # cat("7")
      if(oldval!=sampcomp[thism]){
        pot.bio <- pot.bio - Bio[,oldval] 
        pot.bio <- pot.bio + Bio[,sampcomp[thism]]
      }
      
      if(log){
        old_w <- w
        w = (j * 100) %/% M
        if(old_w != w){
          cat("thism iteration: ", paste0(round(w, 0), "%", "\n"))
          if(w %% 10 == 0){
            cat("\014")
          }
        }
      } 
      j <- j + 1
    }
    allsampcomp[i,] <- sampcomp
    # cat("8")
    if(log){
      old_v <- v;
      v = (i * 100) / it
      if(old_v !=v){
        cat("computing posterior probability, ", paste0(round(v, 0), "%", "\n"))
      }
    } 
  }
  
  
  post <- t(apply(allsampcomp,2,ComputePost, burn=burn, it=it, Nc=Nc))
  if(allSamp){
    out <- list(Post=post, allsampcomp=allsampcomp)
    return(out)
  }else{
    return(post)  
  }
}