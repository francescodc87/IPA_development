# removes all objects from the current workspace (R memory)
rm(list=ls())
library(rbenchmark)
library(mzmatch.R)
mzmatch.init(version.1 = F)

load("~/RworkingDirectory/IPA_development/data/DBs_all_compounds_only1entry.Rdata")  ##database
load("~/RworkingDirectory/IPA_development/data/POSPriors_oneEntry.Rdata")            ##Priors

source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_RT.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_relID.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_corr.R')

Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Sample_Rcpp.cpp')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Multsample_Rcpp.R')

source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Multsample_R.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Compute_post.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_INT_NoPot.R')

source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_INT_NoPot_Rcpp.R')
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_GS_INT_NoPot.cpp')   #'GS' = 'GibbsSampling'

# PeakML.Data <- PeakML.Read("/home/yuqiouyang/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data <- PeakML.Read("C:/Users/oyyqwhuiss/Documents/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")

## These are for debug use
# debug(IPA.sampler.Add.Iso.Bio)
# debug(IPA.sampler.Add.Iso.Bio.Rcpp)

PeakML.Data.Table <- PeakML.Methods.getCompleteTable(PeakML.Data)
dim(PeakML.Data.Table$Intensities) # retrieve dimension

Masses <- apply(PeakML.Data.Table$Masses, 2, mean, na.rm=TRUE)                 # mass values: Numeric Vector (1xM) (1x7842)
RTs <- apply(PeakML.Data.Table$Retentiontimes,2, mean, na.rm=TRUE)             # RT values: Numeric Vector (1xM) (1x7842)
Intensities <- t(PeakML.Data.Table$Intensities)                                # Large matrix (MxG) (7842x38)
colnames(Intensities) <- PeakML.Data$sampleNames                               # set column names for Intensities matrix
Ids <- PeakML.Data$GroupAnnotations$id                                         # Character Vector (1xM) (1x7842)
RelId <- PeakML.Data$GroupAnnotations$relation.id                             # Character Vector (1xM) (1x7842)
PosAllData <- cbind(Masses, RTs, Intensities, Ids,RelId)                    # column bind those 5 datasets

tmp <- PeakML.Data$phenoData                                                   # Character Vector (1xG) (1x38)
tmp <- paste(tmp, c(1:7, 1:7,1:7,1:7,1:7,1:3), sep = "")                       # change element tmp vector: e.g "groupA1" to "groupA11"
colnames(PosAllData) <- c("m/z", "RTs", tmp, "Ids", "RelId")               # Attach column names

# siege the data in Pos.All.data, based on multi-factors
MassKept1<- which(Masses>=90 & Masses<=476.18 & apply(Intensities,1, max, na.rm=T) > 8.59e+05 & RTs>=38.5 & RTs<=440.31)
PosAllData <- PosAllData[MassKept1,]

# rm(Masses, PeakML.Data, PeakML.Data.Table, RTs, Intensities, Ids, RelId, tmp)
Prior <- Prior[MassKept1, ]  # also siege data here, using the same filter
##save(Prior,file="~/RworkingDirectory/IPA_development/data/POSPriors_oneEntry.Rdata")


# adding unknowns in the connectivity matrices, for each connectivity matrix, using cbind and rbind to make both row and column +1
Add <- DB2.structured.POS$Add.M
Add <- cbind(Add, rep(0,nrow(Add)))
Add <- rbind(Add, rep(0,ncol(Add)))

Iso <- DB2.structured.POS$Iso.M
Iso <- cbind(Iso, rep(0,nrow(Iso)))
Iso <- rbind(Iso, rep(0,ncol(Iso)))

Bio<- DB2.structured.POS$Bio.M
Bio <- cbind(Bio, rep(0,nrow(Bio)))
Bio <- rbind(Bio, rep(0,ncol(Bio)))

# removing masses and formulas I don't need to consider
MassKept2 <- which(Prior[,ncol(Prior)]!=1)                              # filter out those masses 100% assigned with unknow
CompKept <- which(colSums(as.matrix(Prior))>0)                       # and filter out those compounds which do not have any mass to assigned with
Prior_filtered <- Prior[MassKept2, CompKept]
Bio <- Bio[CompKept, CompKept]
Iso <- Iso[CompKept, CompKept]
Add <- Add[CompKept, CompKept] 

# making binary Isotopes matrix
# Iso_bin <- Iso
# Iso_bin[Iso_bin>0] <- 1

RT <- as.numeric(PosAllData[,"RTs"])[MassKept2]                       # filter out those masses 100% assigned with unknow

# Intensities of masses: e.g, for one mass, the intensity is the maximum value of the mass being measured in G groups
Int <- as.numeric(apply(PosAllData[,3:40],1, max, na.rm=T)[MassKept2])
rm(PeakML.Data)


## Before running the script, we need to filter out some invalid or negative value (no need)
# Prior_filtered[Prior_filtered < 0 | is.na(Prior_filtered)] <- 0
# Add[Add < 0 | is.na(Add)] <- 0
# Iso_bin[Iso_bin < 0 | is.na(Iso_bin)] <- 0
# Bio[Bio < 0 | is.na(Bio)] <- 0





# ComputePosteriorRcpp_Add_Iso_Bio_Int_NoPot(P=Prior_filtered,Add = Add,
#                                            Iso = Iso, RT = RT,
#                                            Bio = Bio, Int = Int,
#                                            RTwin = 5, it =1100, burn = 100,
#                                            allSamp = F, delAdd =0.4,
#                                            delIso = 0.2, delBio =1, v = T)

# save(Prior_filtered,Add,Iso,RT,Bio,Int, file="~/RworkingDirectory/IPA_development/data/PosteriorEfficiencyTest_Rcpp_180718.Rdata")


system.time({ComputePosteriorRcpp_Add_Iso_Bio_Int_NoPot(P=Prior_filtered,Add = Add,
                                                        Iso = Iso, RT = RT,
                                                        Bio = Bio, Int = Int,
                                                        RTwin = 5, it =1100, burn = 100,
                                                        allSamp = F, delAdd =0.4,
                                                        delIso = 0.2, delBio =1, v = T)})

# system.time({ComputePosteriorR_Add_Iso_Bio_Int_NoPot(P=Prior_filtered,Add = Add,
#                                                         Iso = Iso, RT = RT,
#                                                         Bio = Bio, Int = Int,
#                                                         RTwin = 5, it =1100, burn = 100,
#                                                         allSamp = F, delAdd =0.4,
#                                                         delIso = 0.2, delBio =1, v = T)})


# test_row = nrow(Prior_filtered)
# test_col = ncol(Prior_filtered)
# for (oo in 1:test_row) {
#   for (pp in 1:test_col) {
#     if(!is.finite(Prior_filtered[oo,pp])){
#       cat("na detected: row: ",oo,", column: ",pp,"\n")
#     }
#   }
# }



## Actions after computing
# tmp.Post <- Prior
# tmp.Post[ind.masses.kept, ind.formulas.kept] <-Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$Post
# Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$Post <- tmp.Post
# rm(tmp.Post)
# 
# tmp.allsampcomp <- matrix(0,2000, nrow(Prior))
# tmp.allsampcomp[,ind.masses.kept] <- Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$allsampcomp
# Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$allsampcomp <- tmp.allsampcomp
# rm(tmp.allsampcomp)
# 
# #save(Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen, 
# #     file="~/data/Post_real_oneEntry_ADD_ISO_BIO_Int_noPen2.Rdata")