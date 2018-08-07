# removes all objects from the current workspace (R memory)
# + various initialisations
rm(list = ls())
library(mzmatch.R)
mzmatch.init(version.1 = F)
load("~/RworkingDirectory/IPA_development/data/DBs_all_compounds_only1entry.Rdata")
source('~/RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_R.R')
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_NoRT_Rcpp.cpp')
# PeakML.Data <- PeakML.Read("/home/yuqiouyang/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data <- PeakML.Read("C:/Users/oyyqwhuiss/Documents/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")

#refine PeakML data into a multi element table
PeakML.Data.Table <- PeakML.Methods.getCompleteTable(PeakML.Data)
dim(PeakML.Data.Table$Intensities)

#prepare Pos.All.data
masses <- apply(PeakML.Data.Table$Masses, 2, mean, na.rm=TRUE)                 # mass values: Numeric Vector (1xM) (1x7842)
RTs <- apply(PeakML.Data.Table$Retentiontimes,2, mean, na.rm=TRUE)             # RT values: Numeric Vector (1xM) (1x7842)
intensities <- t(PeakML.Data.Table$Intensities)                                # Large matrix (MxG) (7842x38)
colnames(intensities) <- PeakML.Data$sampleNames                               # set column names for intensities matrix
ids <- PeakML.Data$GroupAnnotations$id                                         # Character Vector (1xM) (1x7842)
relId <- PeakML.Data$GroupAnnotations$relation.id                              # Character Vector (1xM) (1x7842)
posAllData <- cbind(masses, RTs, intensities, ids, relId)                      # column bind those 5 datasets

tmp <- PeakML.Data$phenoData                                                   # Character Vector (1xG) (1x38)
tmp <- paste(tmp, c(1:7, 1:7,1:7,1:7,1:7,1:3), sep = "")                       # change element tmp vector: e.g "groupA1" to "groupA11"
colnames(posAllData) <- c("m/z", "RTs", tmp, "ids", "rel.ids")                 # Attach column names

compMass = DB2.structured.POS$exact.masses.table[,4]

# # code to compute prior in R:
# # 1: Function definition
# # "ComputePriorR" <- function(mass, compMass, ppm, limit =1e-02, unknownProb=0.05,
# #                             pk=rep(1,length(compMass)),compId=NULL,
# #                             RTs=NULL, RTranges=NA, RT_prior_penality=.5,
# #                             v=FALSE, log=100)
# # 2: compute process, it took around 85mins
Prior_R <- ComputePriorR(mass = masses,
                         compMass = compMass, 
                         ppm = 3,
                         unknownProb=0.05, 
                         compId=1:length(compMass), 
                         v=T)
save(Prior_R, file="~/RworkingDirectory/IPA_development/data/Prior_R_oneEntry_180615.Rdata")

# # code to compute prior in Rcpp:
# # 1: Function definition
# # NumericMatrix ComputePriorRcpp(NumericVector mass,
# #                                NumericVector compMass,
# #                                NumericVector pk,
# #                                NumericVector RTs,
# #                                NumericVector ppm,
# #                                NumericVector compId,
# #                                List RTranges, // RTranges: A list, every element: RTmin, RTmax, RT_in, RT_out, is not applied here
# #                                double unknownProb,
# #                                double limit,
# #                                int log,
# #                                bool v = false,
# #                                bool RTrangeList = false)

# # 2: compute process
Prior_Rcpp <- NULL
Prior_Rcpp <- ComputePriorRcpp(masses,
                               as.numeric(compMass),
                               rep(1,length(compMass)),
                               NA,
                               3,
                               NULL,
                               0.05,
                               1e-02,
                               TRUE,
                               FALSE
                               )
save(Prior_Rcpp, file="~/RworkingDirectory/IPA_development/data/Prior_Rcpp_oneEntry_180806.Rdata")

## test agreement, can be achieved by load the data without running the code above
# rm(PeakML.Data)
# load(file="~/RworkingDirectory/IPA_development/data/Prior_R_oneEntry_180615.Rdata")
# load(file="~/RworkingDirectory/IPA_development/data/Prior_Rcpp_oneEntry_180806.Rdata")
# all.equal(Prior_R,Prior_Rcpp, check.attributes = FALSE)

## action after computing, prepare data for posterior computing
# Prior = Prior_Rcpp
# 
# rownames(Prior) <- ids
# 
# Prior.List <- list()
# for(i in 1:nrow(Prior)){
#   ind <- which(Prior[i,]>0)
#   if(length(ind)>0){
#     Prior.List[[i]] <- Prior[i,ind]
#     names(Prior.List[[i]]) <- colnames(Prior)[ind]
#   }else{
#     Prior.List[[i]] <- NA
#   }
# }
# 
# names(Prior.List) <- rownames(Prior)
# 
# save(Prior,Prior.List, file="~/RworkingDirectory/IPA_development/data/POSPriors_oneEntry.Rdata")





















## the following code is not used commonly

# # test agreement using self-construct data
# Rcpp::sourceCpp('scripts/test_Prior_Rcpp.cpp')
# source('~/RworkingDirectory/IPA_development/scripts/test_Prior_R.R')
# newMasses = c(81.5,91.8,106.99,117.76,122.3,146,156,167.79,177)
# newCompMass = c(77.8,87.8,97.8,107.8,117.8,127.8,137.8,147.8,157.8,167.8,177.8,187.8,197.8)

# # code to compute prior in R
# # "ComputePriorR" <- function(mass, compMass, ppm, limit =1e-02, unknownProb=0.05,
# #                             pk=rep(1,length(compMass)),compId=NULL,
# #                             RTs=NULL, RTranges=NA, RT_prior_penality=.5,
# #                             v=FALSE, log=100)
# 
# 
# Prior_R_test <- test_R_IPA.massbased.prior(newMasses,
#                                compounds.mass = newCompMass,ppm = 3000,
#                                unknown.prob=0.05, compounds.id=1:length(newCompMass),v=T)
# 
# # code to compute prior in Rcpp
# # NumericMatrix ComputePriorRcpp(NumericVector mass,
# #                                NumericVector compMass,
# #                                NumericVector pk,
# #                                NumericVector RTs,
# #                                NumericVector ppm,
# #                                NumericVector compId,
# #                                List RTranges, // RTranges: A list, every element: RTmin, RTmax, RT_in, RT_out, is not applied here
# #                                double unknownProb,
# #                                double limit,
# #                                int log,
# #                                bool v = false,
# #                                bool RTrangeList = false)
# 
# Prior_Rcpp_test <- test_computePriorRcpp(newMasses,
#                                     newCompMass,
#                                    rep(1,length(newCompMass)),
#                                    NA,
#                                    rep(3000,length(newCompMass)),
#                                    1:length(newCompMass),
#                                    NULL,
#                                    0.05,
#                                    1e-02,
#                                    100,
#                                    TRUE,
#                                    FALSE
# )












