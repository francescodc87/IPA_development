rm(list = ls())
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/UpdateRelated/IPA_Update_LargeReal.cpp')
library(mzmatch.R)
mzmatch.init(version.1 = F)
# PeakML.Data <- PeakML.Read("/home/yuqiouyang/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data <- PeakML.Read("C:/Users/oyyqwhuiss/Documents/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
load("~/RworkingDirectory/IPA_development/data/DBs_all_compounds_only1entry.Rdata")
load("~/RworkingDirectory/IPA_development/data/Post_real_oneEntry_ADD_ISO_BIO_Int_noPen2.Rdata")


PeakML.Data.Table <- PeakML.Methods.getCompleteTable(PeakML.Data)
Masses <- apply(PeakML.Data.Table$Masses, 2, mean, na.rm=TRUE)                 # mass values: Numeric Vector (1xM) (1x7842)
RTs <- apply(PeakML.Data.Table$Retentiontimes,2, mean, na.rm=TRUE)             # RT values: Numeric Vector (1xM) (1x7842)
Intensities <- t(PeakML.Data.Table$Intensities)                                # Large matrix (MxG) (7842x38)
colnames(Intensities) <- PeakML.Data$sampleNames                               # set column names for Intensities matrix
Ids <- PeakML.Data$GroupAnnotations$id                                         # Character Vector (1xM) (1x7842)
RelId <- PeakML.Data$GroupAnnotations$relation.id                              # Character Vector (1xM) (1x7842)
PosAllData <- cbind(Masses, RTs, Intensities, Ids,RelId)                       # column bind those 5 datasets

tmp <- PeakML.Data$phenoData                                                   # Character Vector (1xG) (1x38)
tmp <- paste(tmp, c(1:7, 1:7,1:7,1:7,1:7,1:3), sep = "")                       # change element tmp vector: e.g "groupA1" to "groupA11"
colnames(PosAllData) <- c("m/z", "RTs", tmp, "Ids", "RelId")                   # Attach column names

# siege the data in Pos.All.data, based on multi-factors
MassKept1<- which(Masses>=90 & Masses<=476.18 & apply(Intensities,1, max, na.rm=T) > 8.59e+05 & RTs>=38.5 & RTs<=440.31)   # 1x3653
PosAllData <- PosAllData[MassKept1,]
# Intensities of masses: e.g, for one mass, the intensity is the maximum value of the mass being measured in G groups
Int <- as.numeric(apply(PosAllData[,3:40],1, max, na.rm=T))
RTs <- RTs[MassKept1]
compMass = DB2.structured.POS$exact.masses.table[,4]

post_spMat <- Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$Post
iso_spMat <- DB2.structured.POS$Iso.M
add_spMat <- DB2.structured.POS$Add.M
bio_spMat <- DB2.structured.POS$Bio.M

ma <- 3
massMZ <- Masses[MassKept1]
fomuMZ <- as.numeric(DB2.structured.POS$exact.masses.table[,4])
fomuCompStringId <- DB2.structured.POS$exact.masses.table[,1]
uniCompStringId <- unique(fomuCompStringId)
fomuCompId <- vector('numeric',length(fomuMZ))
for(i in 1:length(uniCompStringId)){
  stringId <- uniCompStringId[i]
  idx <- which(fomuCompStringId == stringId)
  fomuCompId[idx] <- i
}
# fomuCompId <- 
fomuCompId <- fomuCompId - 1
pkComp <- rep(1,length(uniCompStringId)) / length(uniCompStringId)

  

pkFomu <- rep(1,length(fomuMZ)) / length(fomuMZ)
obsMA <- vector('numeric')
recMA <- vector('numeric')
recMAdistr <- NULL
massInt <- Int
obsIR <- NULL
obsIRlFomu <- vector('numeric')
obsIRrFomu <- vector('numeric')
recIRdistr <- NULL
recIRdistrLFomu <- vector('numeric')
recIRdistrRFomu <- vector('numeric')

fomuId <- seq(length(fomuMZ)) - 1
boolMonoFomu <- as.numeric(DB2.structured.POS$exact.masses.table[,7])
monoFomuIdx <- which(boolMonoFomu == 0)
recMonoFomu <- fomuId[monoFomuIdx]

boolFomuMainAdd <- as.numeric(DB2.structured.POS$exact.masses.table[,6])
recCompMainAddFomu <- which(boolFomuMainAdd == 1) - 1

compId <- seq(length(uniCompStringId)) - 1
recMonoMainAddFomu <- vector('numeric',length(recMonoFomu))
for(i in 1:length(compId)){
  idx <- which(fomuCompId == (i - 1))
  recMonoMainAddFomu[idx] <- recCompMainAddFomu[i]
}

recCompMonoFomu <- list() # fomuCompId, recMonoFomu
monoFomuCompId <- vector('numeric',length(recMonoFomu))
for (i in 1:length(recMonoFomu)){
  monoFomuCompId[i] <- fomuCompId[recMonoFomu[i] + 1]
}

for (i in 1:length(uniCompStringId)){
  monoFomuCompIdx <- which(monoFomuCompId == (i - 1))
  subMonoFomuComp <- recMonoFomu[monoFomuCompIdx]
  recCompMonoFomu[[i]] <- subMonoFomuComp
}

recMonoFomuWithInSameComp <- list()
recMITwithInSameComp <- list()
for (i in 1:length(recMonoFomu)){
  monoCompId <-monoFomuCompId[i]
  subVec <- recCompMonoFomu[[monoCompId + 1]]
  subMITvec <- rep(0,length(subVec))
  mainAddFomu <- recCompMainAddFomu[monoCompId + 1]
  subVecMainAddIdx <- which(subVec == mainAddFomu)
  subMITvec[subVecMainAddIdx] <- subMITvec[subVecMainAddIdx] + 1
  recMonoFomuWithInSameComp[[i]] <- subVec
  recMITwithInSameComp[[i]] <- subMITvec
}


recCompBioLink <- list()
for (i in 1:length(uniCompStringId)){
  bioLinkVecSample <- sample(1:length(uniCompStringId), 2)
  bioLinkVecSample <- bioLinkVecSample - 1
  recCompBioLink[[i]] <- bioLinkVecSample
}
  
massRT <- RTs
obsFomuRT <- NULL
obsRTfomu <- vector('numeric')
recFomuRT <- NULL
recRTfomu <- vector('numeric')
recFomuRTdistr <- NULL
recRTdistrFomu <- vector('numeric')

s_pkAdd <- 0.1
s_pkComp <- 0.1
s_ma <- 0.1
s_iso <- 0.1
t_RT <- 5
t_post <- 0.9

# rm(PeakML.Data)

post_spMat <- post_spMat[1:length(massMZ), 1:length(fomuMZ)]

l <- UpdateMain_LargeReal(post_spMat, iso_spMat, add_spMat, bio_spMat, fomuCompId, pkFomu, pkComp, 
                          ma, massMZ, fomuMZ, obsMA, recMA, recMAdistr, 
                          massInt, obsIR, obsIRlFomu, obsIRrFomu, recIRdistr, recIRdistrLFomu, recIRdistrRFomu,
                          recMonoFomu, recMonoMainAddFomu, recMonoFomuWithInSameComp, recMITwithInSameComp, recCompMainAddFomu, recCompMonoFomu, recCompBioLink,
                          massRT, obsFomuRT, obsRTfomu, recFomuRT, recRTfomu, recFomuRTdistr, recRTdistrFomu,
                          s_pkAdd, s_pkComp, s_ma, s_iso, t_RT, t_post)



# List UpdateMain(arma::sp_mat post_spMat,                // posterior matrix
#                 arma::sp_mat iso_spMat,                 // iso matrix
#                 arma::sp_mat add_spMat,                 // add matrix
#                 arma::sp_mat bio_spMat,                 // bio matrix
#                 NumericVector fomuCompId,               // compound id of fomulas
#                 NumericVector pkFomu,                   // prior knowledge of all fomulas
#                 NumericVector pkComp,                   // prior knowledge of all compounds
#                 double ma,                              // mass accuracy
#                 NumericVector massMZ,                   // mass-charge-ratio values of masses
#                 NumericVector fomuMZ,                   // mass-charge-ratio values of chemical fomulas
#                 NumericVector obsMA,                    // storage of observed mass accuracies in this update time
#                 NumericVector recMA,                    // record of all observed mass accuracies
#                 List recMAdistr,                        // record of all the distributions of mass accuracies
#                 NumericVector massInt,                  // intensity values of masses
#                 List obsIR,                             // storage of observed intensity ratios of all detected isotope pairs in this update time
#                 NumericVector obsIRlFomu,               // storage of all the left fomula numbers of all the stored observed intensity ratios in this update time
#                 NumericVector obsIRrFomu,               // storage of all the right fomula numbers of all the stored observed intensity ratios in this update time
#                 List recIRdistr,                        // record of all the distributions of all recorded intensity ratios
#                 NumericVector recIRdistrLFomu,          // record of all the left fomula numbers of all recorded distributions
#                 NumericVector recIRdistrRFomu,          // record of all the right fomula numbers of all recorded distributions
#                 NumericVector recMonoFomu,              // record of all the fomula numbers of all monoisotopic adducts
#                 NumericVector recMonoMainAddFomu,       // record of all the fomula numbers of all monoisotopic adduct's main adducts
#                 List recMonoFomuWithInSameComp,         // the fomula numbers of all monoisotopic adducts of each monoisotopic compound (share the same compound)
#                 List recMITwithInSameComp,              // the most intense time values of monoisotopic adducts of each monoisotopic compound (share the same compound)
#                 NumericVector recCompMainAddFomu,       // fomula numbers of the main adducts of compounds
#                 List recCompMonoFomu,                   // the monoisotopic adduct fomula numbers of compounds
#                 List recCompBioLink,                    // the compound numbers of compounds (linked with biochemical reactions)
#                 NumericVector massRT,                   // mass retention time values of all the masses
#                 List obsFomuRT,                         // storage of observed retention time values of fomulas in this update time
#                 NumericVector obsRTfomu,                // the fomula numbers of the observed retention time in this update time
#                 List recFomuRT,                         // record all observed retention time values of fomulas
#                 NumericVector recRTfomu,                // the fomula numbers of the recorded retention time values
#                 List recFomuRTdistr,                    // record all the distributions of retention time values of fomulas
#                 NumericVector recRTdistrFomu,           // the fomula numbers of the recorded retention time distributions
#                 double s_pkAdd,                         // update scale(speed) for pkFomu
#                 double s_pkComp,                        // update scale(speed) for pkComp
#                 double s_ma,                            // update scale for ma (0< <=1)
#                 double s_iso,                           // update scale for iso (0< <=1)
#                 double t_RT,                            // threshold used to count iso and add connection when solving overlapped assignments automatically
#                 double t_post                           // threshold used to filter posterior
# 
# )



# data related to compounds:
# iso matrix *
# compound mass values 
# prior knowledge


# posterior matrix (3653 X 25098)


# other data:
# mass accuracy *
# s_pk, s_ma, s_iso *
# threshold *

# data related to masses:
# mass values *
# mass intensity values *
# mass charge ratios (or charge values) *
# mass retention time *

