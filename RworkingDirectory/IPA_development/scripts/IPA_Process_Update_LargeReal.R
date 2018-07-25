rm(list = ls())
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/UpdateRelated/IPA_Update_LargeReal.cpp')
library(mzmatch.R)
mzmatch.init(version.1 = F)
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





# List UpdateMain(arma::sp_mat post_spMat,                // posterior matrix
#                 arma::sp_mat iso_spMat,                 // iso matrix
#                 NumericVector mass,                     // mass values
#                 NumericVector compMass,                 // compound mass values
#                 NumericVector massInt,                  // intensity values of masses
#                 NumericVector pk,                       // prior knowledge
#                 NumericVector massMZ,                   // mass charge ratio of masses
#                 NumericVector massRT,                   // mass retention time
#                 NumericVector obsMA,                    // storage of observed mass accuracies in this update time
#                 NumericVector recMA,                    // record all observed mass accuracies
#                 List recMAdistr,                        // record all the distributions of mass accuracies
#                 List obsIR,            // storage of observed intensity ratios of all detected isotope pairs in this update time
#                 NumericVector obsIRlComp,               // left compound number of stored observed intensity ratio
#                 NumericVector obsIRrComp,               // right compound number of stored observed intensity ratio
#                 List recIRdistr,                // record all distributions of all isotope pairs
#                 NumericVector recIRdistrLComp,          // left compound number of all recorded distributions of all isotope pairs
#                 NumericVector recIRdistrRComp,          // right compound number of all recorded distributions of all isotope pairs
#                 List obsRT,            // storage of observed retention time in this update time
#                 List recRT,            // record all observed retention time
#                 List recRTdistr,                // record all the distributions of retention time
#                 NumericVector obsRTcomp,                // the compound numbers of the observed retention time in this update time
#                 NumericVector recRTdistrComp,           // the compound numbers of recorded retention time distributions
#                 double ma,                              // mass accuracy
#                 double s_pk,                            // update scale(speed) for pk
#                 double s_ma,                            // update scale for ma
#                 double s_iso,                           // update scale for iso
#                 double threshold                        // threshold used to filter posterior
#                 
# )

post_spMat <- Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$Post
iso_spMat <- DB2.structured.POS$Iso.M
mass <- Masses[MassKept1]
compMass <- as.numeric(DB2.structured.POS$exact.masses.table[,4]) 
massInt <- Int
pk <- rep(1,length(compMass))
massMZ <- Masses[MassKept1]
massRT <- RTs
obsMA <- vector('numeric')
recMA <- vector('numeric')
recMAdistr <- NULL
obsIR <- NULL
obsIRlComp <- vector('numeric')
obsIRrComp <- vector('numeric')
recIRdistr <- NULL
recIRdistrLComp <- vector('numeric')
recIRdistrRComp <- vector('numeric')
obsRT <- NULL
recRT <- NULL
recRTdistr <- NULL
obsRTcomp <- vector('numeric')
recRTdistrComp <- vector('numeric')
ma <- 5
s_pk <- 0.1
s_ma <- 0.1
s_iso <- 0.1
threshold <- 0.9
rm(PeakML.Data)

l <- UpdateMainLargeReal(post_spMat, iso_spMat, mass, compMass, massInt, pk, massMZ, massRT, obsMA, recMA, recMAdistr,
                obsIR, obsIRlComp, obsIRrComp, recIRdistr, recIRdistrLComp, recIRdistrRComp, 
                obsRT, recRT, recRTdistr, obsRTcomp, recRTdistrComp, ma, s_pk, s_ma, s_iso, threshold)






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


