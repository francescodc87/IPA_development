rm(list = ls())
library(mzmatch.R)
library(RNeo4j)
library(slam)
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_NoRT_Rcpp.cpp')

source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_RT.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_relID.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_corr.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_INT_NoPot_Rcpp.R')
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_GS_INT_NoPot.cpp')   #'GS' = 'GibbsSampling'



## get data from masses
mzmatch.init(version.1 = F)
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_NoRT_Rcpp.cpp')
# PeakML.Data <- PeakML.Read("/home/yuqiouyang/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data <- PeakML.Read("C:/Users/oyyqwhuiss/Documents/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data.Table <- PeakML.Methods.getCompleteTable(PeakML.Data)
Mass <- apply(PeakML.Data.Table$Masses, 2, mean, na.rm=TRUE)                 # mass values: Numeric Vector (1xM) (1x7842)
RT <- apply(PeakML.Data.Table$Retentiontimes,2, mean, na.rm=TRUE)             # RT values: Numeric Vector (1xM) (1x7842)
Int <- t(PeakML.Data.Table$Intensities) 
Int <- as.numeric(apply(Int,1, max, na.rm=T))         # intensity values: Numeric Vector (1xM) (1x7842)
rm(PeakML.Data)

## set initial ppm, use this wo seach monoisotopic adducts in database
# every mass value should have a formula mass range based on ppm
initPPM <- 5
fomuMassMax <- vector('numeric')
fomuMassMin <- vector('numeric')
massNum <- length(Mass)
for(i in 1:massNum){
  mass <- Mass[i]
  fomuMassMax[i] <- mass / (1 - (initPPM / 1e6))
  fomuMassMin[i] <- mass / (1 + (initPPM / 1e6))
}



query = "
MATCH(n:Adduct)
WHERE n.mz >= {mz_min} AND n.mz <= {mz_max} 
RETURN n.id as id order by n.id
"
l <- cypherToList(graph, query, mz_min = 100, mz_max = 110)



# use that range to search mono adducts in database, return id
# you must have an index on property 'mz' of node 'Adduct' here, if not please create yourself
graph = startGraph("http://localhost:11004/db/data/",username="neo4j", password="oyyq6997")
query = "
MATCH(n:Adduct)-[:is_monoAdd]->(n:Adduct) 
WHERE n.mz >= {mz_min} AND n.mz <= {mz_max} 
RETURN n.id as id order by n.id
"
monoAdds <- vector('integer')
for (i in 1:massNum) {
  mz_min = fomuMassMin[i]
  mz_max = fomuMassMax[i]
  monoAdd <- cypherToList(graph, query, mz_min = mz_min, mz_max = mz_max)
  monoAdds <- c(monoAdds, unlist(monoAdd)) 
}
monoAdds <- unique(monoAdds)

# next seach for the 1st isotopes
query = "
MATCH(m:Adduct)-[:has_iso{isoRank:1}]->(n:Adduct) 
WHERE n.mz >= {mz_min} AND n.mz <= {mz_max}
RETURN m.id as addId, n.id as isoId order by m.id, n.id
"
firstIsos <- vector('integer')
firstIsoMonoAdds <- vector('integer')
for (i in 1:massNum) {
  mz_min = fomuMassMin[i]
  mz_max = fomuMassMax[i]
  firstIso <- cypherToList(graph, query, mz_min = mz_min, mz_max = mz_max)
  if(length(firstIso) != 0){
    for(j in 1:length(firstIso)){
      firstIsoMonoAdds <- c(firstIsoMonoAdds, firstIso[[j]]$addId)
      firstIsos <- c(firstIsos, firstIso[[j]]$isoId)
    }
  }
}
# the mono add which 1st isotope is related must be seen before
intersecMonoAdd <- intersect(monoAdds, unique(firstIsoMonoAdds)) 
# filter the 1st isotope
validFirstIso <- vector('integer')
for (i in 1:length(firstIsoMonoAdds)){
  if(length(which(intersecMonoAdd == firstIsoMonoAdds[i])) != 0){
    validFirstIso <- c(validFirstIso, firstIsos[i])
  }
}
# append adduct id
finalAdds <- c(monoAdds, unique(validFirstIso))




# next seach for the 2nd isotopes
query = "
MATCH (a:Adduct)<-[:has_iso{isoRank:1}]-(m:Adduct)-[:has_iso{isoRank:2}]->(b:Adduct) 
WHERE b.mz >= {mz_min} AND b.mz <= {mz_max}
RETURN a.id as fAddId, b.id as sAddId order by a.id, b.id
"
firstIsos <- vector('integer')
secondIsos <- vector('integer')
for (i in 1:massNum) {
  mz_min = fomuMassMin[i]
  mz_max = fomuMassMax[i]
  secondIso <- cypherToList(graph, query, mz_min = mz_min, mz_max = mz_max)
  if(length(secondIso) != 0){
    for(j in 1:length(secondIso)){
      firstIsos <- c(firstIsos, secondIso[[j]]$fAddId)
      secondIsos <- c(secondIsos, secondIso[[j]]$sAddId)
    }
  }
}
# the 1st isotope which 2nd isotope is related must be seen before
intersecFirstIso <- intersect(validFirstIso, unique(firstIsos)) 
# filter the 2nd isotope
validSecondIso <- vector('integer')
for (i in 1:length(intersecFirstIso)){
  if(length(which(firstIsos == intersecFirstIso[i])) != 0){
    validSecondIso <- c(validSecondIso, secondIsos[i])
  }
}
# append adduct id
finalAdds <- c(finalAdds, unique(validSecondIso))

# next seach for the 3rd isotopes
query = "
MATCH (b:Adduct)<-[:has_iso{isoRank:2}]-(m:Adduct)-[:has_iso{isoRank:3}]->(c:Adduct)
WHERE c.mz >= {mz_min} AND c.mz <= {mz_max}
RETURN b.id as sAddId, c.id as tAddId order by b.id, c.id
"
secondIsos <- vector('integer')
thirdIsos <- vector('integer')
for (i in 1:massNum) {
  mz_min = fomuMassMin[i]
  mz_max = fomuMassMax[i]
  thirdIso <- cypherToList(graph, query, mz_min = mz_min, mz_max = mz_max)
  if(length(thirdIso) != 0){
    for(j in 1:length(thirdIso)){
      secondIsos <- c(secondIsos, thirdIso[[j]]$sAddId)
      thirdIsos <- c(thirdIsos, thirdIso[[j]]$tAddId)
    }
  }
}
# the 2nd isotope which 3rd isotope is related must be seen before
intersecSecondIso <- intersect(validSecondIso, unique(secondIsos)) 
# filter the 3rd isotope
validThirdIso <- vector('integer')
for (i in 1:length(intersecSecondIso)){
  if(length(which(secondIsos == intersecSecondIso[i])) != 0){
    validThirdIso <- c(validThirdIso, thirdIsos[i])
  }
}


rm(firstIsoMonoAdds)
rm(firstIsos)
rm(secondIsos)
rm(thirdIsos)
rm(validFirstIso)
rm(validSecondIso)
rm(validThirdIso)
rm(intersecMonoAdd)
rm(intersecFirstIso)
rm(intersecSecondIso)
rm(monoAdds)


# append adduct id
finalAdds <- c(finalAdds, unique(validThirdIso))

###### we have all the formulas now, next is to construct the data need for IPA
# reoder id
finalAdds <- sort(finalAdds)


## prior computing
# masses and prior knowledge values of formulas
query = "
MATCH (n:Adduct)
WHERE n.id in {addId}
RETURN n.mz as fomuMass, n.pk as fomuPrior order by n.id
"
fomuMasses <- vector('numeric')
fomuPriors <- vector('numeric')
fomuMassInfo <- cypherToList(graph, query, addId = finalAdds)
for (i in 1:length(fomuMassInfo)) {
  fomuMasses <- c(fomuMasses, fomuMassInfo[[i]]$fomuMass)
  fomuPriors <- c(fomuPriors, fomuMassInfo[[i]]$fomuPrior)
}
rm(fomuMassInfo)
ma <- 5
unknowP <- 0.05
limit <- 1e-02
# save(Mass,fomuMasses,fomuPriors,finalAdds, file="C:/Users/oyyqwhuiss/Desktop/prior_data.Rdata")
# load(file = "C:/Users/oyyqwhuiss/Desktop/prior_data.Rdata")
priorRcpp <- ComputePriorRcpp(Mass,                       # mass values of masses
                              fomuMasses,                # mass values of formulas
                              fomuPriors,                # prior knowledge of formulas
                              NA,                        # retention time of masses (not applied yet)
                              ma,                        # mass accuracy
                              NULL,                      # retention time list (not applied yet)
                              unknowP,                   # value for unknown probability
                              limit,                     # limit used when probability value is too small
                              TRUE,                      # print log?
                              FALSE                      # apply retention time into prior computing?
)
# save(file="C:/Users/oyyqwhuiss/Desktop/prior_result06082018.Rdata")
# load(file="C:/Users/oyyqwhuiss/Desktop/prior_result06082018.Rdata")


## posterior computing
# priorRcpp <- priorRcpp[,1:length(finalAdds)]            # ignore unknow
massKept <- which(priorRcpp[,ncol(priorRcpp)] != 1)       # filter out those masses 100% assigned with unknow
fomuKept <- which(col_sums(priorRcpp)>0)                  # and filter out those formulas which do not have any mass to assigned with
finalAdds <- finalAdds[fomuKept]
Prior_filtered <- priorRcpp[massKept, fomuKept]  # ??? normalised?
RT <- RT[massKept]
Int <- Int[massKept]
fomuNum <- length(fomuKept)
# construct adduct matrix
query = "
MATCH (a:Adduct)<-[:is_monoAdd]-(a:Adduct)<-[:has_pos]-(m:Chemical)-[:has_pos]->(b:Adduct), (n:Chemical)-[:has_mainAdd]->(b:Adduct)
WHERE a.id in {addId} and b.id in {addId} and a.id <> b.id and m.id = n.id
RETURN distinct a.id as monoAddId, b.id as mainAddId order by a.id
"
posMonoAddIds <- vector('integer')
posMainAddIds <- vector('integer')
monoMainInfo <- cypherToList(graph, query, addId = finalAdds)
for (i in 1:length(monoMainInfo)) {
  posMonoAddIds <- c(posMonoAddIds, monoMainInfo[[i]]$monoAddId)
  posMainAddIds <- c(posMainAddIds, monoMainInfo[[i]]$mainAddId)
}
query = "
MATCH (a:Adduct)<-[:is_monoAdd]-(a:Adduct)<-[:has_neg]-(m:Chemical)-[:has_neg]->(b:Adduct), (n:Chemical)-[:has_mainAdd]->(b:Adduct)
WHERE a.id in {addId} and b.id in {addId} and a.id <> b.id and m.id = n.id
RETURN distinct a.id as monoAddId, b.id as mainAddId order by a.id
"
negMonoAddIds <- vector('integer')
negMainAddIds <- vector('integer')
monoMainInfo <- cypherToList(graph, query, addId = finalAdds)
for (i in 1:length(monoMainInfo)) {
  negMonoAddIds <- c(negMonoAddIds, monoMainInfo[[i]]$monoAddId)
  negMainAddIds <- c(negMainAddIds, monoMainInfo[[i]]$mainAddId)
}
monoAddIds <- c(posMonoAddIds, negMonoAddIds)
mainAddIds <- c(posMainAddIds, negMainAddIds)
Add <- Matrix(data = 0, nrow = fomuNum, ncol = fomuNum, sparse = TRUE, dimnames = list(as.character(finalAdds), as.character(finalAdds)))
for (i in 1:length(monoAddIds)){
  rowIdxString <- as.character(monoAddIds[i])
  colIdxString <- as.character(mainAddIds[i])
  Add[rowIdxString, colIdxString] <- 1
}
# construct iso matrix
query = "
MATCH (a:Adduct)-[r:has_iso]->(b:Adduct)
WHERE a.id in {addId} and b.id in {addId}
RETURN a.id as monoAddId, b.id as isoId, r.irInit as irMonoIso order by a.id, b.id
"
iso_monoAddIds <- vector('integer')
iso_isoIds <- vector('integer')
iso_irs <- vector('numeric')
isoPairInfo <- cypherToList(graph, query, addId = finalAdds)
for (i in 1:length(isoPairInfo)) {
  iso_monoAddIds <- c(iso_monoAddIds, isoPairInfo[[i]]$monoAddId)
  iso_isoIds <- c(iso_isoIds, isoPairInfo[[i]]$isoId)
  iso_irs <- c(iso_irs, isoPairInfo[[i]]$irMonoIso)
}
Iso <- Matrix(data = 0, nrow = fomuNum, ncol = fomuNum, sparse = TRUE, dimnames = list(as.character(finalAdds), as.character(finalAdds)))
for(i in 1:length(iso_monoAddIds)){
  rowIdxString <- as.character(iso_monoAddIds[i])
  colIdxString <- as.character(iso_isoIds[i])
  irValue <- iso_irs[i]
  Iso[rowIdxString, colIdxString] <- irValue
  Iso[colIdxString, rowIdxString] <- 1 / irValue
}
# construct bio matrix
# first, get valid compounds set
query = "
MATCH (m:Chemical)-[:has_mainAdd]->(n:Adduct)
WHERE n.id in {addId}
RETURN distinct m.id as compId, n.id as mainAddId, m.pk as compPk order by m.id
"
compIds <- vector('integer')
compMainAddIds <- vector('integer')
compPks <- vector('numeric')

compInfo <- cypherToList(graph, query, addId = finalAdds)
for (i in 1:length(compInfo)) {
  compIds <- c(compIds, compInfo[[i]]$compId)
  compMainAddIds <- c(compMainAddIds, compInfo[[i]]$mainAddId)
  compPks <- c(compPks, compInfo[[i]]$compPk)
}
# next, get valid compounds which are connected by reactions
query = "
MATCH p = (b:Chemical)<-[:has_reactant]-(a:Reaction)-[:has_reactant]->(c:Chemical) 
WHERE b.id in {validCompId} and c.id in {validCompId} and b.id <> c.id
RETURN distinct b.id as compId, c.id as linkerCompId order by b.id, c.id
"
bio_compIds <- vector('integer')
bio_linkerCompIds <- vector('integer')
bioCompInfo <- cypherToList(graph, query, validCompId = compIds)
for (i in 1:length(bioCompInfo)) {
  bio_compIds <- c(bio_compIds, bioCompInfo[[i]]$compId)
  bio_linkerCompIds <- c(bio_linkerCompIds, bioCompInfo[[i]]$linkerCompId)
}
# last, fill the bio matrix
Bio <- Matrix(data = 0, nrow = fomuNum, ncol = fomuNum, sparse = TRUE, dimnames = list(as.character(finalAdds), as.character(finalAdds)))
for(i in 1:length(unique(bio_compIds))){
  bio_comp <- bio_compIds[i]
  bio_linkComps <- bio_linkerCompIds[which(bio_compIds == bio_comp)]
  bio_compMainAdds <- compMainAddIds[which(compIds == bio_comp)]
  bio_linkCompMainAdds <- compMainAddIds[which(compIds %in% bio_linkComps)]
  for(j in 1:length(bio_compMainAdds)){
    bio_compMainAdd <- bio_compMainAdds[j]
    Bio[as.character(bio_compMainAdd),as.character(bio_linkCompMainAdds)] <- 1
    Bio[as.character(bio_linkCompMainAdds), as.character(bio_compMainAdd)] <- 1
  }
}
# save(Prior_filtered,Add,Iso,Bio,RT,Int, file="C:/Users/oyyqwhuiss/Desktop/posterior_data.Rdata")

## posterior computing:
RTwin <- 5
it <- 100
burn <- 10
delAdd <- 0.4
delIso <- 0.2
delBio <- 1
ratioToll<- 0.8
## caution: you may crash because of size restriction, go to {your R library}\RcppArmadillo\include\RcppArmadilloConfig.h to do the following change:
# assur this sentence is written in the file (# is included as command):
#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD
#endif
ComputePosteriorRcpp_Add_Iso_Bio_Int_NoPot(P = Prior_filtered, Add = Add,
                                          Iso = Iso, RT = RT,
                                          Bio = Bio, Int = Int,
                                          RTwin = RTwin, it = it, burn = burn,
                                          delAdd = delAdd, delIso = delIso, delBio = delBio,
                                          ratioToll = ratioToll, log = T)

## prepare for update stage
# information from the last experiment
query = "
MATCH (n:Experiment)-[:has_mainAdd]->(n:Adduct)
WHERE n.id in {monoId}
RETURN m.id as compId, m.pk as compPk order by m.id
"





# List UpdateMain(arma::sp_mat post_spMat,               ? // posterior matrix 
#                 arma::sp_mat iso_spMat,                ? // iso matrix
#                 arma::sp_mat add_spMat,                ? // add matrix
#                 arma::sp_mat bio_spMat,                ? // bio matrix
#                 NumericVector pkFomu,                  + // prior knowledge of all fomulas
#                 NumericVector pkComp,                  + // prior knowledge of all compounds
#                 double ma,                             + // mass accuracy
#                 NumericVector massMZ,                  - // mass-charge-ratio values of masses
#                 NumericVector fomuMZ,                  - // mass-charge-ratio values of chemical fomulas
#                 NumericVector obsMA,                   - // storage of observed mass accuracies in this update time
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
#                 double w_int,                           // weight for most intense mass when solving overlapped assignment
#                 double w_ma,                            // weight for min mass accuracy mass when solving overlapped assignment
#                 double w_addLink,                       // weight for adduct connection of mass when solving overlapped assignment
#                 double w_isoLink,                       // weight for isotope connection of of mass when solving overlapped assignment
#                 double t_RT,                            // threshold used to count iso and add connection when solving overlapped assignments automatically
#                 double t_post,                          // threshold used to filter posterior
#                 bool isTest = true                      // is it a test?
# 
# )







