rm(list = ls())
library(mzmatch.R)
library(RNeo4j)
library(slam)
library(Matrix)
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_NoRT_Rcpp.cpp')

source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_RT.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_relID.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_corr.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_INT_NoPot_Rcpp.R')
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_GS_INT_NoPot.cpp')   #'GS' = 'GibbsSampling'

Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/UpdateRelated/IPA_Update.cpp')
# source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Multsample_R.R')
# source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Compute_post.R')
# source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_R.R')
# source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_INT_NoPot.R')


## before running the scripts, please check the indexes of the database, database should have these indexes:
#  index on the 'id' property of 'Chemical' node, if not, please set 'id' as unique property of 'Chemical' node
#  index on the 'id' property of 'Adduct' node, if not, please set 'id' as unique property of 'Adduct' node
#  index on the 'mz' property of 'Adduct' node, if not, please set 'mz' as an index of 'Adduct' node

## user-defined parameter
ppm_query <- 5
ionMode <- "pos"
ma <- 5
unknowP <- 0.05
limit <- 1e-02
rtWin <- 5
it <- 2
burn <- 1
delAdd <- 0.4
delIso <- 0.2
delBio <- 1
ratioToll<- 0.8
s_pkAdd <- 0.1
s_pkComp <- 0.1
s_ma <- 0.1
s_iso <- 0.1
w_int <- 3
w_ma <- 1
w_addLink <- 1
w_isoLink <- 3
t_rt <- 3          # threshold used to count iso and add connection when solving overlapped assignments automatically
t_post <- 0.9        # threshold used to filter posterior
isTest <- TRUE
if(ionMode == "pos"){
  ionQueryString <- "[:has_pos]"
} else{
  ionQueryString <- "[:has_neg]"
}

## get data from masses
mzmatch.init(version.1 = F)
# PeakML.Data <- PeakML.Read("/home/yuqiouyang/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data <- PeakML.Read("C:/Users/oyyqwhuiss/Documents/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data.Table <- PeakML.Methods.getCompleteTable(PeakML.Data)
massMZ <- apply(PeakML.Data.Table$Masses, 2, mean, na.rm=TRUE)                 # mass values: Numeric Vector (1xM) (1x7842)
massRT <- as.numeric(apply(PeakML.Data.Table$Retentiontimes,2, mean, na.rm=TRUE))             # RT values: Numeric Vector (1xM) (1x7842)
massInt <- t(PeakML.Data.Table$Intensities) 
massInt <- as.numeric(apply(massInt,1, max, na.rm=T))         # intensity values: Numeric Vector (1xM) (1x7842)
massMZsort <- sort(massMZ)
rm(PeakML.Data)

dateStart <- Sys.time()

## CHEMICAL FORMULA QUERY
# set initial ppm, use this wo seach monoisotopic adducts in database
# every mass value should have a formula mass range based on ppm
fomuMassMax <- vector('numeric')
fomuMassMin <- vector('numeric')
massNum <- length(massMZ)
for(i in 1:massNum){
  mass <- massMZ[i]
  fomuMassMax[i] <- mass / (1 - (ppm_query / 1e6))
  fomuMassMin[i] <- mass / (1 + (ppm_query / 1e6))
}

# use that range to search mono adducts in database, return id
# you must have an index on property 'mz' of node 'Adduct' here, if not please create yourself
graph <- startGraph("http://localhost:11004/db/data/",username="neo4j", password="oyyq6997") # type your server and authority information
queryStringFront <- "MATCH(a:Adduct)-[:is_monoAdd]->(a:Adduct), (m:Chemical)-"
queryStringBack <- "->(b:Adduct) WHERE a.id = b.id and a.mz >= {mz_min} and a.mz <= {mz_max} RETURN distinct a.id as monoAddId order by a.id"
query <- paste(queryStringFront, ionQueryString, queryStringBack, sep = "")
monoAdds <- vector('integer')
monoAddIds <- vector('integer')
for (i in 1:massNum) {
  mz_min <- fomuMassMin[i]
  mz_max <- fomuMassMax[i]
  monoAdds <- cypherToList(graph, query, mz_min = mz_min, mz_max = mz_max)
  monoAddIds <- c(monoAddIds, unlist(monoAdds, use.names=FALSE)) 
}
rm(monoAdds)
monoAddIds <- unique(monoAddIds)

# next seach for the 1st isotopes
query <- "
MATCH (m:Adduct)-[:has_iso{isoRank:1}]->(n:Adduct) 
WHERE n.mz >= {mz_min} and n.mz <= {mz_max}  
WITH m, n
WHERE m.id in {monoAddId}
RETURN distinct n.id as isoId order by n.id
"
firstIsos <- vector('integer')
for (i in 1:massNum) {
  mz_min <- fomuMassMin[i]
  mz_max <- fomuMassMax[i]
  firstIsoInfo <- cypherToList(graph, query, mz_min = mz_min, mz_max = mz_max, monoAddId = monoAddIds)
  if(length(firstIsoInfo) != 0){
    for(j in 1:length(firstIsoInfo)){
      firstIsos <- c(firstIsos, firstIsoInfo[[j]]$isoId)
    }
  }
}
rm(firstIsoInfo)
firstIsos <- unique(firstIsos)



# next seach for the 2nd isotopes
query <- "
MATCH (a:Adduct)<-[:has_iso{isoRank:1}]-(m:Adduct)-[:has_iso{isoRank:2}]->(b:Adduct) 
WHERE b.mz >= {mz_min} and b.mz <= {mz_max}
WITH a, b 
WHERE a.id in {firstIsoId}
RETURN distinct b.id as isoId order by b.id
"
secondIsos <- vector('integer')
for (i in 1:massNum) {
  mz_min <- fomuMassMin[i]
  mz_max <- fomuMassMax[i]
  secondIsoInfo <- cypherToList(graph, query, mz_min = mz_min, mz_max = mz_max, firstIsoId = firstIsos)
  if(length(secondIsoInfo) != 0){
    for(j in 1:length(secondIsoInfo)){
      secondIsos <- c(secondIsos, secondIsoInfo[[j]]$isoId)
    }
  }
}
rm(secondIsoInfo)
secondIsos <- unique(secondIsos)

# next seach for the 3rd isotopes
query <- "
MATCH (a:Adduct)<-[:has_iso{isoRank:2}]-(m:Adduct)-[:has_iso{isoRank:3}]->(b:Adduct)
WHERE b.mz >= {mz_min} and b.mz <= {mz_max}
WITH a, b
WHERE a.id in {secondIsoId} 
RETURN b.id as isoId order by b.id
"
thirdIsos <- vector('integer')
for (i in 1:massNum) {
  mz_min <- fomuMassMin[i]
  mz_max <- fomuMassMax[i]
  thirdIsoInfo <- cypherToList(graph, query, mz_min = mz_min, mz_max = mz_max, secondIsoId = secondIsos)
  if(length(thirdIsoInfo) != 0){
    for(j in 1:length(thirdIsoInfo)){
      thirdIsos <- c(thirdIsos, thirdIsoInfo[[j]]$isoId)
    }
  }
}
rm(thirdIsoInfo)
finalAdds <- sort(unique(c(monoAddIds, firstIsos, secondIsos, thirdIsos)))

rm(monoAddIds)
rm(firstIsos)
rm(secondIsos)
rm(thirdIsos)


## PRIOR COMPUTING
# masses and prior knowledge values of formulas
query <- "
MATCH (n:Adduct)
WHERE n.id in {addId}
RETURN n.id as fomuId, n.mz as fomuMass, n.pk as fomuPrior order by n.id
"
fomuMZ <- vector('numeric')
fomuPk <- vector('numeric')
fomuInfo <- cypherToList(graph, query, addId = finalAdds)
for (i in 1:length(fomuInfo)) {
  fomuMZ <- c(fomuMZ, fomuInfo[[i]]$fomuMass)
  fomuPk <- c(fomuPk, fomuInfo[[i]]$fomuPrior)
}
rm(fomuInfo)

# get the mass accuracy info
query <- "
MATCH (n:Experiment)
RETURN max(n.id) as maxExprId
"
obsMA <- vector('numeric')
recMA <- vector('numeric')
recMAdistr <- NULL
maxExprId <- cypherToList(graph, query)
if (length(maxExprId[[1]]) != 0){
  expr_maxId <- maxExprId[[1]]$maxExprId
  query <- "
  MATCH (n:Experiment)
  WHERE n.id = {maxId}
  RETURN n.ppmUpdate as initMA, ppmAll as storedMA, ppmMean as maDistrMean, ppmStd as maDistrStd
  "
  exprInfo <- cypherToList(graph, query, maxId = expr_maxId)
  ma <- exprInfo[[1]]$initMA
  maVec <- exprInfo[[1]]$storedMA
  maMean <- exprInfo[[1]]$maDistrMean
  maStd <- exprInfo[[1]]$maDistrStd
  if (length(maVec) >= 2){
    recMA <- unlist(maVec, use.names=FALSE)
    recMAdistr <- list(ma_num = length(maVec), ma_mean = maMean, ma_std = maStd)
    rm(maVec)
  } else if (length(maVec) == 1){
    obsMA <- unlist(maVec, use.names=FALSE)
  } 
} 
rm(exprInfo)

# save(massMZ,fomuMZ,fomuPk,finalAdds, file="C:/Users/oyyqwhuiss/Desktop/prior_data12082018.Rdata")
# load(file = "C:/Users/oyyqwhuiss/Desktop/prior_data12082018.Rdata")
priorRcpp <- ComputePriorRcpp(massMZ,                    # mass values of masses
                              fomuMZ,                    # mass values of formulas
                              fomuPk,                    # prior knowledge of formulas
                              NA,                        # retention time of masses (not applied yet)
                              ma,                        # mass accuracy
                              NULL,                      # retention time list (not applied yet)
                              unknowP,                   # value for unknown probability
                              limit,                     # limit used when probability value is too small
                              TRUE,                      # print log?
                              FALSE                      # apply retention time into prior computing?
)
# save(priorRcpp, file="C:/Users/oyyqwhuiss/Desktop/prior_result12082018.Rdata")
# load(file="C:/Users/oyyqwhuiss/Desktop/prior_result12082018.Rdata")


## POSTERIOR COMPUTING
rowName <- as.character(seq(nrow(priorRcpp))) 
colName <- c(as.character(finalAdds), "unknown")
rownames(priorRcpp) <- rowName
colnames(priorRcpp) <- colName
massKept <- which(priorRcpp[,ncol(priorRcpp)] != 1)       # filter out those masses 100% assigned with unknow
fomuKept <- which(col_sums(priorRcpp)>0)                  # and filter out those formulas which do not have any mass to assigned with
prior_filtered <- priorRcpp[massKept, fomuKept]           # ??? normalised?
massRT <- massRT[massKept]
massInt <- massInt[massKept]
lastCol <- fomuKept[length(fomuKept)]
if(lastCol == ncol(priorRcpp)){                           # unknow has been kept
  finalAdds <- finalAdds[fomuKept[1:(length(fomuKept) - 1)]]
  newColName <- c(as.character(finalAdds), "unknown")
} else{
  finalAdds <- finalAdds[fomuKept]
  newColName <- c(as.character(finalAdds))
}
colNameIdx <- finalAdds
fomuNum <- length(fomuKept)

# get valid compound set+mono adduct set
queryStringFront <- "MATCH (m:Chemical)-"
queryStringBack <- "->(n:Adduct)-[:is_monoAdd]->(n:Adduct) WHERE n.id in {addId} RETURN distinct m.id as compId, n.id as monoAddId order by m.id, n.id"
query <- paste(queryStringFront, ionQueryString, queryStringBack, sep = "")
monoComp_valiCompIds <- vector('integer')
monoComp_valiMonoAddIds <- vector('integer')
compMonoInfo <- cypherToList(graph, query, addId = finalAdds)
for (i in 1:length(compMonoInfo)) {
  monoComp_valiCompIds <- c(monoComp_valiCompIds, compMonoInfo[[i]]$compId)
  monoComp_valiMonoAddIds <- c(monoComp_valiMonoAddIds, compMonoInfo[[i]]$monoAddId)
}
rm(compMonoInfo)
recMonoComp <- monoComp_valiCompIds
recMonoFomu <- monoComp_valiMonoAddIds
monoComp_valiCompIds <- sort(unique(monoComp_valiCompIds))
monoComp_valiMonoAddIds <- sort(unique(monoComp_valiMonoAddIds))

# get valid all adducts+compound set set
queryStringFront <- "MATCH (m:Chemical)-"
queryStringBack <- "->(n:Adduct) WHERE n.id in {addId} RETURN distinct n.id as addId, m.id as compId order by n.id, m.id"
query <- paste(queryStringFront, ionQueryString, queryStringBack, sep = "")
recAllAddFomu <- vector('integer')
recAllAddComp <- vector('integer')
addCompInfo <- cypherToList(graph, query, addId = finalAdds)
for (i in 1:length(addCompInfo)) {
  recAllAddFomu <- c(recAllAddFomu, addCompInfo[[i]]$addId)
  recAllAddComp <- c(recAllAddComp, addCompInfo[[i]]$compId)
}
all_valiCompIds <- sort(unique(recAllAddComp))
rm(addCompInfo)


# get compPK
query <- "
MATCH (n:Chemical)
WHERE n.id in {compId} 
RETURN n.id as compId, n.pk as compPK order by n.id
"
compPK <- vector('numeric')
compPKinfo <- cypherToList(graph, query, compId = recAllAddComp)
for (i in 1:length(compPKinfo)) {
  compPK <- c(compPK, compPKinfo[[i]]$compPK)
}
rm(compPKinfo)

# construct adduct matrix
queryStringFront <- "MATCH (a:Adduct)<-[:is_monoAdd]-(a:Adduct)<-"
queryStringMid <- "-(m:Chemical)-"
queryStringBack <- "->(b:Adduct), (m:Chemical)-[:has_mainAdd]->(b:Adduct) WHERE a.id in {monoAddId} and b.id in {addId} and a.id <> b.id and m.id in {compId} RETURN distinct a.id as monoAddId, b.id as mainAddId order by a.id"
query <- paste(queryStringFront, ionQueryString, queryStringMid, ionQueryString, queryStringBack, sep = "")
add_monoAddIds <- vector('integer')
add_mainAddIds <- vector('integer')
add_monoMainInfo <- cypherToList(graph, query, monoAddId = monoComp_valiMonoAddIds, addId = finalAdds, compId = monoComp_valiCompIds)
for (i in 1:length(add_monoMainInfo)) {
  add_monoAddIds <- c(add_monoAddIds, add_monoMainInfo[[i]]$monoAddId)
  add_mainAddIds <- c(add_mainAddIds, add_monoMainInfo[[i]]$mainAddId)
}
Add <- Matrix(data = 0, nrow = fomuNum, ncol = fomuNum, sparse = TRUE, dimnames = list(newColName, newColName))
for (i in 1:length(add_monoAddIds)){
  rowIdxString <- as.character(add_monoAddIds[i])
  colIdxString <- as.character(add_mainAddIds[i])
  Add[rowIdxString, colIdxString] <- 1
  Add[colIdxString, rowIdxString] <- 1
}
rm(add_monoMainInfo)

# construct iso matrix
query <- "
MATCH (a:Adduct)-[r:has_iso]->(b:Adduct)
WHERE a.id in {monoAddId} and b.id in {addId}
RETURN distinct a.id as monoAddId, b.id as isoId, r.irValue as irMonoIso order by a.id, b.id
"
iso_monoAddIds <- vector('integer')
iso_isoIds <- vector('integer')
iso_irs <- vector('numeric')
iso_pairInfo <- cypherToList(graph, query, monoAddId = monoComp_valiMonoAddIds, addId = {finalAdds})
for (i in 1:length(iso_pairInfo)) {
  iso_monoAddIds <- c(iso_monoAddIds, iso_pairInfo[[i]]$monoAddId)
  iso_isoIds <- c(iso_isoIds, iso_pairInfo[[i]]$isoId)
  iso_irs <- c(iso_irs, iso_pairInfo[[i]]$irMonoIso)
}
Iso <- Matrix(data = 0, nrow = fomuNum, ncol = fomuNum, sparse = TRUE, dimnames = list(newColName, newColName))
for(i in 1:length(iso_monoAddIds)){
  rowIdxString <- as.character(iso_monoAddIds[i])
  colIdxString <- as.character(iso_isoIds[i])
  irValue <- iso_irs[i]
  Iso[rowIdxString, colIdxString] <- irValue
  Iso[colIdxString, rowIdxString] <- 1 / irValue
}
rm(iso_pairInfo)
# construct bio matrix
# first, get valid compounds set
query <- "
MATCH (m:Chemical)-[:has_mainAdd]->(n:Adduct)
WHERE m.id in {compId} and n.id in {addId}
RETURN distinct m.id as compId, n.id as mainAddId order by m.id, n.id
"
bio_compIds <- vector('integer')
bio_compMainAddIds <- vector('integer')
bio_compMainAddInfo <- cypherToList(graph, query, compId = monoComp_valiCompIds, addId = finalAdds)
for (i in 1:length(bio_compMainAddInfo)) {
  bio_compIds <- c(bio_compIds, bio_compMainAddInfo[[i]]$compId)
  bio_compMainAddIds <- c(bio_compMainAddIds, bio_compMainAddInfo[[i]]$mainAddId)
}
rm(bio_compMainAddInfo)

# next, get valid compounds which are connected by reactions
query <- "
MATCH (a:Chemical)<-[:has_reactant]-(:Reaction)-[:has_reactant]->(b:Chemical)
WHERE a.id in {validCompId} and b.id in {validCompId} and a.id <> b.id
RETURN distinct a.id as reacCompId, b.id as reacLinkerCompId order by a.id, b.id
"
bio_reac_compIds <- vector('integer')
bio_reac_linkerCompIds <- vector('integer')
bio_compReacInfo <- cypherToList(graph, query, validCompId = bio_compIds)
for (i in 1:length(bio_compReacInfo)) {
  bio_reac_compIds <- c(bio_reac_compIds, bio_compReacInfo[[i]]$reacCompId)
  bio_reac_linkerCompIds <- c(bio_reac_linkerCompIds, bio_compReacInfo[[i]]$reacLinkerCompId)
}
rm(bio_compReacInfo)
query <- "
MATCH (a:Chemical)<-[:has_cofactor]-(:Reaction)-[:has_cofactor]->(b:Chemical)
WHERE a.id in {validCompId} and b.id in {validCompId} and a.id <> b.id
RETURN distinct a.id as cofCompId, b.id as cofLinkerCompId order by a.id, b.id
"
bio_cof_compIds <- vector('integer')
bio_cof_linkerCompIds <- vector('integer')
bio_compCofInfo <- cypherToList(graph, query, validCompId = bio_compIds)
for (i in 1:length(bio_compCofInfo)) {
  bio_cof_compIds <- c(bio_cof_compIds, bio_compCofInfo[[i]]$cofCompId)
  bio_cof_linkerCompIds <- c(bio_cof_linkerCompIds, bio_compCofInfo[[i]]$cofLinkerCompId)
}
rm(bio_compCofInfo)

# last, fill the bio matrix
Bio <- Matrix(data = 0, nrow = fomuNum, ncol = fomuNum, sparse = TRUE, dimnames = list(newColName, newColName))
for(i in 1:length(bio_reac_compIds)){
  bio_reac_comp <- bio_reac_compIds[i]
  bio_reac_linkcomp <- bio_reac_linkerCompIds[i]
  bio_reac_compMainAdd <- bio_compMainAddIds[which(bio_compIds == bio_reac_comp)]
  bio_reac_linkCompMainAdd <- bio_compMainAddIds[which(bio_compIds == bio_reac_linkcomp)]
  Bio[as.character(bio_reac_compMainAdd),as.character(bio_reac_linkCompMainAdd)] <- 1
}
for(i in 1:length(bio_cof_compIds)){
  bio_cof_comp <- bio_cof_compIds[i]
  bio_cof_linkcomp <- bio_cof_linkerCompIds[i]
  bio_cof_compMainAdd <- bio_compMainAddIds[which(bio_compIds == bio_cof_comp)]
  bio_cof_linkCompMainAdd <- bio_compMainAddIds[which(bio_compIds == bio_cof_linkcomp)]
  Bio[as.character(bio_cof_compMainAdd),as.character(bio_cof_linkCompMainAdd)] <- 1
}
# save(prior_filtered,Add,Iso,Bio,massRT,massInt,newColName, file="C:/Users/oyyqwhuiss/Desktop/posterior_data12082018.Rdata")
# load(file="C:/Users/oyyqwhuiss/Desktop/posterior_data12082018.Rdata")
 
## CAUTION: you may crash because of size restriction, go to {your R library}\RcppArmadillo\include\RcppArmadilloConfig.h to do the following change:
# assure the following command is written in that file (# is included as command):
#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD
#endif

# it spent around 5 minutes to run each iteration in "it" parameter
Posterior <- ComputePosteriorRcpp_Add_Iso_Bio_Int_NoPot(P = prior_filtered, Add = Add,
                                                        Iso = Iso, RT = massRT,
                                                        Bio = Bio, Int = massInt,
                                                        RTwin = rtWin, it = it, burn = burn,
                                                        delAdd = delAdd, delIso = delIso, delBio = delBio,
                                                        ratioToll = ratioToll, log = T)



Posterior <- prior_filtered



# deal with 'unknown' before update
if(lastCol == ncol(priorRcpp)){                           # unknow has been kept
  Posterior <- Posterior[,1:(ncol(Posterior) - 1)]
}

## UPDATE STAGE 1: ALGORITHM 
fomuPk <- fomuPk[fomuKept]
massMZ <- massMZ[massKept]
fomuMZ <- fomuMZ[fomuKept]

# get isotope info: obsIR, obsIRmonoFomu, obsIRisoFomu; 
#                   recIRdistr, recIRdistrMonoFomu, recIRdistrIsoFomu
query <- "
MATCH (m:Adduct)-[r:has_iso]->(n:Adduct)
WHERE m.id in {monoIds} and n.id in {isoIds} and r.irObs <> 0
RETURN m.id as monoId, n.id as isoId, r.irObs as irObs order by m.id, n.id
"
obsIRmonoFomu <- vector('integer')
obsIRisoFomu <- vector('integer')
obsIR <- NULL
obsIRInfo <- cypherToList(graph, query, monoIds = monoComp_valiMonoAddIds, isoIds = iso_isoIds)
if(length(obsIRInfo) != 0){
  obsIR <- list()
  for (i in 1:length(obsIRInfo)) {
    obsIRmonoFomu <- c(obsIRmonoFomu, obsIRInfo[[i]]$monoId)
    obsIRisoFomu <- c(obsIRisoFomu, obsIRInfo[[i]]$isoId)
    obsIR[[length(obsIR) + 1]] <- c(obsIRInfo[[i]]$irObs)
  }
} 
rm(obsIRInfo)

query <- "
MATCH (m:Adduct)-[r:has_iso]->(n:Adduct)
WHERE m.id in {monoIds} and n.id in {isoIds} and r.irNum <> 0
RETURN m.id as monoId, n.id as isoId, r.irNum as irNum, r.irMean as irMean, r.irLogMean as irLogMean, r.irStd as irStd, r.irLogStd as irLogStd order by m.id, n.id
"
recIRdistrMonoFomu <- vector('integer')
recIRdistrIsoFomu <- vector('integer')
recIRirNum <- vector('integer')
recIRirMean <- vector('numeric')
recIRirLogMean <- vector('numeric')
recIRirStd <- vector('numeric')
recIRirLogStd <- vector('numeric')
recIRdistr <- NULL
recIRInfo <- cypherToList(graph, query, monoIds = monoComp_valiMonoAddIds, isoIds = iso_isoIds)
if(length(recIRInfo) != 0){
  for (i in 1:length(recIRInfo)) {
    recIRdistrMonoFomu <- c(recIRdistrMonoFomu, recIRInfo[[i]]$monoId)
    recIRdistrIsoFomu <- c(recIRdistrIsoFomu, recIRInfo[[i]]$isoId)
    recIRirNum <- c(recIRirNum, recIRInfo[[i]]$irNum)
    recIRirMean <- c(recIRirMean, recIRInfo[[i]]$irMean)
    recIRirLogMean <- c(recIRirLogMean, recIRInfo[[i]]$irLogMean)
    recIRirStd <- c(recIRirStd, recIRInfo[[i]]$irStd)
    recIRirLogStd <- c(recIRirLogStd, recIRInfo[[i]]$irLogStd)
  }
  recIRdistr <- list()
  for (i in 1:length()){
    recIRdistr[[i]] <- list(iso_num = recIRirNum[i], iso_mean = recIRirMean[i], iso_logMean = recIRirLogMean[i], iso_std = recIRirStd[i], iso_logStd = recIRirLogStd[i])
  }
} 
rm(recIRInfo)

# recComp, recCompMainAddFomu 
recComp <- monoComp_valiCompIds
queryStringFront <- "MATCH (a:Adduct)<-"
queryStringBack <- "-(m:Chemical), (m:Chemical)-[:has_mainAdd]->(b:Adduct) WHERE m.id in {compId} and a.id = b.id RETURN distinct m.id as compId, a.id as mainAddId order by m.id, a.id"
query <- paste(queryStringFront, ionQueryString, queryStringBack, sep = "")
recComp_compIds <- vector('integer')
recCompMainAddFomu <- vector('integer')
recComp_compMainAddInfo <- cypherToList(graph, query, compId = recComp)
for (i in 1:length(recComp_compMainAddInfo)) {
  recComp_compIds <- c(recComp_compIds, recComp_compMainAddInfo[[i]]$compId)
  recCompMainAddFomu <- c(recCompMainAddFomu, recComp_compMainAddInfo[[i]]$mainAddId)
}
rm(recComp_compMainAddInfo)

# recCompMonoFomu, recMITwithInSameComp
queryStringFront <- "MATCH (m:Chemical)-"
queryStringBack <- "->(n:Adduct)-[r:is_monoAdd]->(n:Adduct) WHERE m.id in {compId} RETURN distinct m.id as compId, n.id as monoAddId, r.mostIntTime as mit order by m.id, n.id"
query <- paste(queryStringFront, ionQueryString, queryStringBack, sep = "")
recComp_dupCompIds <- vector('integer')
recComp_compMonoAddIds <- vector('integer')
recComp_compMonoMIT <- vector('integer')
recComp_compMonoAddInfo <- cypherToList(graph, query, compId = recComp)
for (i in 1:length(recComp_compMonoAddInfo)) {
  recComp_dupCompIds <- c(recComp_dupCompIds, recComp_compMonoAddInfo[[i]]$compId)
  recComp_compMonoAddIds <- c(recComp_compMonoAddIds, recComp_compMonoAddInfo[[i]]$monoAddId)
  recComp_compMonoMIT <- c(recComp_compMonoMIT, recComp_compMonoAddInfo[[i]]$mit)
}
recCompMonoFomu <- list()
recMITwithInSameComp <- list()
for (i in 1:length(recComp)){
  compId <- recComp[i]
  findCompIdx <- which(recComp_dupCompIds == compId)
  recCompMonoFomu[[i]] <- recComp_compMonoAddIds[findCompIdx]
  recMITwithInSameComp[[i]] <- recComp_compMonoMIT[findCompIdx]
}
rm(recComp_compMonoAddInfo)

# obsCompRT, obsRTcomp, recCompRT, recRTcomp, recCompRTdistr, recRTdistrComp
query <- "
MATCH (m:Chemical)
WHERE m.id in {compId}
RETURN m.id as compId, m.rtAll as compAllRT, m.rtMean as compRTmean, m.rtStd as compRTstd order by m.id
"
obsCompRT <- vector('numeric')
obsRTcomp <- vector('integer')
recCompRT <- list()
recRTcomp <- vector('integer')
recCompRTdistr <- list()
recRTdistrComp <- vector('integer')
compRTInfo <- cypherToList(graph, query, compId = all_valiCompIds)
for (i in 1:length(compRTInfo)) {
  rtVec <- compRTInfo[[i]]$compAllRT
  compId <- compRTInfo[[i]]$compId
  rtMean <- compRTInfo[[i]]$compRTmean
  rtStd <- compRTInfo[[i]]$compRTstd
  if (length(rtVec) >= 2){
    recCompRT[[length(recCompRT)+1]] <- unlist(rtVec, use.names=FALSE)
    recRTcomp <- c(recRTcomp, compId)
    recCompRTdistr[[length(recCompRTdistr)+1]] <- list(rt_num = length(rtVec), rt_mean = rtMean, rt_std = rtStd)
    recRTdistrComp <- c(recRTdistrComp, compId)
  } else if (length(rtVec) == 1){
    obsCompRT <- c(obsCompRT, unlist(rtVec, use.names=FALSE))
    obsRTcomp <- c(obsRTcomp, compId)
  } 
}
rm(compRTInfo)


# save(colNameIdx, Posterior, Iso, Add, fomuPk, compPK, ma, massMZ, fomuMZ, obsMA, recMA, recMAdistr,  
#      massInt, obsIR, obsIRmonoFomu, obsIRisoFomu, recIRdistr, recIRdistrMonoFomu, recIRdistrIsoFomu,
#      recMonoFomu, recAllAddFomu, recAllAddComp, recComp, recCompMainAddFomu, recCompMonoFomu, recMITwithInSameComp,
#      massRT, obsCompRT, obsRTcomp, recCompRT, recRTcomp, recCompRTdistr, recRTdistrComp,
#      s_pkAdd, s_pkComp, s_ma, s_iso, w_int, w_ma, w_addLink, w_isoLink, t_rt, t_post, file="C:/Users/oyyqwhuiss/Desktop/update_data12082018.Rdata")
# load(file="C:/Users/oyyqwhuiss/Desktop/update_data12082018.Rdata")

l <- list()
l <- UpdateMain(colNameIdx, Posterior, Iso, Add, fomuPk, compPK, ma, massMZ, fomuMZ, obsMA, recMA, recMAdistr,  
                massInt, obsIR, obsIRmonoFomu, obsIRisoFomu, recIRdistr, recIRdistrMonoFomu, recIRdistrIsoFomu,
                recMonoFomu, recAllAddFomu, recAllAddComp, recComp, recCompMainAddFomu, recCompMonoFomu, recMITwithInSameComp,
                massRT, obsCompRT, obsRTcomp, recCompRT, recRTcomp, recCompRTdistr, recRTdistrComp,
                s_pkAdd, s_pkComp, s_ma, s_iso, w_int, w_ma, w_addLink, w_isoLink, t_rt, t_post, TRUE)




## UPDATE STAGE 2: UPDATE THE DATABASE
# update prior knowledge: finalAdds, all_valiCompIds
fomuPk <- l$pk$pkFomu
compPK <- l$pk$pkComp
map <- list()
query <- "
UNWIND {map} as row
MATCH (n:Adduct)
WHERE n.id = row.addId 
SET n.pk = row.addPK
"
for (i in 1:length(fomuPk)){
  map[[i]] <- list(addId = finalAdds[i], addPK = fomuPk[i])
}
cypher(graph, query, map = map)

map <- list()
query <- "
UNWIND {map} as row
MATCH (n:Chemical)
WHERE n.id = row.compId
SET n.pk = row.compPK
"
for (i in 1:length(compPK)){
  map[[i]] <- list(compId = all_valiCompIds[i], compPK = compPK[i])
}
cypher(graph, query, map = map)

# update iso connectivities
retainIRidx <- l$iso$retainIRidx
removeIRidx <- l$iso$removeIRidx
obsIRmonoFomu <- l$iso$obsIRmonoFomu
obsIRisoFomu <- l$iso$obsIRisoFomu
obsIR <- l$iso$obsIR
query <- "
UNWIND {map} as row
MATCH (m:Adduct)-[r:has_iso]->(n:Adduct)
WHERE m.id = row.monoAddId and n.id = row.isoId
SET r.irObs = row.irObs
"
map <- list() 
if(length(retainIRidx) > 0){
  for (i in 1:length(retainIRidx)){
    idx <- retainIRidx[i]
    map[[i]] <- list(monoAddId = obsIRmonoFomu[idx], isoId = obsIRisoFomu[idx], irObs = obsIR[idx])
  }
  cypher(graph, query, map = map)
}

query <- "
UNWIND {map} as row
MATCH (m:Adduct)-[r:has_iso]->(n:Adduct)
WHERE m.id = row.monoAddId and n.id = row.isoId
SET r.irObs = 0.0
"
map <- list() 
if(length(removeIRidx) > 0){
  for (i in 1:length(removeIRidx)){
    idx <- removeIRidx[i]
    map[[i]] <- list(monoAddId = obsIRmonoFomu[idx], isoId = obsIRisoFomu[idx])
  }
  cypher(graph, query, map = map)
}

newIRmonoFomu <- l$iso$newIRmonoFomu
newIRisoFomu <- l$iso$newIRisoFomu
newIRirValue <- l$iso$newIRirValue
map <- list() 
query <- "
UNWIND {map} as row
MATCH (m:Adduct)-[r:has_iso]->(n:Adduct)
WHERE m.id = row.monoAddId and n.id = row.isoId
SET r.irValue = row.irValue
"
if (length(newIRmonoFomu) > 0){
  for (i in 1:length(newIRmonoFomu)){
    map[[i]] <- list(monoAddId = newIRmonoFomu[i], isoId = newIRisoFomu[i], irValue = newIRirValue[i])
  }
  cypher(graph, query, map = map)
}

recIRdistr <- l$iso$recIRdistr
recIRdistrMonoFomu <- l$iso$recIRdistrMonoFomu
recIRdistrIsoFomu <- l$iso$recIRdistrIsoFomu
map <- list()
query <- "
UNWIND {map} as row
MATCH (m:Adduct)-[r:has_iso]->(n:Adduct)
WHERE m.id = row.monoAddId and n.id = row.isoId
SET r.irNum = row.irNum, r.irMean = row.irMean, r.irLogMean = row.irLogMean, r.irStd = row.irStd, r.irLogStd = row.irLogStd
"
if (length(recIRdistrMonoFomu) > 0){
  for (i in 1:length(recIRdistrMonoFomu)){
    irDistr <- recIRdistr[[i]]
    map[[i]] <- list(monoAddId = recIRdistrMonoFomu[i], isoId = recIRdistrIsoFomu[i], irNum = irDistr$iso_num, irMean = irDistr$iso_mean, irLogMean = irDistr$iso_logMean, irStd = irDistr$iso_std, irLogStd = irDistr$iso_logStd)
  }
  cypher(graph, query, map = map)
}


# update add connectivities (add-add, comp-add)
recCompMonoFomu <- l$add$recCompMonoFomu
recMITwithInSameComp <- l$add$recMITwithInSameComp
map <- list()
query <- "
UNWIND {map} as row
MATCH (n:Adduct)-[r:is_monoAdd]->(n:Adduct)
WHERE n.id = row.monoAddId
SET r.mostIntTime = row.mostIntTime
"
k <- 0
for (i in 1:length(recCompMonoFomu)){
  compMonoSet <- recCompMonoFomu[[i]]
  compMonoMITset <- recMITwithInSameComp[[i]]
  for (j in 1:length(compMonoSet)){
    k <- k + 1
    map[[k]] <- list(monoAddId = compMonoSet[j], mostIntTime = compMonoMITset[j])
  }
}
cypher(graph, query, map = map)

newMainAddCompIdx <- l$add$newMainAddCompIdx
newMainAddId <- l$add$newMainAddId
map1<- list()
query1 <- "
UNWIND {map1} as row
MATCH (m:Chemical)-[r:has_mainAdd]->(a:Adduct)
WHERE m.id = row.compId
DELETE r
"
map2<- list()
query2 <- "
UNWIND {map2} as row
MATCH (m:Chemical), (n:Adduct)
WHERE m.id = row.compId and n.id = row.mainAddId
CREATE (m)-[:has_mainAdd]->(n)
"
if(length(newMainAddCompIdx) > 0){
  for(i in 1:length(newMainAddCompIdx)){
    idx <- newMainAddCompIdx[i]
    newMainAdd <- newMainAddId[i]
    map1[[i]] <- list(compId = recComp[idx])
    map2[[i]] <- list(compId = recComp[idx], mainAddId = newMainAdd)
  }
  cypher(graph, query1, map1 = map1)
  cypher(graph, query2, map2 = map2)
}


# update compound retention time info
retainObsRTidx <- l$rt$retainObsRTidx
obsCompRT <- l$rt$obsCompRT
obsRTcomp <- l$rt$obsRTcomp
map <- list()
query <- "
UNWIND {map} as row
MATCH (n:Chemical)
WHERE n.id = row.compId
SET n.rtAll = row.rtAll
"
if (length(retainObsRTidx) > 0){
  for (i in 1:length(retainObsRTidx)){
    idx <- retainObsRTidx[i]
    map[[i]] <- list(compId = obsRTcomp[idx], rtAll = obsCompRT[idx])
  }
  cypher(graph, query, map = map)
}

recCompRT <- l$rt$recCompRT
recRTcomp <- l$rt$recRTcomp
map <- list()
query <- "
UNWIND {map} as row
MATCH (n:Chemical)
WHERE n.id = row.compId
SET n.rtAll = row.rtAll
"
if (length(recRTcomp) > 0){
  for (i in 1:length(recRTcomp)){
    map[[i]] <- list(compId = recRTcomp[idx], rtAll = recCompRT[idx])
  }
  cypher(graph, query, map = map)
}

recCompRTdistr <- l$rt$recCompRTdistr
recRTdistrComp <- l$rt$recRTdistrComp
map <- list()
query <- "
UNWIND {map} as row
MATCH (n:Chemical)
WHERE n.id = row.compId
SET n.rtMean = row.rtMean, n.rtStd = row.rtStd
"
if (length(recRTdistrComp) > 0){
  for (i in 1:length(recRTdistrComp)){
    rtDistr <- recCompRTdistr[[i]]
    map[[i]] <- list(compId = recRTdistrComp[i], rtMean = rtDistr$rt_mean, rtStd = rtDistr$rt_std)
  }
  cypher(graph, query, map = map)
}


# update mass accuracy info (integrated into an Experiment node)
query <- "
MATCH (n:Experiment)
RETURN max(n.id) as maxExprId
"
expr_maxId <- 0
maxExprId <- cypherToList(graph, query)
if (length(maxExprId[[1]]) != 0){
  expr_maxId <- maxExprId[[1]]$maxExprId
}
obsMA <- l$ma$obsMA
recMA <- l$ma$recMA
maNew <- l$ma$ma
ppmMean <- 0
ppmStd <- 0
recMAdistr <- l$ma$recMAdistr
if (length(obsMA) != 0){
  recMA <- obsMA
}
if (length(recMAdistr) != 0){
  ppmMean <- recMAdistr$ma_mean
  ppmStd <- recMAdistr$ma_std
}
query <- "CREATE (n:Experiment { id:{id}, ionMode:{ionMode}, ppmAll:{ppmAll}, ppmInit:{ppmInit}, ppmNew:{ppmNew}, ppmMean:{ppmMean}, ppmStd:{ppmStd}, 
                                 unknownP:{unknP}, limit:{limit}, rtWin:{rtWin}, it:{it}, burn:{burn}, delAdd:{delAdd}, delIso:{delIso}, delBio:{delBio}, ratioToll:{ratioToll}, 
                                 s_pkAdd:{s_pkAdd}, s_pkComp:{s_pkComp}, s_ma:{s_ma}, s_iso:{s_iso}, w_int:{w_int}, w_ma:{w_ma}, w_addLink:{w_addLink}, w_isoLink:{w_isoLink}, 
                                 t_rt:{t_rt}, t_post:{t_post}, massNum:{massNum}, dateRange:{dateRange}})"
# dateStart <- Sys.time()

cypher(graph, query, id = expr_maxId + 1, ionMode = ionMode, ppmAll = recMA, ppmInit = ma, ppmNew = maNew, ppmMean = ppmMean, ppmStd = ppmStd, 
       unknP = unknowP, limit = limit, rtWin = rtWin, it = it, burn = burn, delAdd = delAdd, delIso = delIso, delBio = delBio, ratioToll = ratioToll,
       s_pkAdd = s_pkAdd, s_pkComp = s_pkComp, s_ma = s_ma, s_iso = s_iso, w_int = w_int, w_ma = w_ma, w_addLink= w_addLink, w_isoLink = w_isoLink,
       t_rt = t_rt, t_post = t_post, massNum = massNum, dateRange = c(dateStart, Sys.time()))

query = "CREATE CONSTRAINT ON (n:Experiment) ASSERT n.id IS UNIQUE"
cypher(graph, query)
