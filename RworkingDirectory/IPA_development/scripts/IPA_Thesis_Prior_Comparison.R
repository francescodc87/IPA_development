# removes all objects from the current workspace (R memory)
# + various initialisations
rm(list = ls())
library(mzmatch.R)
library(RNeo4j)
library(slam)
library(Matrix)
mzmatch.init(version.1 = F)
source('./RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_R.R')
Rcpp::sourceCpp('./RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_NoRT_Rcpp.cpp')

# PeakML.Data <- PeakML.Read("/home/yuqiouyang/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data <- PeakML.Read("C:/Users/oyyqwhuiss/Documents/IPA_development/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data.Table <- PeakML.Methods.getCompleteTable(PeakML.Data)
massMZ <- apply(PeakML.Data.Table$Masses, 2, mean, na.rm=TRUE)                 # mass values: Numeric Vector (1xM) (1x7842)
massRT <- as.numeric(apply(PeakML.Data.Table$Retentiontimes,2, mean, na.rm=TRUE))             # RT values: Numeric Vector (1xM) (1x7842)
massInt <- t(PeakML.Data.Table$Intensities) 
massInt <- as.numeric(apply(massInt,1, max, na.rm=T))         # intensity values: Numeric Vector (1xM) (1x7842)
rm(PeakML.Data)

ppm_query <- 5
ionMode <- "pos"
ma <- 5
u <- 3
limit <- 1e-02
if(ionMode == "pos"){
  ionQueryString <- "[:has_pos]"
} else{
  ionQueryString <- "[:has_neg]"
}


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
  RETURN n.ppmInit as initMA, n.ppmAll as storedMA, n.ppmMean as maDistrMean, n.ppmStd as maDistrStd
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

## Prior probability computing using R
system.time({ComputePriorR(mass = massMZ,
                        compMass = fomuMZ, 
                        ppm = ma,
                        unknownProb=0.05, 
                        compId=1:length(fomuMZ), 
                        v=T)})


## Prior probability computing using Rcpp
unknownProb = 0.05
system.time({ComputePriorRcpp(massMZ,                    # mass values of masses
                              fomuMZ,                    # mass values of formulas
                              fomuPk,                    # prior knowledge of formulas
                              NA,                        # retention time of masses (not applied yet)
                              ma,                        # mass accuracy
                              NULL,                      # retention time list (not applied yet)
                              unknownProb,                         # value for unknown probability
                              limit,                     # limit used when probability value is too small
                              TRUE,                      # print log?
                              FALSE                      # apply retention time into prior computing?
)})

# save(Prior_Rcpp, file="~/RworkingDirectory/IPA_development/data/Prior_Rcpp_oneEntry_180806.Rdata")









