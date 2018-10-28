# removes all objects from the current workspace (R memory)
# + various initialisations
rm(list = ls())
library(mzmatch.R)
library(RNeo4j)
library(slam)
library(Matrix)
mzmatch.init(version.1 = F)
# source('./RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_R.R')
Rcpp::sourceCpp('./RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_NoRT_u_Rcpp.cpp')
# Rcpp::sourceCpp('./RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_NoRT_Rcpp.cpp')

# source('./RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Multsample_R.R')
# source('./RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Compute_post.R')
# source('./RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_R.R')
# source('./RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_INT_NoPot.R')

source('./RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_RT.R')
source('./RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_relID.R')
source('./RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_corr.R')
source('./RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_INT_NoPot_Rcpp.R')
Rcpp::sourceCpp('./RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_GS_INT_NoPot.cpp')   #'GS' = 'GibbsSampling'

# PeakML.Data <- PeakML.Read("/home/yuqiouyang/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data <- PeakML.Read("C:/Users/oyyqwhuiss/Documents/IPA_development/RworkingDirectory/IPA_development/data/allpeaks_filtered.peakml")
PeakML.Data.Table <- PeakML.Methods.getCompleteTable(PeakML.Data)
massMZ <- apply(PeakML.Data.Table$Masses, 2, mean, na.rm=TRUE)                 # mass values: Numeric Vector (1xM) (1x7842)
massRT <- as.numeric(apply(PeakML.Data.Table$Retentiontimes,2, mean, na.rm=TRUE))             # RT values: Numeric Vector (1xM) (1x7842)
massInt <- t(PeakML.Data.Table$Intensities) 
massInt <- as.numeric(apply(massInt,1, max, na.rm=T))         # intensity values: Numeric Vector (1xM) (1x7842)
rm(PeakML.Data)

sampleNum <- 100
massMZ <- massMZ[1:sampleNum]
massRT <- massRT[1:sampleNum]
massInt <- massInt[1:sampleNum]

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

rtWin <- 5
delAdd <- 0.4
delIso <- 0.2
delBio <- 1
ratioToll<- 0.8



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



## Prior probability computing using Rcpp
u = 3
Prior <- NULL
Prior <- ComputePriorRcppThesis(massMZ,                    # mass values of masses
                                fomuMZ,                    # mass values of formulas
                                fomuPk,                    # prior knowledge of formulas
                                NA,                        # retention time of masses (not applied yet)
                                ma,                        # mass accuracy
                                NULL,                      # retention time list (not applied yet)
                                u,               # value for unknown probability
                                limit,                     # limit used when probability value is too small
                                TRUE,                      # print log?
                                FALSE                      # apply retention time into prior computing?
)
# 
# unknownProb = 0.05
# Prior2 <- NULL
# Prior2 <- ComputePriorRcpp(massMZ,                    # mass values of masses
#                           fomuMZ,                    # mass values of formulas
#                           fomuPk,                    # prior knowledge of formulas
#                           NA,                        # retention time of masses (not applied yet)
#                           ma,                        # mass accuracy
#                           NULL,                      # retention time list (not applied yet)
#                           unknownProb,                         # value for unknown probability
#                           limit,                     # limit used when probability value is too small
#                           TRUE,                      # print log?
#                           FALSE                      # apply retention time into prior computing?
# )

## POSTERIOR COMPUTING
Prior[!is.finite(Prior)] <- 0
rowName <- as.character(seq(nrow(Prior))) 
colName <- c(as.character(finalAdds), "unknown")
rownames(Prior) <- rowName
colnames(Prior) <- colName
massKept <- which(Prior[,ncol(Prior)] != 1)              # filter out those masses 100% assigned with unknow
fomuKept <- which(col_sums(Prior, na.rm = TRUE)>0)        # and filter out those formulas which do not have any mass to assigned with
prior_filtered <- Prior[massKept, fomuKept]       
massRT <- massRT[massKept]
massInt <- massInt[massKept]
lastCol <- fomuKept[length(fomuKept)]
if(lastCol == ncol(Prior)){                               # unknow has been kept
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

# construct adduct matrix
queryStringFront <- "MATCH (a:Adduct)<-[:is_monoAdd]-(a:Adduct)<-"
queryStringMid <- "-(m:Chemical)-"
queryStringBack <- "->(b:Adduct), (m:Chemical)-[:has_mainAdd]->(b:Adduct) WHERE a.id in {monoAddId} and b.id in {addId} and a.id <> b.id and m.id in {compId} RETURN distinct a.id as monoAddId, b.id as mainAddId order by a.id"
query <- paste(queryStringFront, ionQueryString, queryStringMid, ionQueryString, queryStringBack, sep = "")
add_monoAddIds <- vector('integer')
add_mainAddIds <- vector('integer')
add_monoMainInfo <- cypherToList(graph, query, monoAddId = monoComp_valiMonoAddIds, addId = finalAdds, compId = monoComp_valiCompIds)
if (length(add_monoMainInfo) != 0){
  for (i in 1:length(add_monoMainInfo)) {
    add_monoAddIds <- c(add_monoAddIds, add_monoMainInfo[[i]]$monoAddId)
    add_mainAddIds <- c(add_mainAddIds, add_monoMainInfo[[i]]$mainAddId)
  }
}
Add <- Matrix(data = 0, nrow = fomuNum, ncol = fomuNum, sparse = TRUE, dimnames = list(newColName, newColName))
if (length(add_monoAddIds) != 0){
  for (i in 1:length(add_monoAddIds)){
    rowIdxString <- as.character(add_monoAddIds[i])
    colIdxString <- as.character(add_mainAddIds[i])
    Add[rowIdxString, colIdxString] <- 1
    Add[colIdxString, rowIdxString] <- 1
  }
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
if (length(iso_pairInfo) != 0){
  for (i in 1:length(iso_pairInfo)) {
    iso_monoAddIds <- c(iso_monoAddIds, iso_pairInfo[[i]]$monoAddId)
    iso_isoIds <- c(iso_isoIds, iso_pairInfo[[i]]$isoId)
    iso_irs <- c(iso_irs, iso_pairInfo[[i]]$irMonoIso)
  }
}

Iso <- Matrix(data = 0, nrow = fomuNum, ncol = fomuNum, sparse = TRUE, dimnames = list(newColName, newColName))
if (length(iso_monoAddIds) != 0){
  for(i in 1:length(iso_monoAddIds)){
    rowIdxString <- as.character(iso_monoAddIds[i])
    colIdxString <- as.character(iso_isoIds[i])
    irValue <- iso_irs[i]
    Iso[rowIdxString, colIdxString] <- irValue
    Iso[colIdxString, rowIdxString] <- 1 / irValue
  }
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
if(length(bio_compMainAddInfo) != 0){
  for (i in 1:length(bio_compMainAddInfo)) {
    bio_compIds <- c(bio_compIds, bio_compMainAddInfo[[i]]$compId)
    bio_compMainAddIds <- c(bio_compMainAddIds, bio_compMainAddInfo[[i]]$mainAddId)
  }
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
if(length(bio_compReacInfo) != 0){
  for (i in 1:length(bio_compReacInfo)) {
    bio_reac_compIds <- c(bio_reac_compIds, bio_compReacInfo[[i]]$reacCompId)
    bio_reac_linkerCompIds <- c(bio_reac_linkerCompIds, bio_compReacInfo[[i]]$reacLinkerCompId)
  }
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
if (length(bio_compCofInfo) != 0){
  for (i in 1:length(bio_compCofInfo)) {
    bio_cof_compIds <- c(bio_cof_compIds, bio_compCofInfo[[i]]$cofCompId)
    bio_cof_linkerCompIds <- c(bio_cof_linkerCompIds, bio_compCofInfo[[i]]$cofLinkerCompId)
  }
}
rm(bio_compCofInfo)

# last, fill the bio matrix
Bio <- Matrix(data = 0, nrow = fomuNum, ncol = fomuNum, sparse = TRUE, dimnames = list(newColName, newColName))
if(length(bio_reac_compIds) != 0){
  for(i in 1:length(bio_reac_compIds)){
    bio_reac_comp <- bio_reac_compIds[i]
    bio_reac_linkcomp <- bio_reac_linkerCompIds[i]
    bio_reac_compMainAdd <- bio_compMainAddIds[which(bio_compIds == bio_reac_comp)]
    bio_reac_linkCompMainAdd <- bio_compMainAddIds[which(bio_compIds == bio_reac_linkcomp)]
    Bio[as.character(bio_reac_compMainAdd),as.character(bio_reac_linkCompMainAdd)] <- 1
  }
}


if (length(bio_cof_compIds) != 0){
  for(i in 1:length(bio_cof_compIds)){
    bio_cof_comp <- bio_cof_compIds[i]
    bio_cof_linkcomp <- bio_cof_linkerCompIds[i]
    bio_cof_compMainAdd <- bio_compMainAddIds[which(bio_compIds == bio_cof_comp)]
    bio_cof_linkCompMainAdd <- bio_compMainAddIds[which(bio_compIds == bio_cof_linkcomp)]
    Bio[as.character(bio_cof_compMainAdd),as.character(bio_cof_linkCompMainAdd)] <- 1
  }
}

# system.time({ComputePosteriorR_Add_Iso_Bio_Int_NoPot(P=prior_filtered2,Add = Add,
#                                                      Iso = Iso, RT = massRT,
#                                                      Bio = Bio, Int = massInt,
#                                                      RTwin = rtWin, it = it, burn = burn,
#                                                      delAdd = delAdd, delIso = delIso, delBio = delBio,
#                                                      ratioToll = ratioToll, allSamp = F, log = T)})

itMax <- 600

resultIt <- vector('integer')
resultSum <- vector('numeric')
for (it in 1:itMax){
  cat("it: ", it)
  burn  = it %/% 10
  Posterior <- ComputePosteriorRcpp_Add_Iso_Bio_Int_NoPot(P = prior_filtered, Add = Add,
                                                          Iso = Iso, RT = massRT,
                                                          Bio = Bio, Int = massInt,
                                                          RTwin = rtWin, it = it, burn = burn,
                                                          delAdd = delAdd, delIso = delIso, delBio = delBio,
                                                          ratioToll = ratioToll, log = T)
  vecSum <- 0
  for (j in 1:nrow(Posterior)) {
    postVec <- Posterior[j,]
    validPostVec <- postVec[postVec>0]
    vecSum <- vecSum - sum(log(validPostVec))
  }
  resultSum <- c(resultSum, vecSum)
  resultIt <- c(resultIt, it)
}


plotSource <- data.frame(resultIt, resultSum)


write.table( x = plotSource, file = "C:\\Users\\oyyqwhuiss\\Desktop\\plotSource.csv", sep=",", col.names=FALSE, row.names = FALSE, append = FALSE)
# plot(resultIt, resultSum, type="b", xlab="it", ylab="value sum")





