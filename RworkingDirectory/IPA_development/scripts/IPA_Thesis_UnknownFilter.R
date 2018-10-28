# removes all objects from the current workspace (R memory)
# + various initialisations
rm(list = ls())
library(mzmatch.R)
library(RNeo4j)
library(slam)
library(Matrix)
mzmatch.init(version.1 = F)

Rcpp::sourceCpp('./RworkingDirectory/IPA_development/class/PriorRelated/IPA_Massbased_Priors_NoRT_u_Rcpp.cpp')
Rcpp::sourceCpp('./RworkingDirectory/IPA_development/class/UpdateRelated/IPA_Update.cpp')

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
limit <- 1e-02
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




sampleNum <- 1000
massMZ <- massMZ[1:sampleNum]
massRT <- massRT[1:sampleNum]
massInt <- massInt[1:sampleNum]


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
fomuPK <- vector('numeric')
fomuInfo <- cypherToList(graph, query, addId = finalAdds)
for (i in 1:length(fomuInfo)) {
  fomuMZ <- c(fomuMZ, fomuInfo[[i]]$fomuMass)
  fomuPK <- c(fomuPK, fomuInfo[[i]]$fomuPrior)
}
rm(fomuInfo)


resultU <- vector('integer')
resultTppm <- vector('numeric')
resultP <- vector('numeric')

for (u in 1:3){
  for (t_ppm in c(1,2,3,4,5,6,7,8,9,10)){
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
    
    finalAddsInit <- finalAdds
    massMZinit <- massMZ
    massRTinit <- massRT
    massIntInit <- massInt
    fomuMZinit <- fomuMZ
    fomuPKinit <- fomuPK
    
    ## Prior probability computing using Rcpp
    Prior <- NULL
    Prior <- ComputePriorRcppThesis(massMZinit,                    # mass values of masses
                                    fomuMZinit,                    # mass values of formulas
                                    fomuPKinit,                    # prior knowledge of formulas
                                    NA,                        # retention time of masses (not applied yet)
                                    ma,                        # mass accuracy
                                    NULL,                      # retention time list (not applied yet)
                                    u,               # value for unknown probability
                                    limit,                     # limit used when probability value is too small
                                    TRUE,                      # print log?
                                    FALSE                      # apply retention time into prior computing?
    )
    
    
    ## POSTERIOR COMPUTING
    Prior[!is.finite(Prior)] <- 0
    rowName <- as.character(seq(nrow(Prior))) 
    colName <- c(as.character(finalAddsInit), "unknown")
    rownames(Prior) <- rowName
    colnames(Prior) <- colName
    massKept <- which(Prior[,ncol(Prior)] != 1)              # filter out those masses 100% assigned with unknow
    fomuKept <- which(col_sums(Prior)>0)                     # and filter out those formulas which do not have any mass to assigned with
    prior_filtered <- Prior[massKept, fomuKept]       
    massRTinit <- massRTinit[massKept]
    massIntInit <- massIntInit[massKept]
    lastCol <- fomuKept[length(fomuKept)]
    if(lastCol == ncol(Prior)){                               # unknow has been kept
      finalAddsInit <- finalAddsInit[fomuKept[1:(length(fomuKept) - 1)]]
      newColName <- c(as.character(finalAddsInit), "unknown")
    } else{
      finalAddsInit <- finalAddsInit[fomuKept]
      newColName <- c(as.character(finalAddsInit))
    }
    colNameIdx <- finalAddsInit
    fomuNum <- length(fomuKept)
    
    
    # get valid compound set+mono adduct set
    queryStringFront <- "MATCH (m:Chemical)-"
    queryStringBack <- "->(n:Adduct)-[:is_monoAdd]->(n:Adduct) WHERE n.id in {addId} RETURN distinct m.id as compId, n.id as monoAddId order by m.id, n.id"
    query <- paste(queryStringFront, ionQueryString, queryStringBack, sep = "")
    monoComp_valiCompIds <- vector('integer')
    monoComp_valiMonoAddIds <- vector('integer')
    compMonoInfo <- cypherToList(graph, query, addId = finalAddsInit)
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
    addCompInfo <- cypherToList(graph, query, addId = finalAddsInit)
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
    add_monoMainInfo <- cypherToList(graph, query, monoAddId = monoComp_valiMonoAddIds, addId = finalAddsInit, compId = monoComp_valiCompIds)
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
    iso_pairInfo <- cypherToList(graph, query, monoAddId = monoComp_valiMonoAddIds, addId = {finalAddsInit})
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
    bio_compMainAddInfo <- cypherToList(graph, query, compId = monoComp_valiCompIds, addId = finalAddsInit)
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
    
    
    Posterior <- prior_filtered
    
    # deal with 'unknown' before update
    if(lastCol == ncol(Prior)){                           # unknow has been kept
      Posterior <- Posterior[,1:(ncol(Posterior) - 1)]
    }
    
    
    ## UPDATE STAGE 1: ALGORITHM 
    fomuPKinit <- fomuPKinit[fomuKept]
    massMZinit <- massMZinit[massKept]
    fomuMZinit <- fomuMZinit[fomuKept]
    
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
    
    
    l <- list()
    l <- UpdateMain(colNameIdx, Posterior, Iso, Add, fomuPKinit, compPK, ma, massMZinit, fomuMZinit, obsMA, recMA, recMAdistr,  
                    massIntInit, obsIR, obsIRmonoFomu, obsIRisoFomu, recIRdistr, recIRdistrMonoFomu, recIRdistrIsoFomu,
                    recMonoFomu, recAllAddFomu, recAllAddComp, recComp, recCompMainAddFomu, recCompMonoFomu, recMITwithInSameComp,
                    massRTinit, obsCompRT, obsRTcomp, recCompRT, recRTcomp, recCompRTdistr, recRTdistrComp,
                    s_pkAdd, s_pkComp, s_ma, s_iso, w_int, w_ma, w_addLink, w_isoLink, t_rt, t_post, TRUE)
    
    
    
    # get observed mass accuracy values after the update stage
    obsMA <- l$ma$obsMA
    recMA <- l$ma$recMA
    if (length(obsMA) != 0){
      recMA <- obsMA
    }
    
    
    # analyse observed mass accuracy values in recMA
    
    failToFilter <- which(recMA >= t_ppm || recMA <= -t_ppm)
    p <- length(failToFilter) / length(recMA)
    
    resultU <- c(resultU, u)
    resultTppm <- c(resultTppm, t_ppm)
    resultP <- c(resultP, p)
  }
}


plotSource <- data.frame(resultU, resultTppm, resultP)

write.table( x = plotSource, file = "C:\\Users\\oyyqwhuiss\\Desktop\\plotSource.csv", sep=",", col.names=FALSE, row.names = FALSE, append = FALSE)

