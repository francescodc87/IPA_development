rm(list = ls())
library(enviPat)
library(RNeo4j)
library(data.table)

# CAUTION THIS SCRIPT IS FOR THE INITIAL BIOCHEM4J DATABASE, RUN THIS AFTER YOU LOAD THE BASIC BIOCHEM4J DATABASE EXAMPLE

# here to connect with database on local server
graph = startGraph("http://localhost:7474/db/data/", username = "neo4j", password = "oyyq6997")

## here is some initial change on database, !!!no need any more!!!
## delete the nodes in Chemical with property formula is empty or NA
## delete the nodes in Chemical with no monoisotopic mass value
# query = "MATCH (n:Chemical) WHERE not exists(n.formula) DETACH DELETE n"
# cypher(graph, query)
# query = "MATCH (n:Chemical) WHERE n.formula = 'NA' DETACH DELETE n"
# cypher(graph, query)
# query = "MATCH (n:Chemical) WHERE not exists(n.monoisotopic_mass) DETACH DELETE n"
# cypher(graph, query)


# deal with dbId and id property in Chemical node
# query = "MATCH (n:Chemical) SET n.dbId = n.id REMOVE n.id"
# cypher(graph, query)
# query = "MATCH (n:Chemical) SET n.id = id(n) + 1"
# cypher(graph, query)

## add a constraint on Chemical node (function as an index), this is important
# query = "CREATE CONSTRAINT ON (n:Chemical) ASSERT n.id IS UNIQUE"
# cypher(graph, query)
##or addConstraint(graph, "Chemical", "id")


## next is to populate the data into the database, the jobs are:
# 1. add prior knowledge ("pk") property for Chemical node, initial value all equal to 1 & add properties (rtAll,rtMean,rtStd) for Chemical node, initialise value to
# query = "MATCH (n:Chemical) SET n.pk = 1, n.rtAll = '', n.rtMean = 0, n.rtStd = 0"
# cypher(graph, query)

# 2. link compounds to adducts & isotopes, add properties for adduct node
# get chemical Information split in multi-times for memory efficiency
# query = "MATCH (n:Chemical) RETURN n.id AS id, n.formula AS formula, n.monoisotopic_mass as mz ORDER BY n.id limit 100000"
# query = "MATCH (n:Chemical) RETURN n.id AS id, n.formula AS formula, n.monoisotopic_mass as mz ORDER BY n.id skip 100000 limit 100000"
# query = "MATCH (n:Chemical) RETURN n.id AS id, n.formula AS formula, n.monoisotopic_mass as mz ORDER BY n.id skip 200000 limit 100000"
# query = "MATCH (n:Chemical) RETURN n.id AS id, n.formula AS formula, n.monoisotopic_mass as mz ORDER BY n.id skip 300000 limit 100000"
# query = "MATCH (n:Chemical) RETURN n.id AS id, n.formula AS formula, n.monoisotopic_mass as mz ORDER BY n.id skip 400000 limit 100000"
# query = "MATCH (n:Chemical) RETURN n.id AS id, n.formula AS formula, n.monoisotopic_mass as mz ORDER BY n.id skip 500000 limit 100000"
# query = "MATCH (n:Chemical) RETURN n.id AS id, n.formula AS formula, n.monoisotopic_mass as mz ORDER BY n.id skip 600000"
chemInfo <- cypher(graph, query)

# chemInfo is a list, turn list into a datatable
chemInfoTable <- setDT(chemInfo)
rm(chemInfo)
# get valid compounds (checked by check_chemform)
data("adducts")
data("isotopes")

# first time filter about formula NA
validFomuNAIdx <- chemInfoTable[formula != "NA", which = TRUE]
chemInfoTable <- chemInfoTable[validFomuNAIdx,]
# second time filter about mono mass value
validMassNAIdx <- chemInfoTable[!is.na(mz), which = TRUE]
chemInfoTable <- chemInfoTable[validMassNAIdx,]
# third time filter about formula valid to isopattern
checkResult <- check_chemform(isotopes,chemInfoTable$formula)
validFomuIdx <- which(checkResult[,"warning"] == FALSE)  # should filter out formula "NA" and empty
chemInfoTable <- chemInfoTable[validFomuIdx,]        
chemInfoTable$formula <- checkResult[validFomuIdx,"new_formula"]

rm(validFomuNAIdx)
rm(validMassNAIdx)
rm(checkResult)
rm(validFomuIdx)

addType <- c("M+H","M+Na","M+2H","M-H","M+Cl","M-2H")
addRowIdx <- vector('numeric')
isoPickRule <- c(2,3,4)
chemNum <- length(chemInfoTable$formula)
for(i in 1: length(addType)){
  addRowIdx[i] <- which(adducts$Name == addType[i])
}
addId <- 0
v <- 0
addTable = data.table(id = rep(0, 3e6), name = rep("",3e6), formula = rep("",3e6), type = rep("", 3e6), mz = rep(0, 3e6),
                      mainAddOf = rep(0,3e6), isPOS = rep(-1,3e6), POSorNEGof = rep(0,3e6), isoOf = rep(0,3e6), isoRank = rep(0,3e6), initIR = rep(0,3e6))
for(i in 1:chemNum){
  # for each compound, get all the monoisotopic adduct firstly
  eachChemFomu <- chemInfoTable$formula[i]
  compId <- chemInfoTable$id[i]
  for(j in 1:length(addRowIdx)){
    # cat("1")
    eachAddFomu <- ""
    idx <- addRowIdx[j]
    if(j == 5){
      eachAddFomu <- mergeform(eachChemFomu, adducts[idx, "Formula_add"])     # the formula of the adduct
      flagPOS <- 0  # deal with has_pos relation
    }else if(adducts[idx, "Ion_mode"] == "positive"){
      # cat("1.2")
      eachAddFomu <- mergeform(eachChemFomu, adducts[idx, "Formula_add"])     # the formula of the adduct
      flagPOS <- 1   # deal with has_pos relation
    }else{
      # cat("1.3")
      Hidx <- gregexpr(pattern ='H',eachChemFomu)[[1]][1]
      if(Hidx >= 1){
        Hnum <- substr(eachChemFomu, start = Hidx + 1, stop = Hidx + 1)
        if(Hnum != "e" && Hnum != "o" && Hnum != "f" && Hnum != "g"){
          Hnum <- as.integer(Hnum) 
          # cat("Hnum: ",Hnum)
          if(Hnum >= as.integer(substr(adducts[idx, "Formula_ded"], start = 2, stop = 2))){
            eachAddFomu <- subform(eachChemFomu, adducts[idx, "Formula_ded"])       # the formula of the adduct
            if(eachAddFomu == "NANA"){
              eachAddFomu = ""
            }
            flagPOS <- 0  # deal with has_neg relation
          }
        }
      }
    }
    # cat("1.5")
    if(eachAddFomu != ""){
      addId <- addId + 1
      eachAddMass <- chemInfoTable$mz[i] + adducts[idx, "Mass"] # the mass of the adduct
      eachAddCharge <- adducts[idx, "Charge"]           # the charge of the adduct
      eachAddMZ <- eachAddMass / abs(eachAddCharge)             # the mz of the adduct
      eachAddType <- addType[j]                            # the type of the adduct
      if(j == 1 || j == 4){ # deal with has_mainAdd relation
        mainAddComp <- compId
      }else{
        mainAddComp <- 0
      }
      eachAddName <- eachAddFomu
      eachAddId <- addId                                # the id of the adduct
      # cat("1.8")
      addTable[addId, c("id","name","formula","type","mz","mainAddOf","isPOS","POSorNEGof") := .(eachAddId, eachAddName, eachAddFomu, eachAddType, eachAddMZ, mainAddComp, flagPOS, compId)]
      # addTable[addId,] <- list("id" = eachAddId, "name" = eachAddName, "formula" = eachAddFomu, "type" = eachAddType, "mz" = eachAddMZ, 
      #                          "mainAddOf" = mainAddOf, "isPOS" = isPOS, "POSorNEGof" = compId, "isoOf" = 0, "isoRank" = 0, "initIR" = 0)
      # cat("2")
      # here start to search for its isotopes
      isoThrld <- 5
      allIsoInfo <- (isopattern(isotopes, eachAddFomu, threshold = isoThrld, charge = eachAddCharge, verbose = FALSE))
      allIsoInfo <- allIsoInfo[[1]]
      if(allIsoInfo[1] != "error"){
        int <- allIsoInfo[, "abundance"]
        while ((length(int) < 4) && (isoThrld != 0)) {
          isoThrld <- isoThrld / 2
          allIsoInfo <- (isopattern(isotopes, eachAddFomu, threshold = isoThrld, charge = eachAddCharge, verbose = FALSE))
          allIsoInfo <- allIsoInfo[[1]]
          int <- allIsoInfo[, "abundance"]
        }
        if(length(int) >= 4){
          intRank <- rank(-int)
          for (k in 1:length(isoPickRule)){
            # cat("3")
            addId <- addId + 1
            rk <- isoPickRule[k] 
            isoIdx <- which(intRank == rk)
            eachIsoMZ <- allIsoInfo[isoIdx, "m/z"]              # mz of the isotope
            eachIsoIR <- 100 / allIsoInfo[isoIdx, "abundance"]  # intensity ratio: calculation: monoisotope intensity / its isotope's intensity
            l <-  ncol(allIsoInfo)
            isoFomuRes <- which(allIsoInfo[isoIdx, 3:l] != 0)
            isoFomuResName <- names(isoFomuRes)
            eachIsoFomuEle <- ""
            for(o in 1:length(isoFomuResName)){
              eleString <- ""
              nameIdx <- isoFomuResName[o]
              eleNum <- allIsoInfo[isoIdx, nameIdx]
              eleString <- paste(nameIdx, ":", toString(eleNum))
              if(eachIsoFomuEle == ""){
                eachIsoFomuEle <- paste(eachIsoFomuEle, eleString)
              }else{
                eachIsoFomuEle <- paste(eachIsoFomuEle, eleString, sep = ", ")
              }
            }
            eachIsoFomu <- paste(eachAddFomu, eachIsoFomuEle, sep = "; ")         # formula of the isotope 
            eachIsoName <- eachIsoFomu                                            # name of the isotope
            eachIsoId <- addId                                                    # id of the isotope
            eachIsoRank <- rk - 1                                                 # rank of the isotope
            addTable[addId, c("id","name","formula","type","mz","isPOS","POSorNEGof", "isoOf", "isoRank", "initIR") := .(eachIsoId, eachIsoName, eachIsoFomu, eachAddType, eachIsoMZ, flagPOS, compId, eachAddId, eachIsoRank, eachIsoIR)]
            # addTable[addId,] <- list("id" = eachIsoId, "name" = eachIsoName, "formula" = eachIsoFomu, "type" = eachAddType, "mz" = eachIsoMZ, 
            #                          "mainAddOf" = 0, "isPOS" = isPOS, "POSorNEGof" = compId, "isoOf" = eachAddId, "isoRank" = rk - 1, "initIR" = eachIsoIR)
            # cat("4")
          }
        }
      }
    }
  }
  old_v <- v
  v <- (i * 1000) %/% chemNum
  if(old_v !=v){
    cat("DataTable generating, ",paste0(round(v, 0), "‰", "\n"))
    if(v %% 100 == 0){
      cat("\014")
    }
  }
}

## store the information into csv files
# save(addTable, file="C:/Users/oyyqwhuiss/Desktop/data_generate_1stpart.Rdata")
# save(addTable, file="C:/Users/oyyqwhuiss/Desktop/data_generate_2ndpart.Rdata")
# save(addTable, file="C:/Users/oyyqwhuiss/Desktop/data_generate_3rdpart.Rdata")
# save(addTable, file="C:/Users/oyyqwhuiss/Desktop/data_generate_4thpart.Rdata")
# save(addTable, file="C:/Users/oyyqwhuiss/Desktop/data_generate_5thpart.Rdata")
# save(addTable, file="C:/Users/oyyqwhuiss/Desktop/data_generate_6thpart.Rdata")
# save(addTable, file="C:/Users/oyyqwhuiss/Desktop/data_generate_7thpart.Rdata")


# delete empty space in datatable and re-arrange the id, then write into .csv file
fileURL <- c("C:/Users/oyyqwhuiss/Desktop/data_generate_1stpart.Rdata",
             "C:/Users/oyyqwhuiss/Desktop/data_generate_2ndpart.Rdata",
             "C:/Users/oyyqwhuiss/Desktop/data_generate_3rdpart.Rdata",
             "C:/Users/oyyqwhuiss/Desktop/data_generate_4thpart.Rdata",
             "C:/Users/oyyqwhuiss/Desktop/data_generate_5thpart.Rdata",
             "C:/Users/oyyqwhuiss/Desktop/data_generate_6thpart.Rdata",
             "C:/Users/oyyqwhuiss/Desktop/data_generate_7thpart.Rdata")

startId <- 1
endId <- 0
for(i in 1:length(fileURL)){
  load(file = fileURL[i])
  validId <- addTable[id != 0, which = TRUE]
  idAdd <- validId[length(validId)]
  startId <- endId + 1
  endId <- endId + idAdd
  # id: start~end
  newId <- c(startId : endId)
  newAddTable <- addTable[1:idAdd,]
  rm(addTable)
  newAddTable[,"id"] <- newId
  csvKeyId <- substr(fileURL[i], start = 43, stop = 43)
  csvURL <- paste("C:/Users/oyyqwhuiss/Desktop/data_generate_", csvKeyId, "part.csv", sep = "")
  fwrite(newAddTable, file = csvURL)
}








## next job is to use neo4j csv bulk import to insert the data into the database
# move all csv files in the <NEO4j_HOME>/import/ directory
# then the cypher command example:
# USING PERIODIC COMMIT
# LOAD CSV WITH HEADERS FROM "file:///data_generate_1part.csv" AS row
# CREATE (n:Adduct { id: toInteger(row.id), srcDb: 'enviPat', name: row.name, formula: row.formula, type: row.type, mz: toFloat(row.mz), pk: 1 })

# # create index for following data populating
# # before that we should already have unique constraint on Chemical:id
# CREATE CONSTRAINT ON (n:Adduct) ASSERT n.id IS UNIQUE

# # after indexes applied, create relations
# # 1. link between Chemical and Adduct node:
# # id, name, formula, type, mz, mainAddOf, isPOS, POSorNEGof, isoOf, isoRank, initIR
# USING PERIODIC COMMIT
# LOAD CSV WITH HEADERS FROM "file:///data_generate_1part.csv" AS row
# MATCH (a:Chemical),(b:Adduct) WHERE a.id = toInteger(row.mainAddOf) AND b.id = toInteger(row.id)
# CREATE (a)-[:has_mainAdd]->(b)
# # link has_pos, has_neg
# USING PERIODIC COMMIT
# LOAD CSV WITH HEADERS FROM "file:///data_generate_1part.csv" AS row
# MATCH (a:Chemical),(b:Adduct) WHERE a.id = toInteger(row.POSorNEGof) AND b.id = toInteger(row.id) AND toInteger(row.isPOS) = 1
# CREATE (a)-[:has_pos]->(b)
# USING PERIODIC COMMIT
# LOAD CSV WITH HEADERS FROM "file:///data_generate_1part.csv" AS row
# MATCH (a:Chemical),(b:Adduct) WHERE a.id = toInteger(row.POSorNEGof) AND b.id = toInteger(row.id) AND toInteger(row.isPOS) = 0
# CREATE (a)-[:has_neg]->(b)
# # 2. link inside Adduct node
# # link is_monoAdd{mostIntTime}
# USING PERIODIC COMMIT
# LOAD CSV WITH HEADERS FROM "file:///data_generate_1part.csv" AS row
# MATCH (a:Adduct) WHERE a.id = toInteger(row.id) AND toInteger(row.isoOf) = 0
# CREATE (a)-[:is_monoAdd{mostIntTime:0}]->(a)
# # update r.mostIntTime to 1 if main adduct
# USING PERIODIC COMMIT
# LOAD CSV WITH HEADERS FROM "file:///data_generate_1part.csv" AS row
# MATCH (a:Adduct)-[r:is_monoAdd]->(a:Adduct) WHERE a.id = toInteger(row.id) AND toInteger(row.mainAddOf) <> 0
# SET r.mostIntTime = 1
# # link has_iso {isoRank,irInit,irMean,irStd}
# USING PERIODIC COMMIT
# LOAD CSV WITH HEADERS FROM "file:///data_generate_1part.csv" AS row
# MATCH (a:Adduct),(b:Adduct) WHERE toInteger(row.isoOf) <> 0 AND a.id = toInteger(row.isoOf) AND b.id = toInteger(row.id)
# CREATE (a)-[:has_iso{isoRank:toInteger(row.isoRank), irInit:toFloat(row.initIR), irMean:toFloat(0), irStd:toFloat(0)}]->(b)


# # build  some indexes for using
# query = "CREATE INDEX ON :Adduct(mz)"
# cypher(graph, query)


# #after populate all data into database, remember keep a backup of the database!!! 
