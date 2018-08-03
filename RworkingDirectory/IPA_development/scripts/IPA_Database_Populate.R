rm(list = ls())
library(enviPat)
library(RNeo4j)

# CAUTION THIS SCRIPT IS FOR THE INITIAL BIOCHEM4J DATABASE, RUN THIS AFTER YOU LOAD THE BIOCHEM4J DATABASE

# here to connect with database on local server
graph = startGraph("http://localhost:7474/db/data/", username = "neo4j", password = "oyyq6997")

# here is some initial change on database
## ???delete the nodes in Chemical with property formula is empty or NA
## ???delete the nodes in Chemical with no monoisotopic mass value
# query = "MATCH (n:Chemical) WHERE not exists(n.formula) DETACH DELETE n"
# cypher(graph, query)
# query = "MATCH (n:Chemical) WHERE n.formula = 'NA' DETACH DELETE n"
# cypher(graph, query)
# query = "MATCH (n:Chemical) WHERE not exists(n.monoisotopic_mass) DETACH DELETE n"
# cypher(graph, query)

# deal with dbId and id property in Chemical node
# query = "MATCH (n:Chemical) SET n.dbId = n.id REMOVE n.id"
# cypher(graph, query)
# query = "MATCH (n:Chemical) SET n.id = id(n)"
# cypher(graph, query)




# next is to populate the data into the database, the jobs are:
#   1. add prior knowledge ("pk") property for Chemical node, initial value all equal to 1 & add properties (rtAll,rtMean,rtStd) for Chemical node, initialise value to
# query = "MATCH (n:Chemical) SET n.pk = 1, n.rtAll = '', n.rtMean = 0, n.rtStd = 0"
# cypher(graph, query)

#   2. link compounds to adducts & isotopes, add properties for adduct node
#   get chemical Information
query = "MATCH (n:Chemical) RETURN n.formula AS formula, n.monoisotopic_mass as mz ORDER BY n.id LIMIT 1"
chemInfo <- cypher(graph, query)
#   ???get valid compounds (checked by check_chemform)
validChemFomu <- chemInfo$formula
data("adducts")
data("isotopes")
addType <- c("M+H","M+Na","M+2H","M-H","M+Cl","M-2H")
addRowIdx <- vector('numeric')
isoPickRule <- c(2,3,4)
for(i in 1: length(addType)){
  addRowIdx[i] <- which(adducts$Name == addType[i])
}
addId <- 0
for(i in 1:length(validChemFomu)){
  # for each compound, get all the monoisotopic adduct firstly
  eachAddFomu <- ""
  eachChemFomu <- validChemFomu[i]
  for(j in 1:length(addRowIdx)){
    idx <- addRowIdx[j]
    if(j == 5){
      eachAddFomu <- mergeform(eachChemFomu, adducts[idx, "Formula_add"])     # the formula of the adduct
      isPOS <- FALSE  # deal with has_pos relation
    }else if(adducts[idx, "Ion_mode"] == "positive"){
      eachAddFomu <- mergeform(eachChemFomu, adducts[idx, "Formula_add"])     # the formula of the adduct
      isPOS <- TRUE   # deal with has_pos relation
    }else{
      Hidx <- gregexpr(pattern ='H',eachChemFomu)[[1]][1]
      Hnum <- as.integer(substr(eachChemFomu, start = Hidx + 1, stop = Hidx + 1)) 
      if(Hnum >= as.integer(substr(adducts[idx, "Formula_ded"], start = 2, stop = 2)))
        eachAddFomu <- subform(eachChemFomu, adducts[idx, "Formula_ded"])       # the formula of the adduct
        isPOS <- FALSE  # deal with has_neg relation
    }
    if(eachAddFomu != ""){
      addId <- addId + 1
      eachAddMass <- chemInfo$mz[i] + adducts[idx, "Mass"] # the mass of the adduct
      eachAddCharge <- abs(adducts[idx, "Charge"])           # the charge of the adduct
      eachAddMZ <- eachAddMass / eachAddCharge             # the mz of the adduct
      eachAddType <- addType[j]                            # the type of the adduct
      if(j == 1 || j == 4){ # deal with has_mainAdd relation
        isMainAdd = TRUE
      }else{
        isMainAdd = FALSE
      }
      # eachAddName <- paste(eachChemFomu, eachAddType, sep = " ")                  # the name of the adduct
      eachAddName <- eachAddFomu
      eachAddId <- addId                                # the id of the adduct
      
      # here to insert monoisotopic adduct information into adduct node
      query = sprintf("CREATE (n:Adduct { id: %s, srcDb: 'enviPat', name: '%s', formula: '%s', type: '%s', mz: %s, pk: 1 })",
                      eachAddId, eachAddName, eachAddFomu, eachAddType, eachAddMZ)
      cypher(graph, query)
      # here to link the ralation between compounds and monoisotopic adducts + link is_monoAdd self-relation on one mono adduct 
      if(isTRUE(isPOS)){
        if(isTRUE(isMainAdd)){
          query = sprintf("MATCH (a:Chemical),(b:Adduct) WHERE a.id = %s AND b.id = %s CREATE (a)-[:has_pos]->(b), (a)-[:has_mainAdd]->(b), (b)-[:is_monoAdd]->(b)",
                          i, eachAddId)
        }else{
          query = sprintf("MATCH (a:Chemical),(b:Adduct) WHERE a.id = %s AND b.id = %s CREATE (a)-[:has_pos]->(b), (b)-[:is_monoAdd]->(b)",
                          i, eachAddId)
        }
      }else{
        if(isTRUE(isMainAdd)){
          query = sprintf("MATCH (a:Chemical),(b:Adduct) WHERE a.id = %s AND b.id = %s CREATE (a)-[:has_neg]->(b), (a)-[:has_mainAdd]->(b), (b)-[:is_monoAdd]->(b)",
                          i, eachAddId)
        }else{
          query = sprintf("MATCH (a:Chemical),(b:Adduct) WHERE a.id = %s AND b.id = %s CREATE (a)-[:has_neg]->(b), (b)-[:is_monoAdd]->(b)",
                          i, eachAddId)
        }
      }
      cypher(graph, query)
      
      # here start to search for its isotopes
      isoThrld <- 5
      allIsoInfo <- (isopattern(isotopes, eachAddFomu, threshold = isoThrld, charge = eachAddCharge))
      allIsoInfo <- allIsoInfo[[1]]
      int <- allIsoInfo[, "abundance"]
      while (length(int) < 4) {
        isoThrld <- isoThrld / 2
        allIsoInfo <- (isopattern(isotopes, eachAddFomu, threshold = isoThrld, charge = eachAddCharge))
        allIsoInfo <- allIsoInfo[[1]]
        int <- allIsoInfo[, "abundance"]
      }
      intRank <- rank(-int)
      for (k in 1:length(isoPickRule)){
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
        # here to insert isotope information into adduct
        query = sprintf("CREATE (n:Adduct { id: %s, srcDb: 'enviPat', name: '%s', formula: '%s', type: '%s', mz: %s, pk: 1 })",
                        eachIsoId, eachIsoName, eachIsoFomu, eachAddType, eachIsoMZ)
        cypher(graph, query)  
        
        # here to link the isotope ralation between a mono adduct and a isotope
        query = sprintf("MATCH (a:Adduct),(b:Adduct) WHERE a.id = %s AND b.id = %s CREATE (a)-[:has_iso { isoRank: %s, irInit: %s, irMean: %s, irStd: %s}]->(b)",
                        eachAddId, eachIsoId, rk - 1, eachIsoIR, 0, 0)
        cypher(graph, query) 
      }
    }
  }
}



# query = "MATCH (n:Chemical) RETURN id(n) AS id, n.pk AS pk, n.formula AS formula ORDER BY id"
# chemInfo <- cypher(graph, query)

# limit <- 1000
# n <- chemNum %/% limit
#   however, there are near 800,000 nodes belong to chemical node, the memory cannot handle such amout of data, thus we need "skip" and "limit" to help with [in SQL it's called OFFSET] 
# if (n == 0){
#   query = "MATCH (n:Chemical) RETURN id(n) AS id, n.pk AS pk, n.formula AS formula ORDER BY id"
#   chemInfo <- cypher(graph, query)
# }else{
#   skip <- 0
#   for (i in 1:n){
#     query = sprintf("MATCH (n:Chemical) RETURN id(n) AS id, n.pk AS pk, n.formula AS formula ORDER BY id SKIP %s LIMIT %s",
#                     skip,
#                     limit)
#     skip <- skip + limit
#   }
# }














# 2. here to populate the data into the database
# MATCH (n:Chemical) SET n.pk = 1