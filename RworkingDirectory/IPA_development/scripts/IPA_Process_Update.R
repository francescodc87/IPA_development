rm(list = ls())
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/UpdateRelated/IPA_Update.cpp')
library("Matrix")

post <- matrix(0, 6, 10)
post[1,1] <- 0.92   
post[2,2] <- 0.94 
post[3,2] <- 0.95 
post[4,5] <- 0.93
post[5,8] <- 0.95
post[6,9] <- 0.96
post_spMat <- as(post, "dgCMatrix")
iso <- matrix(0, 10, 10)
iso[1,2] <- 5
iso[1,8] <- 2
iso[1,9] <- 8
iso[2,8] <- 0.2
iso[2,1] <- 0.2
iso[2,9] <- 2
iso[8,1] <- 0.5
iso[8,2] <- 5
iso[8,9] <- 8
iso[9,1] <- 0.125
iso[9,2] <- 0.5
iso[9,8] <- 0.125
iso_spMat <- as(iso, "dgCMatrix")

mass <- c(40.0001,99.99,100.0001,50.0001,200.0002,400.0008)
compMass <- c(40,100,30,35,50,60,70,200,400,250)
massInt <- c(202,100,50.5,101,252.5,25.25)
pk <- rep(0.1,10)
massMZ <- c(40.0001,99.99,100.0001,50.0001,200.0002,400.0008)
massRT <- c(50,75,100,150,200,250)
obsMA <- vector('numeric')
recMA <- vector('numeric')
recMAdistr <- NULL
obsIR <- NULL
obsIRlComp <- vector('numeric')
obsIRrComp <- vector('numeric')
recIRdistr <- NULL
recIRdistrLComp <- vector('numeric')
recIRdistrRComp <- vector('numeric')
#                 NumericVector monoAddComp,              // compound number of all monoisotopic adducts
#                 List otherAdd,                          // compound number of all other adducts except monoisotopic adducts
#                 List recMostIntAddInfo,                 // observed most intense adduct information of each monoisotopic adduct
#                 NumericVector recMainAdd,               // observed main adduct information of each monoisotopic adduct
monoAddComp <- c(1,4,7,8)
allAdd <- list(c(0,1,2,3,4),c(0,1,2,3,4),c(5,6,7,8,9),c(5,6,7,8,9))
recMostIntAddInfo <- list(list("compNum" = 1, "recTime" = 1),
                          list("compNum" = 1, "recTime" = 1),
                          list("compNum" = 7, "recTime" = 1),
                          list("compNum" = 7, "recTime" = 1))
recMainAdd <- c(1,1,7,7)
obsRT <- NULL
recRT <- NULL
recRTdistr <- NULL
obsRTcomp <- vector('numeric')
recRTdistrComp <- vector('numeric')
ma <- 3
s_pk <- 0.1
s_ma <- 0.1
s_iso <- 0.1
t_mainAdd <- 2
threshold <- 0.9


l <- UpdateMain(post_spMat, iso_spMat, mass, compMass, massInt, pk, massMZ, massRT, obsMA, recMA, recMAdistr,
                obsIR, obsIRlComp, obsIRrComp, recIRdistr, recIRdistrLComp, recIRdistrRComp,
                monoAddComp, allAdd, recMostIntAddInfo, recMainAdd,
                obsRT, recRT, recRTdistr, obsRTcomp, recRTdistrComp, ma, s_pk, s_ma, s_iso, t_mainAdd, threshold)


mass_new <- c(40.0001,99.99,100.0002,50.00005,200.0004,400.0008)
compMass_new <- c(40,100,30,35,50,60,70,200,400,250)
massInt_new <- c(252.5,100,25.25,101,202,50.5)
massRT_new <- c(40,75,80,200,150,200)
l_new <- UpdateMain(post_spMat, iso_spMat, mass_new, compMass_new, massInt_new, l$model$pk, massMZ, massRT_new, l$other$ma$obsMA, l$record$ma$recMA, l$record$ma$recMAdistr,
                    l$other$ir$obsIR, l$other$ir$obsIRlComp, l$other$ir$obsIRrComp, l$record$ir$recIRdistr, l$record$ir$recIRdistrLComp, l$record$ir$recIRdistrRComp,
                    monoAddComp, allAdd, l$record$add$recMostIntAddInfo, l$record$add$recMainAdd,
                    l$other$rt$obsRT, l$record$rt$recRT, l$record$rt$recRTdistr, l$other$rt$obsRTcomp, l$record$rt$recRTdistrComp, l$model$ma, s_pk, s_ma, s_iso, t_mainAdd, threshold)


# // [[Rcpp::export]]
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
#                 List obsIR,                             // storage of observed intensity ratios of all detected isotope pairs in this update time
#                 NumericVector obsIRlComp,               // left compound number of stored observed intensity ratio
#                 NumericVector obsIRrComp,               // right compound number of stored observed intensity ratio
#                 List recIRdistr,                        // record all distributions of all isotope pairs
#                 NumericVector recIRdistrLComp,          // left compound number of all recorded distributions of all isotope pairs
#                 NumericVector recIRdistrRComp,          // right compound number of all recorded distributions of all isotope pairs
#                 NumericVector monoAddComp,              // compound number of all monoisotopic adducts
#                 List otherAdd,                          // compound number of all other adducts except monoisotopic adducts
#                 List recMostIntAddInfo,                 // observed most intense adduct information of each monoisotopic adduct
#                 NumericVector recMainAdd,               // observed main adduct information of each monoisotopic adduct
#                 List obsRT,                             // storage of observed retention time in this update time
#                 List recRT,                             // record all observed retention time
#                 List recRTdistr,                        // record all the distributions of retention time
#                 NumericVector obsRTcomp,                // the compound numbers of the observed retention time in this update time
#                 NumericVector recRTdistrComp,           // the compound numbers of recorded retention time distributions
#                 double ma,                              // mass accuracy
#                 double s_pk,                            // update scale(speed) for pk
#                 double s_ma,                            // update scale for ma
#                 double s_iso,                           // update scale for iso
#                 int t_mainAdd,                          // threshold observed time value for main adduct shuffle
#                 double threshold                        // threshold used to filter posterior
#                 
# )


# k <- 2
# l <- NULL
# 
# massInt_new <- c(50,25,75,125)
# massMZ_new <- c(20.01,39.99,30,40)
# massRT_new <- c(30,40,50,60)
# for(i in 1:k){
#   if(is.null(l)){
#     l <- UpdateMain(post_spMat, iso_spMat, mass, compMass, massInt, pk, massMZ, massRT, obsMA, recMA, recMAdistr,
#                     obsIR, obsIRlComp, obsIRrComp, recIRdistr, recIRdistrLComp, recIRdistrRComp,
#                     monoAddComp, otherAdd, recMostIntAddInfo, recMainAdd,
#                     obsRT, recRT, recRTdistr, obsRTcomp, recRTdistrComp, ma, s_pk, s_ma, s_iso, t_mainAdd, threshold)
#   }else{
#     l_new <- UpdateMain(post_spMat, l$iso, mass, compMass, massInt_new, l$pk, massMZ_new, massRT_new, l$obsMA, l$recMA, l$recMAdistr,
#                         l$obsIR, l$obsIRlComp, l$obsIRrComp, l$recIRdistr, l$recIRdistrLComp, l$recIRdistrRComp, 
#                         l$obsRT, l$recRT, l$recRTdistr, l$obsRTcomp, l$recRTdistrComp, l$ma, s_pk, s_ma, s_iso, threshold)
#   }
#   
# }


# post <-matrix( c(0.03,0.03,0.92,0.02,0.90,0.04,0.03,0.03,0.92,0.04,0.02,0.02,0.01,0.01,0.97,0.01), nrow=4, ncol=4, byrow = TRUE)        # fill matrix by rows 
# post_spMat <- as(post, "dgCMatrix")
# iso <- matrix( c(0,0,0.4,0,0,0,0,0.2,2.5,0,0,0,0,5,0,0), nrow=4, ncol=4, byrow = TRUE)
# iso_spMat <- as(iso, "dgCMatrix")
# mass <- c(20.01,39.99,60,99.99)
# compMass <- c(40,20,100,120)
# massInt <- c(50,25,75,100)
# pk <- c(0.6,0.2,0.1,0.1)
# massMZ <- c(20.01,39.99,30,40)
# massRT <- c(40,65,85,75)
# obsMA <- vector('numeric')
# recMA <- vector('numeric')
# recMAdistr <- NULL
# obsIR <- NULL
# obsIRlComp <- vector('numeric')
# obsIRrComp <- vector('numeric')
# recIRdistr <- NULL
# recIRdistrLComp <- vector('numeric')
# recIRdistrRComp <- vector('numeric')
# obsRT <- NULL
# recRT <- NULL
# recRTdistr <- NULL
# obsRTcomp <- vector('numeric')
# recRTdistrComp <- vector('numeric')
# ma <- 3
# s_pk <- 0.1
# s_ma <- 0.1
# s_iso <- 0.1
# threshold <- 0.9
# 
# 
# k <- 2
# l <- NULL
# 
# massInt_new <- c(50,25,75,125)
# massMZ_new <- c(20.01,39.99,30,40)
# massRT_new <- c(30,40,50,60)
# for(i in 1:k){
#   if(is.null(l)){
#     l <- UpdateMain(post_spMat, iso_spMat, mass, compMass, massInt, pk, massMZ, massRT, obsMA, recMA, recMAdistr,
#                     obsIR, obsIRlComp, obsIRrComp, recIRdistr, recIRdistrLComp, recIRdistrRComp, 
#                     obsRT, recRT, recRTdistr, obsRTcomp, recRTdistrComp, ma, s_pk, s_ma, s_iso, threshold)
#   }else{
#     l <- UpdateMain(post_spMat, l$iso, mass, compMass, massInt_new, l$pk, massMZ_new, massRT_new, l$obsMA, l$recMA, l$recMAdistr,
#                     l$obsIR, l$obsIRlComp, l$obsIRrComp, l$recIRdistr, l$recIRdistrLComp, l$recIRdistrRComp, 
#                     l$obsRT, l$recRT, l$recRTdistr, l$obsRTcomp, l$recRTdistrComp, l$ma, s_pk, s_ma, s_iso, threshold)
#   }
#   
# }
