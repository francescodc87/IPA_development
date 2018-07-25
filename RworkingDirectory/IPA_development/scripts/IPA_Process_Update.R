<<<<<<< HEAD
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
=======
rm(list = ls())
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/UpdateRelated/IPA_Update.cpp')
library("Matrix")

fomuCompId <- c(0,0,0,0,0,1,1,1,1,1)

post <- matrix(0, 6, 10)
post[1,1] <- 0.92   
post[2,2] <- 0.94 
post[3,2] <- 0.95 
post[4,5] <- 0.93
post[5,8] <- 0.95
post[6,9] <- 0.96
post[2,4] <- 0.9
post[3,4] <- 0.9
post[4,7] <- 0.9
post[5,7] <- 0.9
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

add <- matrix(0, 10, 10)
add[2,2] <- 1
add[3,2] <- 1
add[5,2] <- 1
add[6,8] <- 1
add[8,8] <- 1
add[9,8] <- 1
add_spMat <- as(add, "dgCMatrix")

bio <- matrix(0, 10, 10)
bio[2,8] <- 1
bio[8,2] <- 1
bio_spMat <- as(bio, "dgCMatrix")

pkFomu <- rep(0.1,10)
pkComp <- rep(0.5,2)

ma <- 3
massMZ <- c(40.0001,99.99,100.0001,49.9999,200.0002,400.0008)
fomuMZ <- c(40,100,30,35,50,60,70,200,400,250)
obsMA <- vector('numeric')
recMA <- vector('numeric')
recMAdistr <- NULL

massInt <- c(202,100,50.5,101,252.5,25.25)
obsIR <- NULL
obsIRlFomu <- vector('numeric')
obsIRrFomu <- vector('numeric')
recIRdistr <- NULL
recIRdistrLFomu <- vector('numeric')
recIRdistrRFomu <- vector('numeric')

recMonoFomu <- c(1,2,4,5,7,8) 
recMonoMainAddFomu <- c(1,1,1,7,7,7)
recMonoFomuWithInSameComp <- list(c(1,2,4),c(1,2,4),c(1,2,4),c(5,7,8),c(5,7,8),c(5,7,8))
recMITwithInSameComp <- list(c(1,0,0),c(1,0,0),c(1,0,0),c(0,1,0),c(0,1,0),c(0,1,0))
recCompMainAddFomu <- c(1,7)
recCompMonoFomu <- list(c(1,2,4),c(5,7,8))
recCompBioLink <- list(c(1),c(0))

massRT <- c(50,75,100,150,200,250)
obsFomuRT <- NULL
obsRTfomu <- vector('numeric')
recFomuRT <- NULL
recFomuRTdistr <- NULL
recRTdistrFomu <- vector('numeric')

s_pkAdd <- 0.1
s_pkComp <- 0.1
s_ma <- 0.1
s_iso <- 0.1
threshold <- 0.9




l <- UpdateMain(post_spMat, iso_spMat, add_spMat, bio_spMat, fomuCompId, pkFomu, pkComp, 
                ma, massMZ, fomuMZ, obsMA, recMA, recMAdistr, 
                massInt, obsIR, obsIRlFomu, obsIRrFomu, recIRdistr, recIRdistrLFomu, recIRdistrRFomu,
                recMonoFomu, recMonoMainAddFomu, recMonoFomuWithInSameComp, recMITwithInSameComp, recCompMainAddFomu, recCompMonoFomu, recCompBioLink,
                massRT, obsFomuRT, obsRTfomu, recFomuRT, recFomuRTdistr, recRTdistrFomu,
                s_pkAdd, s_pkComp, s_ma, s_iso, threshold)
 
massMZ_new <- c(39.9999,99.99,100.0002,50.00005,200.0004,400.0008)
fomuMZ_new <- c(40,100,30,35,50,60,70,200,400,250)
massInt_new <- c(252.5,100,25.25,101,202,50.5)
massRT_new <- c(40,75,80,200,150,200)
l_new <- UpdateMain(post_spMat, l$model$iso, l$model$add, l$model$bio, fomuCompId, l$model$pkFomu, l$model$pkComp, 
                    l$model$ma, massMZ_new, fomuMZ_new, l$other$ma$obsMA, l$record$ma$recMA, l$record$ma$recMAdistr, 
                    massInt_new, l$other$iso$obsIR, l$other$iso$obsIRlFomu,l$other$iso$obsIRrFomu, l$record$iso$recIRdistr, l$record$iso$recIRdistrLFomu, l$record$iso$recIRdistrRFomu,
                    l$record$add$recMonoFomu, l$record$add$recMonoMainAddFomu, l$record$add$recMonoFomuWithInSameComp, l$record$add$recMITwithInSameComp, l$record$bio$recCompMainAddFomu, l$record$bio$recCompMonoFomu, l$record$bio$recCompBioLink,
                    massRT_new, l$other$rt$obsFomuRT, l$other$rt$obsRTfomu, l$record$rt$recFomuRT, l$record$rt$recFomuRTdistr, l$record$rt$recRTdistrFomu,
                    s_pkAdd, s_pkComp, s_ma, s_iso, threshold)


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
#                 List recFomuRTdistr,                    // record all the distributions of retention time values of fomulas
#                 NumericVector recRTdistrFomu,           // the fomula numbers of the recorded retention time distributions
#                 double s_pkAdd,                         // update scale(speed) for pkFomu
#                 double s_pkComp,                        // update scale(speed) for pkComp
#                 double s_ma,                            // update scale for ma (0< <=1)
#                 double s_iso,                           // update scale for iso (0< <=1)
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
>>>>>>> a94fb0c2a142d258bb253d79559f22b40ae89587
