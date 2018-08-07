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
add[3,2] <- 1
add[2,3] <- 1
add[5,2] <- 1
add[2,5] <- 1
add[6,8] <- 1
add[8,6] <- 1
add[9,8] <- 1
add[8,9] <- 1
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
recRTfomu <- vector('numeric')
recFomuRTdistr <- NULL
recRTdistrFomu <- vector('numeric')

s_pkAdd <- 0.1
s_pkComp <- 0.1
s_ma <- 0.1
s_iso <- 0.1
w_int <- 0.1
w_ma <- 0.1
w_addLink <- 0.4
w_bioLink <- 0.4
t_RT <- 49
t_post <- 0.9




l <- UpdateMain(post_spMat, iso_spMat, add_spMat, bio_spMat, fomuCompId, pkFomu, pkComp, 
                ma, massMZ, fomuMZ, obsMA, recMA, recMAdistr, 
                massInt, obsIR, obsIRlFomu, obsIRrFomu, recIRdistr, recIRdistrLFomu, recIRdistrRFomu,
                recMonoFomu, recMonoMainAddFomu, recMonoFomuWithInSameComp, recMITwithInSameComp, recCompMainAddFomu, recCompMonoFomu, recCompBioLink,
                massRT, obsFomuRT, obsRTfomu, recFomuRT, recRTfomu, recFomuRTdistr, recRTdistrFomu,
                s_pkAdd, s_pkComp, s_ma, s_iso, w_int, w_ma, w_addLink, w_bioLink, t_RT, t_post)

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
#                 NumericVector recRTfomu,                // the fomula numbers of the recorded retention time values
#                 List recFomuRTdistr,                    // record all the distributions of retention time values of fomulas
#                 NumericVector recRTdistrFomu,           // the fomula numbers of the recorded retention time distributions
#                 double s_pkAdd,                         // update scale(speed) for pkFomu
#                 double s_pkComp,                        // update scale(speed) for pkComp
#                 double s_ma,                            // update scale for ma (0< <=1)
#                 double s_iso,                           // update scale for iso (0< <=1)
#                 double t_RT,                            // threshold used to count iso and add connection when solving overlapped assignments automatically
#                 double t_post                           // threshold used to filter posterior
#                 
# )

