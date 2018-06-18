library(rbenchmark)
rm(list=ls())
load("~/RworkingDirectory/IPA_development/data/dataset_prova.Rdata")
source('~/RworkingDirectory/IPA_development/R/PosteriorRelated/IPA_Posterior_Checking_RT.R')
source('~/RworkingDirectory/IPA_development/R/PosteriorRelated/IPA_Posterior_Checking_relID.R')
source('~/RworkingDirectory/IPA_development/R/PosteriorRelated/IPA_Posterior_Checking_corr.R')

Rcpp::sourceCpp('R/PosteriorRelated/IPA_Posterior_Sample_Rcpp.cpp')
Rcpp::sourceCpp('R/PosteriorRelated/IPA_Posterior_GS_ADD_ISO_BIO.cpp')   #'GS' = 'GibbsSampling'
source('~/RworkingDirectory/IPA_development/R/PosteriorRelated/IPA_Posterior_Multsample_Rcpp.R')
source('~/RworkingDirectory/IPA_development/R/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_Rcpp.R')

source('~/RworkingDirectory/IPA_development/R/PosteriorRelated/IPA_Posterior_Multsample_R.R')
source('~/RworkingDirectory/IPA_development/R/PosteriorRelated/IPA_Posterior_Compute_post.R')
source('~/RworkingDirectory/IPA_development/R/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_R.R')

## These two are for debug use
# debug(IPA.sampler.Add.Iso.Bio)
# debug(IPA.sampler.Add.Iso.Bio.Rcpp)




# using dataset_prova.Rdata to test if prior computing in R and Rcpp agreed
# postR <- ComputePosteriorR_Add_Iso_Bio(P = Prior, Add = adducts, Iso = isotopes, Bio = biotransforamtions, RT = RTs, relId = rel.id, corrMat = Corr.matrix)
# postRcpp <- ComputePosteriorRcpp_Add_Iso_Bio(P = Prior, Add = adducts, Iso = isotopes, Bio = biotransforamtions, RT = RTs, relId = rel.id, corrMat = Corr.matrix)
# tmp <- NULL
# for(i in 1:100){
#   postR <- ComputePosteriorR_Add_Iso_Bio(P = Prior, Add = adducts, Iso = isotopes, Bio = biotransforamtions, RT = RTs, relId = rel.id, corrMat = Corr.matrix)
#   postRcpp <- ComputePosteriorRcpp_Add_Iso_Bio(P = Prior, Add = adducts, Iso = isotopes, Bio = biotransforamtions, RT = RTs, relId = rel.id, corrMat = Corr.matrix)
#   tmp <- rbind(tmp, as.vector(t(postR)-t(postRcpp)))
# }
# boxplot(tmp)

# # using dataset_prova.Rdata to do efficiency comparison
# benchmark(replications = 1,
#           samplerR = ComputePosteriorR_Add_Iso_Bio(P = Prior, Add = adducts, Iso = isotopes, Bio = biotransforamtions, RT = RTs, relId = rel.id, corrMat = Corr.matrix),
#           samplerRcpp = ComputePosteriorRcpp_Add_Iso_Bio(P = Prior, Add = adducts, Iso = isotopes, Bio = biotransforamtions, RT = RTs, relId = rel.id, corrMat = Corr.matrix),
#           columns = c('test','elapsed','relative','replications'),
#           order = c('relative'),
#           relative = 'elapsed')












######### following codes are not commonly used

## Here not using dataset_prova.Rdata, for test in large amount of data (because temporarily I don't have large real-world data)
# Add2 <- matrix(sample(0:1,50*50, replace=TRUE),50,50)
# Iso2 <- matrix(sample(0:1,50*50, replace=TRUE),50,50)
# Bio2 <- matrix(sample(0:1,50*50, replace=TRUE),50,50)
# Add2[lower.tri(Add2)] <- t(Add2)[lower.tri(Add2)]
# Iso2[lower.tri(Iso2)] <- t(Iso2)[lower.tri(Iso2)]
# Bio2[lower.tri(Bio2)] <- t(Bio2)[lower.tri(Bio2)]
# m <- matrix(sample(0:10,27*50, replace=TRUE),27,50)
# Prior2 <- m/rowSums(m)[row(m)]
# RTs2 <- c(RTs, RTs, RTs)
# benchmark(replications = 100,
#           samplerR = ComputePosteriorR_Add_Iso_Bio(P = Prior2, Add = Add2, Iso = Iso2, Bio = Bio2, RT = RTs2),
#           samplerRcpp = ComputePosteriorRcpp_Add_Iso_Bio(P = Prior2, Add = Add2, Iso = Iso2, Bio = Bio2, RT = RTs2),
#           columns = c('test','elapsed','relative','replications'),
#           order = c('relative'),
#           relative = 'elapsed')

## checking whether function multsample and function multsampleRcpp agreed
# po <- rep(5,100)
# po[1] <- 500
# po[6] <- 500
# po <- po/sum(po)
# tmp.R <- NULL
# tmp.Rcpp <- NULL
# for(i in 1:20000){
#   tmp.R <- c(tmp.R, multsample(po))
#   tmp.Rcpp <- c(tmp.Rcpp, multsampleRcpp(po))
# }
# boxplot(cbind(tmp.Rcpp, tmp.R))

## checking whether function multsample and function sampleRcppExport agreed (sampleRcppExport is in Gibbs_Sampling_Main.cpp, need manual add // [[Rcpp::export]] if want to test)
# po <- rep(5,100)
# po[1] <- 500
# po <- po/sum(po)
# 
# tmp.R <- NULL
# tmp.Rcpp <- NULL
# for(i in 1:50000){
#   tmp.R <- c(tmp.R, multsample(po))
#   tmp.Rcpp <- c(tmp.Rcpp, sampleRcppExport(1:length(po),1,TRUE,po))
# }
# boxplot(cbind(tmp.Rcpp, tmp.R))














