rm(list=ls())
library(rbenchmark)
library(Matrix)
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_RT.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_relID.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_corr.R')

source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Multsample_R.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Compute_post.R')

source('~/RworkingDirectory/IPA_development/class/Other/IPA_Posterior_Efficiency_BeforeGs.R')
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/Other/IPA_Posterior_Efficiency_GS.cpp')   #'GS' = 'GibbsSampling'

source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_ADD_ISO_BIO_INT_NoPot.R')

load("~/RworkingDirectory/IPA_development/data/dataset_prova.Rdata")

Iso <- matrix(sample(0:100,18*18, replace=TRUE),18,18)
Iso[lower.tri(Iso)] <- t(Iso)[lower.tri(Iso)]
idx <- which(Iso < 50)
Iso[idx] <- 0
Int = c(65,70,75,80,85,88,90,92,95)

Prior_spMat <- as(Prior, "dgCMatrix")
Add_spMat <- as(adducts, "dgCMatrix")
Iso_spMat <- as(Iso, "dgCMatrix")
Bio_spMat <- as(biotransforamtions, "dgCMatrix")

# postR <- ComputePosteriorR_Add_Iso_Bio_Int_NoPot(P = Prior, Add = adducts, Iso = Iso, Int = Int, Bio = biotransforamtions, RT = RTs, relId = rel.id, corrMat = Corr.matrix)
# postRcpp <- Efficiency_PosteriorRcpp_Add_Iso_Bio_Int_NoPot(P = Prior_spMat, Add = Add_spMat, Iso = Iso_spMat, Int = Int, Bio = Bio_spMat, RT = RTs, relId = rel.id, corrMat = Corr.matrix)

 
# # RUNNING TIME COMPARE
# benchmark(replications = 10,
#           postR = ComputePosteriorR_Add_Iso_Bio_Int_NoPot(P = Prior, Add = adducts, Iso = Iso, Int = Int, Bio = biotransforamtions, RT = RTs, relId = rel.id, corrMat = Corr.matrix),
#           postRcpp = Efficiency_PosteriorRcpp_Add_Iso_Bio_Int_NoPot(P = Prior_spMat, Add = Add_spMat, Iso = Iso_spMat, Int = Int, Bio = Bio_spMat, RT = RTs, relId = rel.id, corrMat = Corr.matrix),
#           columns = c('test','elapsed','relative','replications'),
#           order = c('relative'),
#           relative = 'elapsed')

# RESULT COMPARE
tmp <- NULL
for(i in 1:100){
postR <- ComputePosteriorR_Add_Iso_Bio_Int_NoPot(P = Prior, Add = adducts, Iso = Iso, Int = Int, Bio = biotransforamtions, RT = RTs, relId = rel.id, corrMat = Corr.matrix)
postRcpp <- Efficiency_PosteriorRcpp_Add_Iso_Bio_Int_NoPot(P = Prior_spMat, Add = Add_spMat, Iso = Iso_spMat, Int = Int, Bio = Bio_spMat, RT = RTs, relId = rel.id, corrMat = Corr.matrix)
tmp <- rbind(tmp, as.vector(t(postR)-t(postRcpp)))
}
boxplot(tmp)







# PostRcpp = Efficiency_PosteriorRcpp_Add_Iso_Bio_Int_NoPot(P=Prior,Add = Add,
#                                                           Iso = Iso, RT = RTs,
#                                                           Bio = Bio, Int = Int,
#                                                           RTwin = 5, it =1, burn = 0,
#                                                           allSamp = F, delAdd =0.4,
#                                                           delIso = 0.2, delBio =1, v = T)


# Add <- matrix(sample(0:100,10000*10000, replace=TRUE),10000,10000)
# Iso <- matrix(sample(0:1,10000*10000, replace=TRUE),10000,10000)
# Bio <- matrix(sample(0:1,10000*10000, replace=TRUE),10000,10000)
# Add[lower.tri(Add)] <- t(Add)[lower.tri(Add)]
# Iso[lower.tri(Iso)] <- t(Iso)[lower.tri(Iso)]
# Bio[lower.tri(Bio)] <- t(Bio)[lower.tri(Bio)]
# m <- matrix(sample(0:1000,1000*10000, replace=TRUE),1000,10000)
# Prior <- m/rowSums(m)[row(m)]
# RT <- c(sample(40:160,1000,replace=TRUE))
# Int <- c(sample(40:1000,1000,replace=TRUE))
# 
# # Prior[which(is.na(Prior))] <- 0
# 
# samplerRcpp = ComputePosteriorRcpp_Add_Iso_Bio_Int_NoPot(P=Prior,Add = Add,
#                                            Iso = Iso, RT = RT,
#                                            Bio = Bio, Int = Int,
#                                            RTwin = 5, it =1, burn = 0,
#                                            allSamp = F, delAdd =0.4,
#                                            delIso = 0.2, delBio =1, v = T)