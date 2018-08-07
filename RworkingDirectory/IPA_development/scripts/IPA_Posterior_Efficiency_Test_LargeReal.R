rm(list=ls())
library(Matrix)
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_RT.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_relID.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Checking_corr.R')

source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Multsample_R.R')
source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Compute_post.R')

source('~/RworkingDirectory/IPA_development/class/Other/IPA_Posterior_Efficiency_BeforeGs.R')
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/Other/IPA_Posterior_Efficiency_GS.cpp')   #'GS' = 'GibbsSampling'


load("~/RworkingDirectory/IPA_development/data/PosteriorEfficiencyTest_Rcpp_180718.Rdata")

PostRcpp = Efficiency_PosteriorRcpp_Add_Iso_Bio_Int_NoPot(P=Prior_filtered,Add = Add,
                                                          Iso = Iso, RT = RT,
                                                          Bio = Bio, Int = Int,
                                                          RTwin = 5, it =1, burn = 0,
                                                          allSamp = F, delAdd =0.4,
                                                          delIso = 0.2, delBio =1, v = T)








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