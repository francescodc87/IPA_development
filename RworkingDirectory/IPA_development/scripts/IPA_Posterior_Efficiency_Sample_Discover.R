library(rbenchmark)
rm(list=ls())




source('~/RworkingDirectory/IPA_development/class/PosteriorRelated/IPA_Posterior_Multsample_R.R')
Rcpp::sourceCpp('~/RworkingDirectory/IPA_development/class/Other/IPA_Posterior_Efficiency_sampleRcpp.cpp')


m <- 10000
c <- 10000000

X <- sample(0:10,1*c, replace=TRUE)
seqMass <- seq(m)
seqComp <- seq(c)

# cat("m: ", m, "; c: ", c)

benchmark(replications = 1,
          sampleR = sample(m),
          sampleRcpp = sampleRcppExport(seqMass, m, FALSE),
          columns = c('test','elapsed','relative','replications'),
          order = c('relative'),
          relative = 'elapsed')


# efficiency comparison of multsample & sampleRcppExport
benchmark(replications = 1,
          multSampleR = multsample(X/sum(X)),
          multSampleRcpp = sampleRcppExport(seqComp, 1, TRUE, X/sum(X)),
          columns = c('test','elapsed','relative','replications'),
          order = c('relative'),
          relative = 'elapsed')





# efficiency comparison of sample & sampleRcppExport
