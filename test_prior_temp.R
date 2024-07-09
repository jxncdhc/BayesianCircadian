simdata1 <- rjmcmc_getSimData(n=24)
yy <- simdata1$Y
tt <- simdata1$time

bayes_out <- bayes_rhythmicity(yy=yy,tt=tt,M=400)
qvalue <- bayes_out$qvalues
prior_info <- ifelse(qvalue<0.05,1,0)

#set.seed(i*27)
simdata2 <- rjmcmc_getSimData(n=24)
yy <- simdata2$Y
tt <- simdata2$time

bayes_rhythmicity_prior(yy=yy,tt=tt,M=400)
