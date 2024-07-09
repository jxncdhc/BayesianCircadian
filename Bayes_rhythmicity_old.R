##' Bayesian Method for Circadian Rhythmicity Detection. 
##'
##' Test the significance of circadian curve fitting using Bayesian method.
##' @title Bayesian Method for Detecting Circadian Rhythmicity.
##' @param tt Time vector.
##' @param yy Gene expression value data matrix. Rows are genes' expression values, columns are observation time points. 
##' @param period Period of your data. Default is 24.
##' @param M number of MCMC iterations. Default is 3000.
##' @param alpha00 hyper parameter. Default is 0. 
##' @param beta10 hyper parameter. Default is 0. 
##' @param beta20 hyper parameter. Default is 0. 
##' @param tau2a hyper parameter. Default is 10.
##' @param tau2b1 hyper parameter. Default is 10.
##' @param tau2b2 hyper parameter. Default is 10.
##' @param mu0 hyper parameter. Default is 0. 
##' @param nu0 hyper parameter. Default is 0. 
##' @param a_p hyper parameter. Default is 1. 
##' @param b_p hyper parameter. Default is 1. 
##' 
##' @return A list of Bayesian q-values and running time in seconds. 
##' \item{qvalues}{A vector of q-values for each gene.}
##' \item{alpha0}{A matrix contains all non-circadian genes alpha values for each iteration.}
##' \item{sigma20}{A matrix contains all non-circadian genes sigma square values for each iteration.}
##' \item{alpha1}{A matrix contains all circadian genes alpha values for each iteration.}
##' \item{beta1}{A matrix contains all circadian genes beta1 values for each iteration.}
##' \item{beta2}{A matrix contains all circadian genes beta2 values for each iteration.}
##' \item{sigma21}{A matrix contains all circadian genes sigma square values for each iteration.}
##' \item{r}{Gene indicator. 0 is non-circadian gene, 1 is circadian gene.}
##' \item{k}{Sum of circadian genes in each MCMC iteration.}
##' \item{p}{Posterior probability of number of circadian gene.}
##' \item{timeUsed}{Running time used in seconds.}
##' @author Zhiguang Huo, Haocheng Ding, Lingsong Meng
##' @export
##' @examples
##' set.seed(32611)
bayes_rhythmicity <- function(yy, tt, period=24, M=3000, alpha00 = 0, beta10 = 0, beta20 = 0, tau2a = 10, tau2b1 = 10, tau2b2 = 10, mu0 = 0, nu0 = 0, a_p = 1, b_p = 1){  
  
  w <- 2*pi/period
  n <- length(tt)
  G <- nrow(yy)
  
  # Prior hyper-parameters
  alpha00 <- alpha00
  beta10 <- beta10
  beta20 <- beta20
  tau2a <- tau2a
  tau2b1 <- tau2b1
  tau2b2 <- tau2b2
  mu0 <- mu0
  nu0 <- nu0
  a_p <- a_p
  b_p <- b_p
  
  # Setting up starting values
  alpha0 <- alpha1 <- beta1 <- beta2 <- beta01 <- beta02 <- rep(0,G); sig20 <- sig21 <- rep(.5,G)
  
  # Gibbs sampler
  p <- 0.5
  r <- rbinom(G,1,p)
  k <- sum(r)
  draw_alpha0 <- draw_alpha1 <- draw_beta1 <- draw_beta2 <- draw_sig20 <- draw_sig21 <- draw_r <- matrix(0,M,G)
  draw_p <- rep(0,M)
  draw_k <- rep(0,M)
  
  # give initial parameter value 
  draw_alpha1[1,] <- alpha1
  draw_alpha0[1,] <- alpha0
  draw_beta1[1,] <- beta1
  draw_beta2[1,] <- beta2
  draw_sig21[1,] <- sig21
  draw_sig20[1,] <- sig20
  ##############
  
  y <- yy
  x1 <- cos(w*tt)
  x2 <- sin(w*tt)
  xMat <- rbind(rep(1,n),x1,x2)
  betaMat1 <- rbind(alpha1,beta1,beta2)
  betaMat0 <- rbind(alpha0,beta01,beta02)
  
  ######### MCMC iteration ##########
  
  start_time <- Sys.time()
  tY <- t(y)
  
  
  for(i in 2:M){
    # Step 1
    k <- sum(r)
    p <- rbeta(1, k+a_p, G-k+b_p)
    #########
    
    # Step 2
    #sum (y-a-b1-b2)^2
    temp_r1 <- colSums((tY-crossprod(xMat,betaMat1))^2)
    #sum (y-a)^2
    temp_r0 <- colSums((tY-crossprod(xMat,betaMat0))^2)
    
    log_A <- log(p)-1/2*n*log(2*pi*sig21) - 1/(2*sig21) * temp_r1
    log_B <- log(1-p)-1/2*n*log(2*pi*sig20) - 1/(2*sig20) * temp_r0
    
    
    pr <- 1/(1+exp(log_B-log_A))
    
    r <- rbinom(G,1,pr)
    ##################
    
    # Step 3 update beta1 ###
    var_beta1 = 1/(1/tau2b1 + c(tcrossprod(xMat[2,,drop=F],xMat[2,,drop=F]))/sig21)
    mean_beta1 = var_beta1*(xMat[2,,drop=F]%*%(tY-crossprod(xMat[c(1,3),],betaMat1[c(1,3),]))/sig21+beta10/tau2b1)
    beta1 <- rnorm(G,mean_beta1,sqrt(var_beta1))
    betaMat1['beta1',] <- beta1
    ####################
    
    # Step 4 update beta2 ###
    var_beta2 = 1/(1/tau2b2 + c(tcrossprod(xMat[3,,drop=F],xMat[3,,drop=F]))/sig21)
    mean_beta2 = var_beta2*(xMat[3,,drop=F]%*%(tY-crossprod(xMat[c(1,2),],betaMat1[c(1,2),]))/sig21+beta20/tau2b2)
    beta2 <- rnorm(G,mean_beta2,sqrt(var_beta2))
    betaMat1['beta2',] <- beta2
    ####################
    
    ### Step 5 update alpha r=0 ###
    var_alpha0 <- 1/(1/tau2a + n/sig20)
    mean_alpha0 <- var_alpha0*(colSums(tY)/sig20 + alpha0/tau2a)
    alpha0 <- rnorm(G,mean_alpha0,sqrt(var_alpha0))
    betaMat0['alpha0',] <- alpha0
    ####################
    
    # Step 6 update alpha r=1
    var_alpha1 <- 1/(1/tau2a + n/sig21)
    mean_alpha1 <- var_alpha1*(colSums(tY-crossprod(xMat[-1,],betaMat1[-1,]))/sig21 + alpha00/tau2a)
    alpha1 <- rnorm(G,mean_alpha1,sqrt(var_alpha1))
    betaMat1['alpha1',] <- alpha1
    ##################
    
    # Step 7 update sigma2 r=0 ###
    ##
    #sum(y-a)^2
    temp0 <- colSums((tY-crossprod(xMat[1,,drop=F],betaMat0[1,,drop=F]))^2)
    ##
    post_shape0 <-  nu0+n/2
    post_rate0 <- mu0+temp0*1/2
    sig20 <- 1/rgamma(G, post_shape0, post_rate0)
    ####################
    
    # Step 8 update sigma2 r=1 ###
    ##
    #sum(y-a-b1-b2)^2
    temp1 <- colSums((tY-crossprod(xMat,betaMat1))^2)
    post_shape1 <- nu0+n/2
    post_rate1 <- mu0+temp1*1/2
    ##
    sig21 <- 1/rgamma(G, post_shape1, post_rate1)  
    ####################
    
    # Step 9 acceptance probability
    
    # p(Yg|theta_Lg', Lg')/p(Yg|theta_Lg, Lg)
    Y1sampleMean <- crossprod(xMat,betaMat1)
    sig21_temp <- rep(sig21,each=n)
    Y1sample <- dnorm(tY,Y1sampleMean,sig21_temp,log=T)
    Y1sampleMat <- matrix(Y1sample,n,G)
    pYg_upper <- apply(Y1sampleMat, 2, sum)
    
    Y0sampleMean <- crossprod(xMat,betaMat0)
    sig20_temp <- rep(sig20,each=n)
    Y0sample <- dnorm(tY,Y0sampleMean,sig20_temp, log=T)
    Y0sampleMat <- matrix(Y0sample,n,G)
    pYg_lower <- apply(Y0sampleMat, 2, sum)
    
    pYg_ratio <- pYg_upper-pYg_lower
    
    #p(theta_Lg')/p(theta_Lg)
    ptheta_upper <- dnorm(betaMat1['alpha1',],alpha00,sqrt(tau2a),log=T)+dnorm(betaMat1['beta1',],beta10,sqrt(tau2b1),log=T)+dnorm(betaMat1['beta2',],beta20,sqrt(tau2b2),log=T)-log(sig21)
    ptheta_lower <- dnorm(betaMat0['alpha0',],alpha00,sqrt(tau2a),log=T)-log(sig20)
    ptheta_ratio <- ptheta_upper-ptheta_lower
    
    # p(Lg')/p(Lg) = p/(1-p)
    pLg_ratio <- log(p)-log(1-p)
  
    
    # J(Lg' to Lg)/J(Lg to Lg')
    prob_logA <- log(p)+pYg_upper
    prob_logB <- log(1-p)+pYg_lower
    J_ratio <- prob_logB-prob_logA

    
    # q(uLg')/q(uLg)
    q_uLg_upper <- dnorm(betaMat0['alpha0',],mean_alpha0,sqrt(var_alpha0),log=T)+dinvgamma(sig20, post_shape0, post_rate0,log=T)
    q_uLg_lower <- dnorm(betaMat1['alpha1',],mean_alpha1,sqrt(var_alpha1),log=T)+dnorm(betaMat1['beta1',],mean_beta1,sqrt(var_beta1),log=T)+dnorm(betaMat1['beta2',],mean_beta2,sqrt(var_beta2),log=T)+dinvgamma(sig21,post_shape1,post_rate1,log=T)
    q_uLg_ratio <- q_uLg_upper-q_uLg_lower
    
    
    # accept prob
    log_accept_prob_0to1_ratio <- pYg_ratio+ptheta_ratio+pLg_ratio+J_ratio+q_uLg_ratio
    log_accept_prob_0to1 <- pmin(log_accept_prob_0to1_ratio,0)
    
    log_accept_prob_1to0_ratio <- -log_accept_prob_0to1_ratio
    log_accept_prob_1to0 <- pmin(log_accept_prob_1to0_ratio,0)
    
    model_change <- r-draw_r[(i-1),]
    
    # 0 to 1 
    index_01 <- which(model_change == 1)
    index_10 <- which(model_change == -1)
    
    length_01 <- length(log_accept_prob_0to1[index_01])
    change_index_01 <- index_01[log_accept_prob_0to1[index_01]<=log(runif(length_01))]
    r[change_index_01] <- 0
    
    alpha1[change_index_01] <- draw_alpha1[(i-1),change_index_01]
    alpha0[change_index_01] <- draw_alpha0[(i-1),change_index_01]
    beta1[change_index_01] <- draw_beta1[(i-1),change_index_01] 
    beta2[change_index_01] <- draw_beta2[(i-1),change_index_01]
    sig21[change_index_01] <- draw_sig21[(i-1),change_index_01]
    sig20[change_index_01] <- draw_sig20[(i-1),change_index_01]
    
    
    length_10 <- length(log_accept_prob_1to0[index_10])
    change_index_10 <- index_10[log_accept_prob_1to0[index_10]<=log(runif(length_10))]
    r[change_index_10] <- 1
    
    
    alpha1[change_index_10] <- draw_alpha1[(i-1),change_index_10]
    alpha0[change_index_10] <- draw_alpha0[(i-1),change_index_10]
    beta1[change_index_10] <- draw_beta1[(i-1),change_index_10] 
    beta2[change_index_10] <- draw_beta2[(i-1),change_index_10]
    sig21[change_index_10] <- draw_sig21[(i-1),change_index_10]
    sig20[change_index_10] <- draw_sig20[(i-1),change_index_10]
    
    ###################
   
    draw_alpha1[i,] <- alpha1
    draw_alpha0[i,] <- alpha0
    draw_beta1[i,] <- beta1
    draw_beta2[i,] <- beta2
    draw_sig21[i,] <- sig21
    draw_sig20[i,] <- sig20
    draw_r[i,] <- r
    draw_k[i] <- k
    draw_p[i] <- p
  }
  end_time <- Sys.time()
  timeUsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  #################################
  
  
  Ding_p_L1 <- colSums(draw_r[(M/2+1):M,])/(M/2)
  
  p_g <- 1-Ding_p_L1
  BayesianFDR_qvalues <- numeric(0)
  for(i in 1:length(p_g)){
    BayesianFDR_qvalues[i] <- sum((p_g<=p_g[i])*p_g)/sum(p_g<=p_g[i])
  }

  #########################
  return(list(qvalues = BayesianFDR_qvalues, alpha0 = draw_alpha0, sigma20 = draw_sig20, alpha1 = draw_alpha1, beta1 = draw_beta1, beta2 = draw_beta2, sigma21 = draw_sig21, r = draw_r, k = draw_k, p = draw_p , timeUsed=timeUsed))
}