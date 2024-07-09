##' Calculate the false discovery rate (FDR) using the result from Bayesian circadian detection.
##' @title False Discovery Rate (FDR) Calculation for Bayesian Circadian Rhythmicity Detection.
##' @param r Gene indicator matrix. Each row is the result of 1 MCMC iteration. Each column is gene indicator for each gene, 1 = circadian gene, 0 = non-circadian gene.
##' @param G Number of total genes. 
##' @param G1 Number of true circadian genes. 
##' @param M number of MCMC iterations in your Bayesian circadian detection method. 

##' 
##' @return A list of Bayesian q-values and running time in seconds. 
##' \item{qvalues}{A vector of q-values for each gene.}
##' \item{alpha0}{A matrix contains all non-circadian genes alpha values for each iteration.}
##' \item{sigma20}{A matrix contains all non-circadian genes sigma square values for each iteration.}
##' 
##' @author Zhiguang Huo, Haocheng Ding, Lingsong Meng
##' @export
##' @examples
##' set.seed(32611)
bayesFDR <- function(r, G, G1, M){
  draw_r <- r
  
  p_L1 <- colSums(draw_r[(M/2+1):M,])/(M/2)
  p_g <- 1-p_L1
  
  BayesianFDR_qvalues <- numeric(0)
  for(i in 1:length(p_g)){
    BayesianFDR_qvalues[i] <- sum((p_g<=p_g[i])*p_g)/sum(p_g<=p_g[i])
  }
  
  # True positive
  TP_BayFDR <- sum(BayesianFDR_qvalues[1:G1]<0.05)
  # False positive
  FP_BayFDR <- sum(BayesianFDR_qvalues[(G1+1):G]<0.05)
  
  # Bayesian_FDR FDR
  Bayesian_FDR <- FP_BayFDR/(FP_BayFDR+TP_BayFDR)
  #########################
  return(list(FDR = Bayesian_FDR, TP = TP_BayFDR, FP = FP_BayFDR))
}