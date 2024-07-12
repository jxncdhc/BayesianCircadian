# BayesianCircadian
Omics data circadian analysis with Bayesian approaches.


## Install This Package from github
* In R console

```{R}
library(devtools)
install_github("jxncdhc/BayesianCircadian") 
```

## Citation

## Full tutorial

## Short tutorial for Bayesian circadian detection

```{R}
set.seed(32611)
G=1000 
G1=100 
G0=G-G1
A1=3 
Mesor1=3
sigma21=1 
A0=0 
Mesor0=3 
sigma20=1
n=24
tt <- runif(n,0,24)
# generate circadian gene samples #
phi1 <- runif(G1,0,24)
w <- 2*pi/24
EE1 <- A1*cos(w*phi1)
FF1 <- A1*sin(w*phi1) 
yy1 <- matrix(0,G1,n)
for(i in 1:G1){
  yy1[i,] <- rnorm(n, EE1[i]*cos(w*tt) + FF1[i]*sin(w*tt) + Mesor1, sqrt(sigma21))
}  

# generate non-circadian gene samples #
phi0 <- runif(G0,0,24)
w <- 2*pi/24
EE0 <- A0*cos(w*phi0)
FF0 <- A0*sin(w*phi0) 
yy0 <- matrix(0,G0,n)
for(i in 1:G0){
  yy0[i,] <- rnorm(n, EE0[i]*cos(w*tt) + FF0[i]*sin(w*tt) + Mesor0, sqrt(sigma20))
}

yy <- rbind(yy1,yy0)

bayes_rhythmicity(yy,tt,M=1000)
```

## Short tutorial for Bayesian circadian detection with prior information
```{R}
set.seed(32611)
G=1000 
G1=100 
G0=G-G1
A1=3 
Mesor1=3
sigma21=1 
A0=0 
Mesor0=3 
sigma20=1
n=24
tt <- runif(n,0,24)
# generate circadian gene samples #
phi1 <- runif(G1,0,24)
w <- 2*pi/24
EE1 <- A1*cos(w*phi1)
FF1 <- A1*sin(w*phi1) 
yy1 <- matrix(0,G1,n)
for(i in 1:G1){
  yy1[i,] <- rnorm(n, EE1[i]*cos(w*tt) + FF1[i]*sin(w*tt) + Mesor1, sqrt(sigma21))
}  
# generate non-circadian gene samples #
phi0 <- runif(G0,0,24)
w <- 2*pi/24
EE0 <- A0*cos(w*phi0)
FF0 <- A0*sin(w*phi0) 
yy0 <- matrix(0,G0,n)
for(i in 1:G0){
  yy0[i,] <- rnorm(n, EE0[i]*cos(w*tt) + FF0[i]*sin(w*tt) + Mesor0, sqrt(sigma20))
}
yy <- rbind(yy1,yy0)
# Given prior information
prior_info <- c(rep(1,G1), rep(0,G0))
bayes_rhythmicity_prior(yy,tt,M=1000, prior_gene = prior_info)
```
