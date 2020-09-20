#-------Example 3----------
# library packages
rm(list = ls())
library(MASS)
library(splines)
library(alabama)


# input functions 
filepathfun = '~/Functions/'
source(paste(filepathfun,'three_step.R',sep=''))
source(paste(filepathfun,'myknots.R',sep=''))
source(paste(filepathfun,'sim2data.R',sep=''))
source(paste(filepathfun,'tune_select.R',sep=''))


# sample size
n = 500
nsim = 500
q = 3
k = floor(n^(1/(2*q+1)))
k1 =  floor(n^(1/(2*q+1)))


# result output
mse_est_beta1f = rep(0, nsim)
mse_est_beta2f = rep(0, nsim)
mse_est_beta3f = rep(0, nsim)
mse_est_beta4f = rep(0, nsim)
mse_est_beta5f = rep(0, nsim)
mse_est_beta6f = rep(0, nsim)
mse_est_g1f = rep(0, nsim)
mse_est_g2f = rep(0, nsim)

mse_est_Pbeta1f = rep(0, nsim)
mse_est_Pbeta2f = rep(0, nsim)
mse_est_Pbeta3f = rep(0, nsim)
mse_est_Pbeta4f = rep(0, nsim)
mse_est_Pbeta5f = rep(0, nsim)
mse_est_Pbeta6f = rep(0, nsim)
mse_est_Pg1f = rep(0, nsim)
mse_est_Pg2f = rep(0, nsim)
est_varyind = matrix(0,6,nsim)
optune =  rep(0, nsim)
for(ns in 1:nsim){
  set.seed(1000*ns)
  Dat = simpdata(n)
  Y=Dat$Y;X=Dat$X;Z=Dat$Z;Beta=Dat$Beta;U=Dat$U;g=Dat$g;dz=ncol(Z)
  #three-step spline estimation of varying-coefficient additive model
  L = Spest1(q, k, k1, y=Y, z=Z, x=X, u=U)

  mse_est_beta1f[ns] = sqrt(mean((Beta[,1] - L$Betanew[,1])^2))
  mse_est_beta2f[ns] = sqrt(mean((Beta[,2] - L$Betanew[,2])^2))
  mse_est_beta3f[ns] = sqrt(mean((Beta[,3] - L$Betanew[,3])^2))
  mse_est_beta4f[ns] = sqrt(mean((Beta[,4] - L$Betanew[,4])^2))
  mse_est_beta5f[ns] = sqrt(mean((Beta[,5] - L$Betanew[,5])^2))
  mse_est_beta6f[ns] = sqrt(mean((Beta[,6] - L$Betanew[,6])^2))
  mse_est_g1f[ns] = sqrt(mean((g[,1] - L$gnew[,1])^2))
  mse_est_g2f[ns] = sqrt(mean((g[,2] - L$gnew[,2])^2))
  
  #model identification
  tuneseq = seq(0.05, 1, by=0.01)
  PL = OptimPenalest(tuneseq,  q, k, k1, y=Y, z=Z, x=X, u=U)
  optune[ns] = PL$optune
  est_varyind[,ns] = c(PL$varyind,rep(0, ncol(Z)-length(PL$varyind)))
  mse_est_Pbeta1f[ns] = sqrt(mean((Beta[,1] - PL$Betanew[,1])^2))
  mse_est_Pbeta2f[ns] = sqrt(mean((Beta[,2] - PL$Betanew[,2])^2))
  mse_est_Pbeta3f[ns] = sqrt(mean((Beta[,3] - PL$Betanew[,3])^2))
  mse_est_Pbeta4f[ns] = sqrt(mean((Beta[,4] - PL$Betanew[,4])^2))
  mse_est_Pbeta5f[ns] = sqrt(mean((Beta[,5] - PL$Betanew[,5])^2))
  mse_est_Pbeta6f[ns] = sqrt(mean((Beta[,6] - PL$Betanew[,6])^2))
  mse_est_Pg1f[ns] = sqrt(mean((g[,1] - PL$gnew[,1])^2))
  mse_est_Pg2f[ns] = sqrt(mean((g[,2] - PL$gnew[,2])^2))
  print(ns)
}





true_index = c(1,3,5);wrong_index = c(2,4,6)
low_est = 0
true_est = 0
high_est = 0
ar = rep(0, nsim)
ai = rep(0, nsim)
for(ns in 1:nsim){
  choose = est_varyind[,ns][est_varyind[,ns]!=0]
  if(length(choose) < length(true_index)){any(choose %in% true_index);
    low_est = low_est + 1
  }else if(length(choose) == length(true_index)){all(choose %in% true_index);
    true_est = true_est + 1
  }else if(length(choose) > length(true_index)){all(choose %in% true_index);
    high_est = high_est + 1}
  ar[ns] = length(intersect(choose, true_index))
  ai[ns] = length(intersect(choose, wrong_index))
}


low_est/nsim
true_est/nsim
high_est/nsim
mean(ar)
mean(ai)

