#-------Example 1----------
# library packages
rm(list=ls())
library(MASS)
library(splines)
library(alabama)

# input functions 
filepathfun = '~/Functions/'
source(paste(filepathfun,'three_step.R',sep=''))
source(paste(filepathfun,'myknots.R',sep=''))
source(paste(filepathfun,'sim1data.R',sep=''))


# sample size
n = 300
nsim = 500
q = 3
kseq = seq(2, 4, by=1)

# result output
mse_est_beta1f = rep(0, nsim)
mse_est_beta2f = rep(0, nsim)
mse_est_beta3f = rep(0, nsim)
mse_est_g1f = rep(0, nsim)
mse_est_g2f = rep(0, nsim)



for(ns in 1:nsim){
  Dat = simdata(n)
  Y=Dat$Y;X=Dat$X;Z=Dat$Z;Beta=Dat$Beta;U=Dat$U;g=Dat$g
  kk = myknots(kseq, q, y=Y, z=Z, x=X, u=U)
  L = Spest1(q, k=kk$optk, k1=kk$optk1, y=Y, z=Z, x=X, u=U)
  mse_est_beta1f[ns] = sqrt(mean((Beta[,1] -  L$Betanew[,1])^2))
  mse_est_beta2f[ns] = sqrt(mean((Beta[,2] -  L$Betanew[,2])^2))
  mse_est_beta3f[ns] = sqrt(mean((Beta[,3] -  L$Betanew[,3])^2))
  mse_est_g1f[ns] = sqrt(mean((g[,1] - L$gnew[,1])^2))
  mse_est_g2f[ns] = sqrt(mean((g[,2] - L$gnew[,2])^2))
  print(ns)
} 



mean(mse_est_beta1f);sd(mse_est_beta1f)
mean(mse_est_beta2f);sd(mse_est_beta2f)
mean(mse_est_beta3f);sd(mse_est_beta3f)
mean(mse_est_g1f);sd(mse_est_g1f)
mean(mse_est_g2f);sd(mse_est_g2f)




