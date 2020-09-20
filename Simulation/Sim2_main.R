#-------Example 2----------
# library packages
rm(list=ls())
library(MASS)
library(splines)
library(alabama)


# input functions 
filepathfun = '~/Functions/'
source(paste(filepathfun,'three_step.R',sep=''))
source(paste(filepathfun,'myknots.R',sep=''))
source(paste(filepathfun,'betaesti.R',sep=''))
source(paste(filepathfun,'betaknot.R',sep=''))
source(paste(filepathfun,'sim2data.R',sep=''))



# sample size
B = 500;
nsim = 300
n = 700
q = 3
kseq = seq(1, 4,by=1)

# result output
pcu = rep(0,nsim)


for(ns in 1:nsim){
  set.seed(1000*ns)
  Dat = data_test(n,rr=0)
  Y=Dat$Y;X=Dat$X;Z=Dat$Z;U=Dat$U
  kk = myknots(kseq, q, y=Y, z=Z, x=X, u=U)
  Lfun = Spest1(q, k=kk$optk, k1=kk$optk1, y=Y, z=Z, x=X, u=U)
  beta_fun = Lfun$Betanew
  haterr = Lfun$res - mean(Lfun$res)
  bk = betaknot(kseq, q, y=Y, z=Z, x=X)
  Lcons = betaesti(q, k1=bk$optk1, y=Y, z=Z, x=X)
  beta_cons = Lcons$betaold
  beta_consmm = sapply(beta_cons, rep, n)
  TS_sample =  sum(apply((beta_fun - beta_consmm)^2, 2 ,mean))
  TS_boot = rep(0,B)
  for(b in 1:B){  
    err = haterr * rnorm(n,0,1)
    booty =  Lcons$fit + err
    Lfunb = Spest1(q, k=kk$optk, k1=kk$optk1, y=booty, z=Z, x=X, u=U)
    beta_funb = Lfunb$Betanew
    Lconsb = betaesti(q, k1=bk$optk1, y=booty, z=Z, x=X)
    beta_consb = Lconsb$betaold
    beta_consmmb = sapply(beta_consb, rep, n)
    TS_boot[b] =  sum(apply((beta_funb - beta_consmmb)^2, 2 ,mean))
  }
  pcu[ns] = sum(TS_boot > TS_sample) / B
}

sum(pcu<=0.05) / nsim  
 
