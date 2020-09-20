#-------Body fat data----------
# library packages
rm(list = ls())
library(MASS)
library(splines)
library(alabama)

# input functions 
filepathfun = '~/Functions/'
source(paste(filepathfun,'three_step.R',sep=''))
source(paste(filepathfun,'myknots.R',sep=''))
source(paste(filepathfun,'betaesti.R',sep=''))
source(paste(filepathfun,'betaknot.R',sep=''))



bostondata<-readr::read_csv("boston.csv")
D1 = bostondata$MEDV
D2 = bostondata$LSTAT
D3 = bostondata$PTRATIO
D4 = bostondata$NOX
D5 = bostondata$TAX
D6 = bostondata$AGE
D7 = bostondata$CRIM
D8 = bostondata$RM
D9 = bostondata$DIS 


n = length(D1)
Y = D1 - mean(D1)
U = D9
Z = cbind((D8-mean(D8))/sd(D8), (D2-mean(D2))/sd(D2), (D5-mean(D5))/sd(D5))
X1 = log(D7)

index = order(U)
SU = U[index]
SY = Y[index]
SZ = Z[index, ]

SX = cbind(rep(1, n), X1[index])

q = 3
kseq = seq(2,5,by=1)
kk = myknots(kseq, q, y=SY, z=SZ, x=SX, u=SU)
Lboston= Spest1(q, k=kk$optk, k1=kk$optk1, y=SY, z=SZ, x=SX, u=SU)
sd(Lboston$res) 


qqnorm(Lboston$res,main='QQ plot of Sample Data versus Standard Normal',cex.main=1)
qqline(Lboston$res)



#--------------------model check
bk = betaknot(kseq, q, y=SY, z=SZ, x=SX)
Lcons = betaesti(q, k1=bk$optk1, y=SY, z=SZ, x=SX)

beta_fun = Lboston$Betanew
haterr = Lboston$res - mean(Lboston$res)
beta_cons = Lcons$betaold
beta_consmm = sapply(beta_cons, rep, n)
TS_sample =  sum(apply((beta_fun - beta_consmm)^2, 2 ,mean))

B = 500
t = Sys.time()
TS_boot = rep(0,B)

for(b in 1:B){  
  set.seed(100*b)
  err = haterr * rnorm(n,0,1)
  booty =  Lcons$fit + err
  Lfunb = Spest1(q, k=kk$optk, k1=kk$optk1, y=booty, z=SZ, x=SX, u=SU)
  beta_funb = Lfunb$Betanew
  Lconsb = betaesti(q, k1=bk$optk1, y=booty, z=SZ, x=SX)
  beta_consb = Lconsb$betaold
  beta_consmmb = sapply(beta_consb, rep, n)
  TS_boot[b] =  sum(apply((beta_funb - beta_consmmb)^2, 2 ,mean))
}
Sys.time() -t 

sum(TS_boot > TS_sample) / B  # 0.058




