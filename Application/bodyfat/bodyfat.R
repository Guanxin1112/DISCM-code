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


bodyfatdata<-readr::read_csv("bodyfat.csv")
D1 = bodyfatdata$density
D2 = bodyfatdata$perfat 
D3 = bodyfatdata$weight
D10 = bodyfatdata$age

D4 = bodyfatdata$chest
D5 = bodyfatdata$abdomen
D6 = bodyfatdata$hip
D7 = bodyfatdata$thigh 
D8 = bodyfatdata$FOREARM 
D9 = bodyfatdata$WRIST 
D11 = bodyfatdata$neck


n = length(D1)
Brozek_perfat = rep(457,n)/D1 - rep(414.2,n)
ind = which(Brozek_perfat<2| Brozek_perfat>40 | D3 > 250)

Brozek_perfat = Brozek_perfat[-ind]
Y = log(Brozek_perfat) 
n = n - length(ind)


D4t = (D4[-ind] - mean(D4[-ind])) / sd(D4[-ind])
D5t = (D5[-ind] - mean(D5[-ind])) / sd(D5[-ind])
D6t = (D6[-ind] - mean(D6[-ind])) / sd(D6[-ind])
D7t = (D7[-ind] - mean(D7[-ind])) / sd(D7[-ind])

D9t = (D9[-ind] - mean(D9[-ind])) / sd(D9[-ind])
D11t = (D11[-ind] - mean(D11[-ind])) / sd(D11[-ind])


Z = cbind(D5t, D4t, D6t) # sort by correlation with y
X1 = rep(1, n)

X2 = (rep(1,n)-Brozek_perfat/100) * D3[-ind]
X21 = (X2 - mean(X2)) / sd(X2)



U = X21
#sort data
SU = sort(U);index = order(U)
SX = matrix(X1, ncol=1)
SZ = Z[index,]
SY = Y[index]




q = 3
kseq = seq(2,5,by=1)
kk = myknots(kseq, q, y=SY, z=SZ, x=SX, u=SU)
# kk$optk =2;kk$optk=2
Lbody = Spest1(q, k=kk$optk, k1=kk$optk1, y=SY, z=SZ, x=SX, u=SU)
shapiro.test(Lbody$res)
sd(Lbody$res)  # 0.2140551

qqnorm(Lbody$res,main='QQ plot of Sample Data versus Standard Normal',cex.main=1)
qqline(Lbody$res)

#---------------------------------------------------------------------
bk = betaknot(kseq, q, y=SY, z=SZ, x=SX)
Lcons = betaesti(q, k1=bk$optk1, y=SY, z=SZ, x=SX)


#--------------------model check
beta_fun = Lbody$Betanew
haterr = Lbody$res - mean(Lbody$res)
beta_cons = Lcons$betaold
beta_consmm = sapply(beta_cons, rep, n)
TS_sample =  sum(apply((beta_fun - beta_consmm)^2, 2 ,mean))

B = 500
t = Sys.time()
TS_boot = rep(0,B)
TS_boot1 = rep(0,B)
TS_boot2 = rep(0,B)
TS_boot3 = rep(0,B)
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

sum(TS_boot > TS_sample) / B



