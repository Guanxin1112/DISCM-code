#data generation for example 1

#index function
g1f<-function(u) 10*exp(5*u)/(1+exp(5*u))
g2f<-function(u) 5*sin(pi*u)

#varying-coefficient function
beta1f<-function(u)  1+u^2
beta2f<-function(u)  (1-u)^2
beta3f<-function(u)  -1+4*(u-0.5)^2


simdata<-function(n){
  xsigma = matrix(c(1,0.3,0.3,1),2,2)
  X = mvrnorm(n, mu=rep(0,2), Sigma=xsigma)
  dx = 2
  
  zsigma = matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),3,3)
  sZ = mvrnorm(n, mu=rep(0,3), Sigma=zsigma)
  Z1 = pnorm(sZ[,1])-0.5
  Z2 = pnorm(sZ[,2])-0.5
  Z3 = pnorm(sZ[,3])-0.5
  Z = cbind(Z1,Z2,Z3)
  dz = 3
  
  U = sort(runif(n,min=0,max=1))
  
  #standardized---identifiability condition
  temp1 = beta1f(U)
  temp2 = beta2f(U)
  temp3 = beta3f(U)
  
  temp11 = cbind(temp1, temp2, temp3)
  norm = sqrt(sum(colMeans(temp11^2)))
  Beta = temp11 / norm
  
  V = rowSums(Z * Beta)
  #centralization
  g1 = g1f(V)-mean(g1f(V))
  g2 = g2f(V)-mean(g2f(V))
  g = cbind(g1, g2)
  my = rowSums(g * X)
  
  ##noise
  esigma = sqrt(0.2*var(my))
  eps = rnorm(n, 0, esigma)
  Y = my + eps
  
  list(Y=Y,X=X,Z=Z,U=U,Beta=Beta,g=g,my=my,eps=eps)
}

