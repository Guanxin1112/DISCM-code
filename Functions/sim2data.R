#index function
g1f<-function(u) 10*exp(5*u)/(1+exp(5*u)) 
g2f<-function(u) 5*sin(pi*u)


#varying coefficient function
beta1f<-function(u)  1 + u^2
beta2f<-function(u)  (1-u)^2 
beta3f<-function(u)  -1+4*(u-0.5)^2


data_test<-function(n, rr){
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
  
  U = seq_along(1:n)/n
  
  c1 = integrate(f = beta1f, 0, 1)$value
  c2 = integrate(f = beta2f, 0, 1)$value
  c3 = integrate(f = beta3f, 0, 1)$value
  
  beta1ft<-function(u)  rr * beta1f(u) + c1
  beta2ft<-function(u)  rr * beta2f(u) + c2
  beta3ft<-function(u)  rr * beta3f(u) + c3
  
  #standardized---identifiability condition
  temp1 = beta1ft(U)
  temp2 = beta2ft(U)
  temp3 = beta3ft(U)
  
  temp11 = cbind(temp1, temp2, temp3)
  ident = sqrt(sum(colMeans(temp11^2)))
  Beta = temp11 / ident
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
  
  
  list(Y=Y,X=X,Z=Z,U=U,Beta=Beta,g=g,eps=eps)
}

