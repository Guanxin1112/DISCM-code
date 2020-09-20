betaesti<-function(q, k1, y, z, x){
  #parametric estimation and nonparametric function estimation of single index varying-coefficient model
  #q: order of B-spline for single index varying model
  #k1: knots for index function
  dz = ncol(z); dx = ncol(x); n = nrow(z)
  zx = x
  for(l in 1:dx){
    zx1 = z * x[,l]
    zx = cbind(zx, zx1)
  }
  betaini1 = as.vector(solve(t(zx) %*% zx)%*% t(zx) %*% y)
  betaini2 = matrix(betaini1[-(1:dx)], ncol = dx, byrow= F)
  betaini = rowMeans(betaini2)
  betaini = betaini / sqrt(sum(betaini^2))
  betaini = betaini * sign(betaini[1])
  betaini1 = betaini[-1]
  
  fn<-function(rr){
    ttt = sqrt(abs(1-sum(rr^2)))
    mbeta = c(ttt, rr)
    wblock = as.vector(z %*% mbeta)
    ww = rep(0, k1)  
    for(j in 1:k1){ ww[j] = quantile(sort(wblock), j/(k1+1))}
    b = max(wblock) + 10^(-10)
    a = min(wblock) - 10^(-10)
    kw = c(rep(a,q),t(ww),rep(b,q))
    Bq  = splineDesign(knots = kw, x=wblock, ord = q)
    C = NULL
    for(l in 1:dx){ C = cbind(C, Bq * x[,l]) }
    lambdaold = as.vector(ginv(t(C) %*% C) %*% t(C) %*% y)
    res = as.vector(y - C %*% lambdaold)
    Ln =  as.numeric(0.5*t(res) %*% (res))
    Ln
  }
  
  gr<-function(rr){
    ttt = sqrt(abs(1-sum(rr^2)))
    mbeta = c(ttt, rr)
    wblock = as.vector(z %*% mbeta)
    ww = rep(0, k1)  
    for(j in 1:k1){ ww[j] = quantile(sort(wblock), j/(k1+1))}
    b = max(wblock) + 10^(-10)
    a = min(wblock) - 10^(-10)
    kw = c(rep(a,q),t(ww),rep(b,q))
    Bq  = splineDesign(knots = kw, x=wblock, ord = q)
    C = NULL
    for(l in 1:dx){ C = cbind(C, Bq * x[,l]) }
    lambdaold = as.vector(ginv(t(C) %*% C) %*% t(C) %*% y)
    res = as.vector(y - C %*% lambdaold)
    
    eta = ginv(t(C) %*% C) %*% t(C) %*% z
    ztilde = z - C %*% eta
    J11 = -t(rr) / ttt
    Jmatrix = rbind(J11, diag(length(J11)))
    
    znew = ztilde %*% Jmatrix
    lambdaold1 = matrix(lambdaold, ncol=dx, byrow = F)
    dev = matrix(1, nrow=n, ncol=1)
    bpwdev  = splineDesign(knots = kw, x=wblock, ord = q, derivs=dev)
    ghatdev = rep(0,n)
    for(l in 1:dx){
      ghatdev0 =  as.vector(bpwdev %*% lambdaold1[,l])
      ghatdev = ghatdev + ghatdev0 * x[,l]
    }
    dLn1 = 0
    for(i in 1:n){dLn1 = dLn1 + res[i] * ghatdev[i] * znew[i,] }
    dLn1 = - dLn1
    dLn1
  }
  
  #constraint matrix
  hin<-function(rr){
    ss = sum(rr^2)
    1-ss
  }
  
  L = constrOptim.nl(par=betaini1, fn=fn, gr=gr, hin=hin) 
  conver = L$convergence
  betaold0 = L$par
  betaold = c(sqrt(abs(1-sum(betaold0^2))),betaold0)
  
  wblock = as.vector(z %*% betaold)
  ww = rep(0, k1)  
  for(j in 1:k1){ ww[j] = quantile(sort(wblock), j/(k1+1))}
  b = max(wblock) + 10^(-10)
  a = min(wblock) - 10^(-10)
  kw = c(rep(a,q),t(ww),rep(b,q))
  nbpw  = splineDesign(knots = kw, x=wblock, ord = q)
  D = NULL
  for(l in 1:dx){D = cbind(D, nbpw * x[,l])}
  lambdanew = as.vector(ginv(t(D) %*% D)%*% t(D) %*% y)
  lambdanew1 = matrix(lambdanew, ncol=dx, byrow = F)
  
  gnew = matrix(0, n, dx)
  for(l in 1:dx){gnew[,l] = as.vector(nbpw %*% lambdanew1[,l])}
  fit = rowSums(gnew * x)
  res = y - fit
  
  list(betaold = betaold, spg = gnew, conver=conver, fit = fit, res=res)
  
}
