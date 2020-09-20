bfun_ini<-function(q, k, y, z, x, u){
  #obtain the initial estimate of varying-coefficient function by assuming that g() are linear functions
  #q: order of B-spline for both varying-coefficient function and index function
  #k: knots for varying-coefficient function
  dz = ncol(z);dx = ncol(x);n = nrow(z)
  tmp = rep(0, k)  # internal knots
  for(j in 1:k){ tmp[j] = quantile(sort(u), j/(k+1))}
  b = max(u) + 10^(-10)
  a = min(u) - 10^(-10)
  ku = c(rep(a,q),t(tmp),rep(b,q))
  bpu  = splineDesign(knots = ku, x=u, ord = q)
  ##Initialization step  estimation
  Phi = NULL
  rex = NULL
  for(l in 1:dz){ Phi = cbind(Phi, bpu * z[, l]) }
  for(l in 1:dx){ rex = cbind(rex, Phi * x[, l]) }
  rex = cbind(x, rex)
  atheta = as.vector(ginv(t(rex) %*% rex)%*% t(rex) %*% y) #ridge regression
  
  part_atheta = atheta[-(1:dx)]
  matrix_atheta = matrix(part_atheta, ncol=dx, byrow= F)
  delta0 = rowMeans(matrix_atheta)
  deltaold = c(sort(delta0[1:(q+k)]), delta0[-(1:(q+k))])
 
  return(deltaold)
} 


Est_nlfun1<-function(deltaold, q, k, k1, y, z, x, u){
  #It contains the iterative estimation process and get the final estimation of functions
  #deltaold: the initial estimate of the varying-coefficient function
  #q: order of B-spline for both varying-coefficient function and index function
  #k: knots for varying-coefficient function
  #k1: knots for index function
  dz = ncol(z); dx = ncol(x); n = nrow(z)
  tmp = rep(0, k)  # internal knots
  for(j in 1:k){ tmp[j] = quantile(sort(u), j/(k+1))}
  b = max(u) + 10^(-10)
  a = min(u) - 10^(-10)
  ku = c(rep(a,q),t(tmp),rep(b,q))
  bpu  = splineDesign(knots = ku, x=u, ord = q)
  Phi = NULL
  for(l in 1:dz){Phi = cbind(Phi, bpu * z[, l])}
  
  fn<-function(ddd){
    ddd1 = matrix(ddd, ncol = dz, byrow= F)
    #standardized
    Betaold = matrix(0, n, dz)
    for(l in 1:dz){Betaold[,l] =  as.vector(bpu %*% ddd1[,l])}
    norm = sqrt(sum(colMeans(Betaold^2)))
    Betaold = Betaold / norm
    #
    wblock = rowSums(Betaold * z)
    ww = rep(0, k1)  #xx1 : internal knots
    for(j in 1:k1){ww[j] = quantile(sort(wblock),j/(k1+1)) }
    b = max(wblock) + 10^(-10)
    a = min(wblock) - 10^(-10)
    kw = c(rep(a,q),t(ww),rep(b,q))
    bpwf  = splineDesign(knots = kw, x=wblock, ord = q)
    D = NULL
    for(l in 1:dx){D = cbind(D, bpwf * x[,l])}
    lambdaold = as.vector(ginv(t(D) %*% D)%*% t(D) %*% y) 
    res = as.vector(y - D %*% lambdaold)
    # minimize objective function
    L1 = as.numeric(0.5*t(res) %*% (res))
    return(L1)
  }
  
  gn<-function(ddd){
    ddd1 = matrix(ddd, ncol = dz, byrow= F)
    #standardized
    Betaold = matrix(0, n, dz)
    for(l in 1:dz){Betaold[,l] =  as.vector(bpu %*% ddd1[,l])}
    norm = sqrt(sum(colMeans(Betaold^2)))
    Betaold = Betaold / norm
    #
    wblock = rowSums(Betaold * z)
    ww = rep(0, k1)  #xx1 : internal knots
    for(j in 1:k1){ww[j] = quantile(sort(wblock),j/(k1+1)) }
    b = max(wblock) + 10^(-10)
    a = min(wblock) - 10^(-10)
    kw = c(rep(a,q),t(ww),rep(b,q))
    bpwf  = splineDesign(knots = kw, x=wblock, ord = q)
    D = NULL
    for(l in 1:dx){D = cbind(D, bpwf * x[,l])}
    lambdaold = as.vector(ginv(t(D) %*% D)%*% t(D) %*% y) 
    res = as.vector(y - D %*% lambdaold)
    
    
    eta = ginv(t(D) %*% D)%*% t(D) %*% Phi
    tildePhi =  Phi- D %*% eta
    lambdaold1 = matrix(lambdaold, ncol=dx, byrow = F)
    dev = matrix(1, nrow=n, ncol=1)
    bpwdev  = splineDesign(knots = kw, x=wblock, ord = q, derivs=dev)
    ghatdev = rep(0,n)
    for(l in 1:dx){
      ghatdev0 =  as.vector(bpwdev %*% lambdaold1[,l])
      ghatdev = ghatdev + ghatdev0 * x[,l]
    }
    dLn1 = 0
    for(i in 1:n){dLn1 = dLn1 + res[i] * ghatdev[i] * tildePhi[i,]}
    dLn1 = - dLn1
    dLn1
  }
  
  #beta1(u) is restricted to be nondecreasing
  Amat = matrix(0, nrow = dz*(q+k), ncol = dz*(q+k))
  for(i in 2:(q+k)){Amat[i,(i-1)] = -1}
  for(i in 2:(q+k)){Amat[i,i] = 1}
  bvec = rep(0, dz*(q+k))  - 1e-10
  #the constrained optimization problem
  R = constrOptim(theta=deltaold, f=fn, grad = gn, ui = Amat, ci = bvec, method = "BFGS")
  conver = R$convergence
  deltanew = R$par
  
  ###standardized Beta
  deltanew1 = matrix(deltanew, ncol = dz, byrow= F)
  Betanew = matrix(0, n, dz)
  for(l in 1:dz){Betanew[,l] =  as.vector(bpu %*% deltanew1[,l])}
  normnew = sqrt(sum(colMeans(Betanew^2)))
  Betanew = Betanew /normnew
  
  wblocknew = rowSums(Betanew * z)
  nww = rep(0, k1)  #xx1 : internal knots
  for(j in 1:k1){nww[j] = quantile(sort(wblocknew),j/(k1+1))}
  b = max(wblocknew) + 10^(-10)
  a = min(wblocknew) - 10^(-10)
  kw = c(rep(a,q),t(nww),rep(b,q))
  nbpw  = splineDesign(knots = kw, x=wblocknew, ord = q)
  nD = NULL
  for(l in 1:dx){nD = cbind(nD, nbpw * x[,l])}
  lambdanew = as.vector(ginv(t(nD) %*% nD)%*% t(nD) %*% y)
  lambdanew1 = matrix(lambdanew, ncol=dx, byrow = F)
  gnew = matrix(0, n, dx)
  for(l in 1:dx){gnew[,l] = as.vector(nbpw %*% lambdanew1[,l])}
  
  list(Betanew = Betanew, gnew=gnew, conver = conver,deltanew = deltanew,lambdanew=lambdanew)
}


Spest1<-function(q, k, k1, y, z, x, u){
  #get the final estimate of varying-coefficient function and index function
  #step 0
  deltaini = bfun_ini(q, k, y, z, x, u) 
  #step1-3
  Res_step3 = Est_nlfun1(deltaini, q, k, k1, y, z, x, u)
  conver = Res_step3$conver
  Betanew = Res_step3$Betanew 
  gnew = Res_step3$gnew
  deltanew = Res_step3$deltanew
  lambdanew = Res_step3$lambdanew
  fit = rowSums(gnew * x)
  res = y - fit
  
  list(Betanew = Betanew, gnew=gnew, conver = conver, res=res, 
       fit=fit,deltanew = deltanew,lambdanew=lambdanew)
}


