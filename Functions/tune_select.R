SCAD<-function(u,tune){
  #the SCAD penalization function 
  ca1 = (0 <= u) * (u <= tune) * tune * u 
  ca2 = -(u > tune) * (u < 3.7*tune) * (u^2 - 2*3.7*tune*u + tune^2)/(2*2.7)
  ca3 = (u >= 3.7*tune) * 4.7 * tune^2 / 2
  ca1 + ca2 + ca3
}


devSCAD<-function(u,tune){
  #compute the first derivative of SCAD penalization function 
  tune*(u <= tune) + (u > tune)*max(3.7*tune-u, 0)/2.7
}



Penalest<-function(tune, deltaini, q, k, k1, y, z, x, u){
  #obtain the penalized estimation and identify constant coefficient functions
  #get the penalized function estimation and the norm of the first derivative function of (\beta_{k}) function
  dz = ncol(z); dx = ncol(x); n = nrow(z)
  tmp = rep(0, k)  # internal knots
  for(j in 1:k){ tmp[j] = quantile(sort(u), j/(k+1))}
  b = max(u) + 10^(-10)
  a = min(u) - 10^(-10)
  ku = c(rep(a,q),t(tmp),rep(b,q))
  bpu  = splineDesign(knots = ku, x=u, ord = q)
  Rmatrix = t(bpu) %*% bpu / n  
  
  #initial value while looping
  Pcoeff = deltaini
  coeff = Pcoeff
  diff = 1 
  zeroind  = 0; varyind = seq(1,dz,by=1)
  while(diff > 1e-2){
    oldcoeff = matrix(0, nrow = (q+k), ncol=dz)
    oldcoeff1 = matrix(coeff, nrow = (q+k), byrow= F)
    oldcoeff[,varyind] = oldcoeff1
    Betaold = matrix(0, n, dz)
    for(l in varyind){Betaold[,l] =  as.vector(bpu %*% oldcoeff[,l])}
    oldnorm = sqrt(sum(colMeans(Betaold^2)))
    Betaold = Betaold / oldnorm
    
    Pfn<-function(ddd){
      #standardized estimation of (\beta_{k}) function
      ddd1 = matrix(ddd, nrow = (q+k), byrow= F)
      Betamm = matrix(0, n, dim(ddd1)[2])
      for(l in 1:dim(ddd1)[2]){Betamm[,l] =  as.vector(bpu %*% ddd1[,l])}
      norm = sqrt(sum(colMeans(Betamm^2)))
      Betamm = Betamm / norm
      betanorm = NULL
      for(l in 2:dim(ddd1)[2]){betanorm = c(betanorm, sqrt(mean(Betamm[,l]^2)))}
      
      if(0 %in% zeroind) {newz = z
      }else if(length(zeroind)==0){newz = z   #when zeroind = NULL
      }else {newz = z[,-zeroind]}
      
      wblock = rowSums(Betamm * newz)
      ww = rep(0, k1)  #xx1 : internal knots
      for(j in 1:k1){ww[j] = quantile(sort(wblock),j/(k1+1)) }
      b = max(wblock) + 10^(-10)
      a = min(wblock) - 10^(-10)
      kw = c(rep(a,q),t(ww),rep(b,q))
      bpwf  = splineDesign(knots = kw, x=wblock, ord = q)
      D = NULL
      for(l in 1:dx){D = cbind(D, bpwf * x[,l])}
      lambdaold = ginv(t(D) %*% D)%*% t(D) %*% y 
      res = y - D %*% lambdaold
      
      # minimize objective function
      PL1 = as.numeric(0.5 * t(res) %*% res)  + n * sum(SCAD(betanorm,tune))
      return(PL1)
    }
    
    
    Pgn<-function(ddd){
      ddd1 = matrix(ddd, nrow = (q+k), byrow= F)
      Betamm = matrix(0, n, dim(ddd1)[2])
      for(l in 1:dim(ddd1)[2]){Betamm[,l] =  as.vector(bpu %*% ddd1[,l])}
      norm = sqrt(sum(colMeans(Betamm^2)))
      Betamm = Betamm / norm
      betanorm = NULL
      for(l in 1:dim(ddd1)[2]){betanorm = c(betanorm, sqrt(mean(Betamm[,l]^2)))}
      
      if(0 %in% zeroind) {newz = z
      }else if(length(zeroind)==0){newz = z   #when zeroind = NULL
      }else {newz = z[,-zeroind]}
      
      wblock = rowSums(Betamm * newz)
      ww = rep(0, k1)  #xx1 : internal knots
      for(j in 1:k1){ww[j] = quantile(sort(wblock),j/(k1+1)) }
      b = max(wblock) + 10^(-10)
      a = min(wblock) - 10^(-10)
      kw = c(rep(a,q),t(ww),rep(b,q))
      bpwf  = splineDesign(knots = kw, x=wblock, ord = q)
      D = NULL
      for(l in 1:dx){D = cbind(D, bpwf * x[,l])}
      lambdaold = ginv(t(D) %*% D)%*% t(D) %*% y 
      res = y - D %*% lambdaold
      
      Phi = NULL
      for(l in 1:dim(newz)[2]){Phi = cbind(Phi, bpu * newz[, l])}
      
      eta = ginv(t(D) %*% D)%*% t(D) %*% Phi
      tildePhi =  Phi- D %*% eta
      lambdaold1 = matrix(lambdaold, ncol=dx, byrow = F)
      dev = matrix(1, nrow=n, ncol=1)
      bpwdev  = splineDesign(knots = kw, x=wblock, ord = q, derivs=dev)
      ghatdev = rep(0,n)
      for(l in 1:dx){
        ghatdev0 =  bpwdev %*% lambdaold1[,l]
        ghatdev = ghatdev + ghatdev0 * x[,l]
      }
      dLn1 = 0
      for(i in 1:n){dLn1 = dLn1 + res[i] * ghatdev[i] * tildePhi[i,]}
      dLn1 = - dLn1 
      dvec = rep(0, (q+k))
      for(l in 2:dim(ddd1)[2]){
        tp = n * (devSCAD(betanorm[l], tune) / betanorm[l]) * Rmatrix %*%  ddd1[,l]
        dvec = c(dvec, tp)
      }
      dLn1 + dvec
    }
    
    
    Amat = matrix(0, nrow = length(coeff), ncol = length(coeff))
    for(i in 2:(q+k)){Amat[i,(i-1)] = -1}
    for(i in 2:(q+k)){Amat[i,i] = 1}
    bvec = rep(0, length(coeff)) - 1e-10    #ui %*% theta - ci >=0 
    R = constrOptim(theta=coeff, f=Pfn, grad = Pgn, ui = Amat, ci = bvec) 
    temp = R$par    
     
    #the construction of temp is used to make sure that zeroind is nondecreasing
    tempm = matrix(0, nrow = (q+k), ncol=dz)
    tempm1 = matrix(temp, nrow = (q+k), byrow= F)
    tempm[,varyind] = tempm1
    Betamm1 = matrix(0, n, dz)
    for(l in varyind){Betamm1[,l] =  as.vector(bpu %*% tempm[,l])}
    norm1 = sqrt(sum(colMeans(Betamm1^2)))
    Betamm1 = Betamm1 / norm1
    
    betanorm1 = NULL
    for(l in 1:dz){betanorm1 = c(betanorm1, sqrt(mean(Betamm1[,l]^2)))}
    #new zeroind and varyind by testing if some beta is zero 
    zeroind1 = which(betanorm1 < 0.01)
    if(0 %in% zeroind1){varyind1 = seq(1,dz)
    }else if(length(zeroind1) == 0){varyind1 = seq(1,dz)
    }else{varyind1 = seq(1,dz)[-zeroind1]}
    
    #update zeroind,varyind and coeff
    zeroind = zeroind1
    varyind = varyind1
    if(0 %in% zeroind){coeff = as.vector(tempm)
    }else if(length(zeroind) == 0) {coeff = as.vector(tempm)
    }else{coeff = as.vector(tempm[,-zeroind])}
    
    ##test the convergence condition,i.e update the diffremce norm
    tempm2 = matrix(0, nrow = (q+k), ncol=dz)
    tempm21 = matrix(coeff, nrow = (q+k), byrow= F)
    tempm2[,varyind] = tempm21
    Betanew = matrix(0, n, dz)
    for(l in varyind){Betanew[,l] = as.vector(bpu %*% tempm2[,l])}
    newnorm = sqrt(sum(colMeans(Betanew^2)))
    Betanew = Betanew / newnorm
    
    diffBeta = Betanew - Betaold
    diff = sqrt(sum(colMeans(diffBeta^2)))
    if(length(varyind)==1) break
  }
  
  varyterm = length(varyind)
  wblock = rowSums(Betanew * z)
  ww = rep(0, k1)  #xx1 : internal knots
  for(j in 1:k1){ww[j] = quantile(sort(wblock),j/(k1+1)) }
  b = max(wblock) + 10^(-10)
  a = min(wblock) - 10^(-10)
  kw = c(rep(a,q),t(ww),rep(b,q))
  nbpw  = splineDesign(knots = kw, x=wblock, ord = q)
  D = NULL
  for(l in 1:dx){D = cbind(D, nbpw * x[,l])}
  lambdanew = ginv(t(D) %*% D)%*% t(D) %*% y 
  lambdanew1 = matrix(lambdanew, ncol=dx, byrow = F)
  gnew = matrix(0, n, dx)
  for(l in 1:dx){gnew[,l] = nbpw %*% lambdanew1[,l]}
  
  fit = rowSums(gnew * x)
  rss = sum((y - fit)^2) / n
  mybic = log(rss) + varyterm * (q+k) * log(n)/n
  
  list(Betanew = Betanew, gnew=gnew, varyind = varyind, mybic = mybic)
}






OptimPenalest<-function(tuneseq,q, k, k1, y, z, x, u){
  #choose optimal tuning parameter
  mybic = rep(0, length(tuneseq))
  deltaold = bfun_ini(q, k, y, z, x, u)
  for(i in 1:length(tuneseq)){
    mybic[i] = Penalest(tuneseq[i],deltaold,q, k, k1, y, z, x, u)$mybic
  }
  sind = order(mybic)
  optune = tuneseq[sind[1]]
  PP = Penalest(tune=optune,deltaold,q, k, k1, y, z, x, u)
  Betanew = PP$Betanew
  varyind = PP$varyind
  gnew = PP$gnew
  list(optune=optune,varyind=varyind,Betanew=Betanew,gnew=gnew)
}



