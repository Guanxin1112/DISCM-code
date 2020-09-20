myknots<-function(kseq, q, y, z, x, u){
  #choose optimal knots for dynamic single index varying-coefficient model
  #kesq: potential knots
  dx = ncol(x); dz = ncol(z); n = nrow(z)
  bicseq = matrix(0, length(kseq), length(kseq))
  for(i in 1:length(kseq)){
    for(j in 1:length(kseq)){
      Ures = Spest1(q=q, k=kseq[i], k1= kseq[j], y, z, x, u)$res
      rss = sum(Ures^2)/n
      para = dx*(q+kseq[j]) + dz*(q+kseq[i])
      bicseq[i,j] = log(rss) + para * log(n)/n
    }
  }
  optind1 = which(bicseq == min(bicseq),arr.ind=T)
  optk = kseq[optind1[1]]; optk1 = kseq[optind1[2]]
  list(optk=optk,optk1=optk1,bicseq=bicseq)
}

