betaknot<-function(kseq, q, y, z, x){
  #choose optimal knots for single index varying-coefficient model
  #kesq: potential knots
  dx = ncol(x);dz = ncol(z);n = nrow(z)
  bicseq = rep(0, length(kseq))
  for(i in 1:length(kseq)){
    Ures = betaesti(q, kseq[i], y, z, x)$res
    rss = mean(Ures^2)
    para = dx*(q+kseq[i])
    bicseq[i] = log(rss) + para * log(n)/n
  }
  optind1 = which.min(bicseq)
  optk1 = kseq[optind1[1]]
  list(optk1=optk1,bicseq=bicseq)
}
