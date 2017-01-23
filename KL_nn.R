library(flexclust)

## one saples
KL.mv<- function(p.samp, q.dens, k = 1){
  d<- ncol(p.samp)
  n<- nrow(p.samp)
  
  eudist<- as.matrix(dist(p.samp, method = "euclidean"))
  diag(eudist)<- NA
  
  if(k == 1){
    r<- apply(eudist, 1, function(x) min(x, na.rm = T))
  }else{
    r<- apply(eudist, 1, function(x) sort(x)[k])
  }
  cont<- k/(n-1) * gamma(d/2 + 1) / pi^(d/2)
  
  mean(
    log(
      (cont/r^d)/q.dens(p.samp)
    )
  ) + digamma(k)
}

## two samples
KL.nn<- function(q.samp, pi.samp, k = 1){
  d.q<- ncol(q.samp)
  n.q<- nrow(q.samp)
  
  d.pi<- ncol(pi.samp)
  n.pi<- nrow(pi.samp)
  
  if(d.q != d.pi){
    stop("dimention error!")
  }else{
    d = d.q
  }
  
  ## density estimation of q
  eudist.q<- as.matrix(dist(q.samp, method = "euclidean"))
  diag(eudist.q)<- NA
  
  if(k == 1){
    r<- apply(eudist.q, 1, function(x) min(x, na.rm = T))
  }else{
    r<- apply(eudist.q, 1, function(x) sort(x)[k])
  }
  constant<- k/(n.q-1) * gamma(d/2 + 1) / pi^(d/2)
  dens.q<- constant/r^d
  
  ## density estimation of pi
  eudist.pi<- dist2(q.samp, pi.samp, method = "euclidean")
  
  if(k == 1){
    r<- apply(eudist.pi, 1, function(x) min(x, na.rm = T))
  }else{
    r<- apply(eudist.pi, 1, function(x) sort(x)[k])
  }
  constant<- k/(n.pi) * gamma(d/2 + 1) / pi^(d/2)
  dens.pi<- constant/r^d
    
  mean(
    log(
      dens.q/dens.pi
    )
  )
}