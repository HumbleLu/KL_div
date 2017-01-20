library(flexclust)

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