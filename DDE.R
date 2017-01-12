library(compiler)
library(MASS)
library(flexclust)

## G matrix
G_mtx<- cmpfun(
  function(X, sigma){
    d<- ncol(X)
    dist_mtx<- as.matrix(dist(X))
    (pi * sigma^2)^(d/2) * exp(dist_mtx^2/(-4*sigma^2))
  }
)

## psi: Guassian model
psi_func<- function(X, sigma){
  cmpfun(
    function(x){
      if(class(x) == "matrix"){
        exp((dist2(x, X)^2)/(-2 * sigma^2))
      }else{
        exp(apply(X, 1, function(y) sum((x-y)^2))/(-2 * sigma^2))
      }
    }
  )
}

## phi_j: here only implement the trace of Hessian matrix
phi_func<- function(X, j, sigma){
  ## jth element in the diagonal of Hessian matrix
  
  cmpfun(
    function(x){
      if(class(x) == "matrix"){
        psi_func(X, sigma)(x) * (dist2(x[, j, drop = F], X[, j, drop = F])^2 - sigma^2)/(sigma^4)
      }else{
        psi_func(X, sigma)(x) * ((x[j] - X[,j])^2 - sigma^2)/(sigma^4)
      }
    }
  )
}

# test
# X<- matrix(rnorm(100), 5, 20)
# G_mtx(X)
# xx<- matrix(rnorm(3 * 3), 3, 3)
# yy<- matrix(rnorm(5 * 3), 5, 3)
# g_func(yy, 1, 1, 1)(xx)
# psi_func(yy, 1)(xx)
# psi2_func(yy, 1, 1)(xx)

## Second derivative estimation
theta<- cmpfun(
  function(X, j, sigma, lambda){
    d<- ncol(X)
    b<- nrow(X)
    
    # constructing h_j
    # mtx<- matrix(rep(X[, j, drop = F], n), ncol = n)
    # dist_mtx<- as.matrix(dist(X))
    # h<- apply(exp(dist_mtx^2/-2*sigma^2)/-4*sigma^2 * ((mtx - t(mtx))^2 - sigma^2), 2, mean)
    h<- apply(phi_func(X, j, sigma)(X), 1, mean)
    
    # get theta_j
    solve(G_mtx(X, sigma) + lambda * diag(b)) %*% h
  }
)

g_func<- function(X, j, sigma, lambda, theta_j = NULL){
  # get theta_j

  cmpfun(
    function(x){
      psi<- psi_func(X, sigma)(x)
      
      theta_j<- theta(X, j, sigma, lambda)
      
      if(class(psi) == "matrix"){
        t(theta_j) %*% t(psi)
      }else{
        t(theta_j) %*% t(t(psi))
      }
    } 
  )
}

# corss validation (hold-out)
# TODO 
# CV

CV_t<- function(X_t, x_t, j, sigma, lambda){
  d<- ncol(X_t)

  # theta_hat<- theta(X_t, j, sigma, lambda)
  theta_hat<- theta(x_t, j, sigma, lambda)
  
  # t(theta_hat) %*% G_mtx(X_t, sigma) %*% theta_hat - (2/n) * sum(g_func(X_t, j, sigma, lambda, theta_hat)(x_t))
    
  t(theta_hat) %*% G_mtx(x_t, sigma) %*% theta_hat - 2 * mean(t(theta_hat) %*% t(phi_func(x_t, j, sigma)(X_t)))
  # t(theta_hat) %*% G_mtx(X_t, sigma) %*% theta_hat - 2 * mean(t(theta_hat) %*% t(phi_func(X_t, j, sigma)(x_t)))
}


CV<- function(X, j, sigma, lambda, t){
  # t: number of partitions
  n<- nrow(X)
  start<- seq(1, n, n/t)
  end<- start + n/t - 1
  
  mean(
    mapply(function(s, e) CV_t(X[-c(s:e),,drop = F], X[s:e,,drop = F], j, sigma, lambda),
           s = start,
           e = end)
  )
}

