library(MASS)

## G matrix
G_mtx<- function(X, sigma){
  d<- nrow(X)
  (pi * sigma^2)^(d/2) * exp(-1 * as.matrix(dist(X))^2/4*sigma^2)
}


## psi: Guassian model
psi_gauss<- function(x, X, sigma){
  exp(apply(X, 1, function(y) sum((x - y)^2))/(-2 * sigma^2))
}


## psi_j: here the trace of Hessian matrix
psi_2<- function(X, j, sigam){
  function(x){
    psi_gauss(x, X, 2) * ((x[j] - X[,j])^2 - sigma^2)/(-sigma^4)
  }
}

# test
# X<- matrix(rnorm(100), 5, 20)
# G_mtx(X)
# x<- rnorm(5)
# psi_2(X, 5, x)(x)

## Second derivative estimation
g_func<- function(X, j, sigma, lambda){
  n<- nrow(X)
  h<- apply(apply(X, 1, function(x) psi_2(X, j, sigma)(x)), 1, mean)
  theta_hat<- ginv(G_mtx(X, sigma) + lambda * diag(n)) %*% h
  function(x) (t(theta_hat) %*% psi_gauss(x, X, sigma))[1,1]
}

# corss validation (hold-out)
# TODO 
# CV

X_t<- X[-c(1:5),]
x_t<- X[1:5,]

CV_t<- function(X_t, x_t, j, sigma, lambda){
  d<- ncol(X_t)
  n<- nrow(X_t)
  
  pairs<- as.data.frame(expand.grid(x = 1:n, y = 1:n))
  pairs[pairs$x != pairs$y,]
  
  h<- apply(apply(X_t, 1, function(x) psi_2(X_t, j, sigma)(x)), 1, mean)
  theta_hat<- ginv(G_mtx(X_t, sigma) + lambda * diag(n)) %*% h
  
  sum((theta_hat)^2 * (sqrt(pi)*sigma)^d) +
    sum(
      mapply(
        function(i, j) 2 * theta_hat[i] * theta_hat[j] * (sqrt(pi)*sigma)^d * exp(sum((X[i,,drop = F] - X[j,,drop = F])^2/(-4*sigma^2))),
        i = pairs$x,
        j = pairs$y
      )  
    ) + 
    sum(sapply(x_t, function(x) g_func(X_t, j, sigma, lambda)(x)))
}

CV<- function(X, j, sigma, lambda, t){
  n<- nrow(X)
  start<- seq(1, n, t)
  end<- seq(5, n, t)
  
  mean(
    mapply(function(s, e) CV_t(X[-c(s:e),,drop = F], X[s:e,,drop = F], j, sigma, lambda),
           s = start,
           e = end)
  )
}



## cross validation
X<- matrix(rnorm(5000*5), 5000, 5)

sigma_cand<- 10^seq(-3, 1, length.out = 10)
lambda_cand<- 10^seq(-1, 1, length.out = 10)
pars<- expand.grid(x = sigma_cand, y = lambda_cand)

sigma_opt<- sigma_cand[1]
lambda_opt<- lambda_cand[1]
cv_opt<- CV(X, 1, sigma_opt, lambda_opt, 5)

for(i in 2: nrow(pars)){
  cv_cur<- CV(X, 1, pars[i,1], pars[i,2], 5)
  if(cv_cur< cv_opt){
    cv_opt<- cv_cur
    sigma_opt<- sigma_cand[i]
    lambda_opt<- lambda_cand[i]
  }
}

