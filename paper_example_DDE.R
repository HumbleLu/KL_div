source("https://raw.githubusercontent.com/HumbleLu/KL_div/master/DDE.R")

## cross validation, model selection
cross_val<-function(sigma_cand, lambda_cand){
  pars<- expand.grid(x = sigma_cand, y = lambda_cand)
  
  sigma_opt<- sigma_cand[1]
  lambda_opt<- lambda_cand[1]
  cv_opt<- CV(X, 1, sigma_opt, lambda_opt, 10)
  print(paste("sigma:", sigma_cand[1], ",lambda:", lambda_cand[1], ",CV: ", cv_opt))

  for(i in 2: nrow(pars)){
    gc()
    cv_cur<- CV(X, 1, pars[i,1], pars[i,2], 10)
    print(paste("sigma: ", pars[i,1], ",lambda: ", pars[i,2], ",CV: ", cv_cur))
    if(cv_cur< cv_opt){
      cv_opt<- cv_cur
      sigma_opt<- pars[i,1]
      lambda_opt<- pars[i,2]
    }
  }
  
  list(sigma = sigma_opt, lambda = lambda_opt)
}

## true value
tv_2nd<- function(x, sigma){
  dnorm(x, 0, sigma) * (x^2 - sigma^2) / sigma^4
}

## generate data
X<- matrix(rnorm(500, 0, 1), 500, 1)

## get optimal parameters
sigma_cand<- 10^seq(-.3, 1, length.out = 9)
lambda_cand<- 10^seq(-1, 1, length.out = 9)
opt_par<- cross_val(sigma_cand, lambda_cand)

## plot result
x<- seq(-4, 4, .01)
plot(x, tv_2nd(x, 1), type = "l") ## true values
g_est<- g_func(X, 1, opt_par$sigma, opt_par$lambda)
lines(x, g_est(t(t(x))), lty = 2) ## estimations