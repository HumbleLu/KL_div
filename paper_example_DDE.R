source("https://raw.githubusercontent.com/HumbleLu/KL_div/master/DDE.R")

## true value
tv_2nd<- function(x, sigma){
  dnorm(x, 0, sigma) * (x^2 - sigma^2) / sigma^4
}

## generate data
X<- matrix(rnorm(500, 0, 1), 500, 1)

## get optimal parameters
sigma_cand<- 10^seq(-.3, 1, length.out = 9)
lambda_cand<- 10^seq(-1, 1, length.out = 9)
opt_par<- cross_val(X, sigma_cand, lambda_cand, 1)

## plot result
x<- seq(-4, 4, .01)
plot(x, tv_2nd(x, 1), type = "l") ## true values
g_est<- g_func(X, 1, opt_par$sigma, opt_par$lambda)
lines(x, g_est(t(t(x))), lty = 2) ## estimations