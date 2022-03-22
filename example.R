##Simple Example
source('prcox.R')
library(MASS)
library(survival)
library(MASS)
library(ncvreg2)
library(gbm)
library(dplyr)
seed <- 123
set.seed(seed)
n <- 300
p <- 8
beta <- c(1, -2, 1, rep(0, p-3))
## generate x matrix
mu.x <- rep(0, p)
Sigma.x <- matrix(rep(0.1,p*p),p,p)+diag(0.9,p)
X <- mvrnorm(n, mu=mu.x, Sigma = Sigma.x)
## generate failure time with long-survival outliers
num.out <- n*0.05
gamma <- c(rep(0, n-num.out),  rep(-5,ceiling(num.out)))
k.shape <- 1
h0 <- 2
r.u <- runif(n)
ftime <- (-log(r.u)/(h0^k.shape*exp(X %*% beta + gamma )))^(1/k.shape)
delta <- c(rep(0, n*0.2), rep(1,n*0.8 ))
y <- Surv(ftime, delta) 

# fit 
res <- prcoxreg(y,X, seed=seed, alpha= 0.5)
#vanilla PAWPH estimator
res$betaHat

# PAWPH estimator
res$betaHat_re

# estiamted log(w)
res$gammaHat
