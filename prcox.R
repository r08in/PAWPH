lik.surv <- function(t.train, delta.train,lp.train, t.test, delta.test,lp.test){
  lik.mat<- -log(sapply(1:ncol(lp.train), function(j){
    if(sum(!is.na(lp.train[,j]))==0){## all are NA
      rep(0, length(t.test))
    }else{
      cum.basehaz <- basehaz.gbm(t.train, delta.train, f.x=lp.train[,j], t.eval=t.test,
                                 cumulative = TRUE)
      val.basehaz <- basehaz.gbm(t.train, delta.train, f.x=lp.train[,j], t.eval=t.test,
                                 smooth = TRUE, cumulative = FALSE)
      St.base <- exp(-cum.basehaz)
      St <- St.base^(exp(lp.test[,j]))
      #liklihood
      if_else(delta.test==1, St * val.basehaz * exp(lp.test[,j]), St)
    }
  }))
  lik.mat[lik.mat==Inf] <- NA
  lik.mat
}

prcox <- function(y, X, penalty.factor=c(rep(0,p),  rep(1,n)), lambda.min=0.05, lambda, ...){
  p <- ncol(X)
  n <- nrow(X)
  X.aug <- cbind(X, diag(n))
  # lasso fit
  if(missing(lambda)){
    ncvsurv(X.aug, y, penalty="lasso", lambda.min = lambda.min, 
            penalty.factor = penalty.factor, returnX = FALSE, ...)
  }else{
    ncvsurv(X.aug, y, penalty="lasso",  penalty.factor = penalty.factor,
            lambda=lambda, returnX = FALSE, ...)
  }
}

cvf.prcox <- function(i,XX, y, fold, cv.args){
  # browser()
  p <- ncol(XX)
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i,]
  cv.args$permu.weight <- cv.args$permu.weight[fold!=i]
  cv.args$penalty.factor <- cv.args$penalty.factor[c(1:p, p+which(fold!=i))]
  fit.i <- do.call("prcox", cv.args)
  prop.out <- apply(fit.i$beta,2, function(beta) {
    sum(beta[-(1:p)]!=0)/length(cv.args$y)
  })
  X.aug <- cbind(cv.args$X, diag(nrow(cv.args$X)))
  lp.train <- X.aug %*% fit.i$beta
  test.X <- XX[fold==i, , drop=FALSE]
  test.y <- y[fold==i,]
  lp.test <- test.X %*% fit.i$beta[1:p,]
  lik <- lik.surv(t.train=cv.args$y[,1], delta.train=cv.args$y[,2], lp.train=lp.train,
                  t.test=test.y[,1], delta.test=test.y[,2], lp.test=lp.test)
  colnames(lik) <- fit.i$lambda
  return(list(lik=lik, prop.out=prop.out))
}
cv.prcox <- function(y,X, nfolds=10, seed, 
                     penalty.factor=c(rep(0,p),  rep(1,n)),pout.max=NULL, dfmax=100,
                     ...){
  n <- nrow(X)
  p <- ncol(X)
  fit <- prcox(y, X, penalty.factor = penalty.factor,...)
  
  # begin cross-validaton
  fail <- y[,2]
  # set up folds
  if (!missing(seed)) set.seed(seed)
  ind1 <- which(fail==1)
  ind0 <- which(fail==0)
  n1 <- length(ind1)
  n0 <- length(ind0)
  fold1 <- 1:n1 %% nfolds
  fold0 <- (n1 + 1:n0) %% nfolds
  fold1[fold1==0] <- nfolds
  fold0[fold0==0] <- nfolds
  fold <- integer(n)
  fold[fail==1] <- sample(fold1)
  fold[fail==0] <- sample(fold0)
  #likelihood on validation set
  lik.val<- matrix(Inf, nrow=n, ncol=length(fit$lambda))
  colnames(lik.val) <- fit$lambda
  #ags
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$penalty.factor <- penalty.factor
  for(i in 1:nfolds){
    res <- cvf.prcox(i, X, y, fold, cv.args)
    lik.val[fold==i,colnames(res$lik)] <- res$lik
  }
  cve <- apply(lik.val, 2, function(x) {
    xx <- x[x<quantile(x, 0.8, na.rm=T)]
    mean(xx, na.rm = T)
  })
  if(!is.null(pout.max)){
    prop.out <- apply(fit$beta[-(1:p),], 2, function(b) sum(b!=0))/n
    lam.ind <- which(prop.out<=pout.max) 
    cve[-lam.ind] <- Inf
  }
  # exclude tuning parameter producing NA
  cve[which(apply(fit$beta, 2, anyNA))] <- Inf
  lam <- fit$lambda[which.min(cve)]
  betaHat <- coef(fit, lambda=lam)
  gammaHat <- betaHat[(p+1):(n+p)]
  betaHat <- betaHat[1:p]
  list(betaHat=betaHat, gammaHat=gammaHat, fit=fit, cve=cve, cve.min=min(cve, na.rm = T),
       opt.lam=lam, lik.val=lik.val, lambda=fit$lambda)
}
# 

cv.prcox.adap <- function(y, X, seed, alpha=1, penalty.factor=rep(1, p+n), lambda.min=0.05,
                          pout.max=NULL, ...){
  p <- ncol(X)
  n <- nrow(X)
  res_list <- array(list(), length(alpha))
  for(i in 1:length(alpha)){
    res0 <- cv.prcox(y, X,seed = seed,
                     penalty.factor=c(rep(1-alpha[i],p),rep(alpha[i],n))*penalty.factor,
                     lambda.min=lambda.min, pout.max=pout.max, ...)
    gammaHat <- res0$gammaHat
    betaHat <- res0$betaHat
    
    res_list[[i]] <- cv.prcox(y, X, 
                              seed=seed,
                              penalty.factor=c((1-alpha[i])/pmax(abs(betaHat), 1e-5),
                                               alpha[i]/pmax(ifelse(gammaHat>0, gammaHat, -gammaHat),1e-5))*penalty.factor,
                              lambda.min=lambda.min, pout.max=pout.max, ...)
    
  }
  cve <- sapply(res_list, getElement,"cve.min")
  i.alpha <- which.min(cve)
  res <- res_list[[i.alpha]]
  res$opt.alpha <- alpha[i.alpha]
  return(res)
}


prcoxreg <- function(y, X, seed, alpha=1, cond.CI=FALSE,pout.max=0.2, permu.weight=rep(1, n), ...){
  p <- ncol(X)
  n <- nrow(X)
  delta <- y[,2]
  is.vs <- !(all.equal(alpha,1)==T)
  is.cens <- sum(delta)!=n
  prop.cens <- sum(1-delta)/n
  pout.max2 <- NULL
  if(!is.null(pout.max)){
    pout.max2 <- pout.max+prop.cens
  }
  if(is.vs){
    res.vs <- cv.prcox.adap(y, X, seed, alpha=alpha, ...)
    betaHat.vs <- res.vs$betaHat
    # all as event for outlier detection
    y[,2] <- rep(1, n)
    if(sum(betaHat.vs!=0)>0){
      res <- cv.prcox.adap(y, as.matrix(X[,which(betaHat.vs!=0)]), seed, 
                           alpha=1, pout.max=pout.max2)
      res$betaHat <- betaHat.vs
      res$opt.alpha <- res.vs$opt.alpha
    }else{
      res <- cv.prcox.adap(y, X, seed, alpha=alpha, pout.max=pout.max2,...)
    }
    
  }else{
    y[,2] <- rep(1, n)
    res<- cv.prcox.adap(y, X, seed, alpha=alpha, pout.max=pout.max2,permu.weight=permu.weight)
  }
  ind.cens <- which(delta==0)
  res$gammaHat[ind.cens] <- ifelse(res$gammaHat[ind.cens]>-2, 0, res$gammaHat[ind.cens])
  y[,2] <- delta
  ## refit
  if(sum(res$betaHat!=0)>0){
    fit <- ncvsurv(as.matrix(X[,which(res$betaHat!=0)]), y, penalty="lasso",intercept =res$gammaHat, lambda=0, returnX = FALSE,permu.weight = permu.weight)
    betaHat_re <- rep(0,p)
    betaHat_re[which(res$betaHat!=0)] <- as.numeric(fit$beta)
    res$betaHat_re <- betaHat_re
  }else{
    res$betaHat_re <- res$betaHat
  }
  names(res$betaHat_re) <- names(res$betaHat)
  res$gammaHat <- ifelse(abs(res$gammaHat)<0.5, 0, res$gammaHat)
  
  ## CI
  if(cond.CI){
    se.res <- SE.boot(y, X, pout.max=pout.max, alpha=ifelse(is.vs, res$opt.alpha, 1), seed=seed,...)
    res$CI <- sapply(1:ncol(X), function(i) getCI(se.res$betaHat[,i], b=res$betaHat_re[i]))
    rownames(res$CI) <- c("SE", "lower_n", "upper_n", "lower_q", "upper_q", "z", "pvalue")
    colnames(res$CI) <- names(res$betaHat)
    res$se.res <- se.res
  }
  return(res)
}

SE.boot <- function(orig_y, orig_X, B=100, seed=1234, ...){
  set.seed(seed)
  n <- nrow(orig_X)
  p <- ncol(orig_X)
  boot.args <- list(...)
  cl <- makeCluster(16)
  registerDoParallel(cl)
  res <-  foreach(ind=replicate(B, sample(n,n,replace = T), simplify = FALSE),
                  .packages = c("ncvreg2","gbm", "dplyr","survival", "glmnet"))%dopar%{
                    source("prcox.R")
                    if(!is.null(boot.args$penalty.factor)){
                      boot.args$penalty.factor <- boot.args$penalty.factor[c(1:p, p+ind)]
                    }
                    boot.args$X <- orig_X[ind,]
                    boot.args$y <- orig_y[ind,]
                    boot.args$seed <- seed
                    do.call("prcoxreg", boot.args)
                  }
  betaHat <- t(sapply(res, getElement, "betaHat_re"))
  stopCluster(cl)
  registerDoSEQ()
  return(list(betaHat=betaHat, res=res))
}


getCI <- function(betas, b, alpha=0.05){
  se <- sd(betas)
  c(SE=se, 
    lower_n=b-qnorm(1-alpha/2)*se,
    upper_n=b+qnorm(1-alpha/2)*se,
    lower_q=quantile(betas, probs = alpha/2),
    upper_q=quantile(betas, probs = 1-alpha/2),
    z.stat= b/se,
    pvalue=pnorm(abs(b/se), lower.tail = F)*2)
  
}