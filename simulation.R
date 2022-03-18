library(survival)
library(MASS)
library(glmnet)
library(gbm)
library(dplyr)
library(doParallel)
library(foreach)
library(coxrobust)
source("generateData.R")
source('prcox.R')

simulation <- function(m, n, beta, prop.out=0.1, prop.cens=0.2, 
                       method=c("n.cox", "r.cox", "o.cox","op.cox","prcox", "np.cox", "coxr"),
                       scenario = c("NU1", "NU2", "NN1", "NN2", "NC1", "NC2"), seed=2021, 
                       res_rda=NULL,...){
  cl <- makeCluster(parallel::detectCores())
  registerDoParallel(cl)
  #registerDoSEQ()
  p <- length(beta)
  ind_beta <- which(beta!=0)
  numout <- n * prop.out
  if(numout==0){
    ind.o<- NULL
  }else{
    ind.o <- (n-numout+1): n
  }
  
  ns <- length(scenario)
  nm <- length(method) # number of methods
  nres <- array(list(), ns)
  betaHat <- gammaHat <- alpha <- NULL
  betaHat_re <-  matrix(0, nrow=p, ncol=m)
  
  ## for each scenario
  for(ks in 1:ns){
    #load data
    data_list <-  loadDataByScenario(n, beta, prop.out = prop.out, prop.cens = prop.cens, scenario = scenario[ks])
    datas <- data_list$datas
    datas <- datas[1:m]
    X_list <- lapply(datas, getElement, "X")
    y_list <- lapply(datas, getElement,"y")
    deltaMat <- sapply(datas, getElement, "delta")
    nres[[ks]] <- array(list(), nm)
    names(nres)[ks] <- data_list$name
    for(km in 1:nm){
      if(!is.null(res_rda)){
        betaHat <- res_rda[[km]]$betaHat
        betaHat_re <- res_rda[[km]]$betaHat_re
        gammaHat <- res_rda[[km]]$gammaHat
        alpha <- res_rda[[km]]$alpha
      }else if(method[km]=="n.cox"){
        # naive cox
        betaHat <- foreach(y=y_list, X=X_list, .packages = "survival") %dopar% {
          fit <- coxph(y~X)
          summary(fit)$coefficients[,"coef"]
        }
        betaHat <- sapply(betaHat, function(a) a)
        gammaHat <- NULL
      }else if(method[km]=="coxr"){
        # robust cox
        betaHat <- foreach(y=y_list, X=X_list, .packages = c("survival", "coxrobust")) %dopar% {
          fit <- coxr(y~X, trunc = 0.95)
          fit$coefficients
        }
        #browser()
        betaHat <- sapply(betaHat, function(a) a)
        gammaHat <- NULL
      }else if(method[km]=="np.cox"){
        betaHat <- foreach(y=y_list, X=X_list, .packages = c("survival", "ncvreg2")) %dopar% {
          fit <- cv.ncvsurv(X, y, penalty="lasso",returnX = FALSE, seed = seed)
          fit <- cv.ncvsurv(X, y, penalty="lasso",penalty.factor = 1/pmax(abs(coef(fit)), 1e-5), returnX = FALSE,
                            seed=seed)
          coef(fit)
        }
        betaHat <- sapply(betaHat, function(a) a)
        gammaHat <- NULL
      }else if(method[km]=="r.cox"){
        # estimator by removing outlier
        betaHat <- foreach(y=y_list, X=X_list, .packages = "survival") %dopar% {
          summary(coxph(y[1:(n-numout)]~X[1:(n-numout),]))$coefficients[,"coef"]
        }
        betaHat <- sapply(betaHat, function(a) a)
        gammaHat <- NULL
      }else if(method[km]=="rp.cox"){
        yy <- y_list[[1]][1:(n-numout)]
        XX <- X_list[[1]][1:(n-numout),]
        ## alasso
        fit <- cv.ncvsurv(XX, yy, seed=seed, penalty="lasso",returnX = FALSE)
        fit <- cv.ncvsurv(XX, yy, seed=seed, penalty="lasso",penalty.factor = 1/pmax(abs(coef(fit)), 1e-5), returnX = T)
        beta.alasso <- coef(fit)
        #local_mfdr(fit$fit, lambda=fit$lambda.min)
        SE.alasso <- SE.alasso(yy, XX, B=100)
        # estimator by removing outlier
        betaHat <- foreach(y=y_list, X=X_list, .packages = c("survival", "ncvreg2")) %dopar% {
          fit <- cv.ncvsurv(X[1:(n-numout),], y[1:(n-numout)], penalty="lasso",returnX = FALSE, seed = seed)
          fit <- cv.ncvsurv(X[1:(n-numout),], y[1:(n-numout)], penalty="lasso",
                            penalty.factor = 1/pmax(abs(coef(fit)), 1e-5), returnX = FALSE, seed=seed)
          coef(fit)
        }
        betaHat <- sapply(betaHat, function(a) a)
        gammaHat <- NULL
        
      }else if(method[km]=="o.cox"){
        # oracle estimator
        betaHat <- foreach(y=y_list, X=X_list, .packages = "survival") %dopar% {
          X.aug <- cbind(X, diag(nrow(X)))
          beta_tmp <- summary(coxph(y~X.aug[, c(1:p, p+ind.o)]))$coefficients[,"coef"]
          beta_tmp[1:p]
        }
        betaHat <- sapply(betaHat, function(a) a)
        gammaHat <- NULL
      }else if(method[km]=="op.cox"){
        # penalized oracle estimator
        if(numout==0)
          betaHat=NULL
        else{
          betaHat <- foreach(y=y_list, X=X_list, .packages = "ncvreg2") %dopar% {
            X.aug <- cbind(X, diag(nrow(X)))
            fit <- ncvsurv(X.aug[, c(1:p, p+ind.o)], y, penalty="lasso", penalty.factor = c(rep(0,p), rep(1,numout)),
                           seed=seed, returnX = FALSE)
            lam <- fit$lambda[which.min(BIC(fit))]
            betaHat.op <- coef(fit, lambda=lam)
            betaHat.op[1:p]
          }
          betaHat <- sapply(betaHat, function(a) a)
          gammaHat <- NULL
        }
        
      }else if(method[km]=="prcox"){
        ## prcox
        res <- foreach(y=y_list, X=X_list, .packages = c("ncvreg2","gbm", "dplyr","survival", "glmnet")) %dopar%{
          source("prcox.R")
          prcoxreg(y=y, X=X, seed = seed, ...)
        }
        # browser()
        # for(i in 1:m){
        #   res <- prcoxreg(y=y_list[[i]], X=X_list[[i]], seed=seed, cond.CI=F,...)
        # }
        betaHat <- sapply(res, getElement, "betaHat")
        betaHat_re <- sapply(res, getElement, "betaHat_re")
        gammaHat <- sapply(res, getElement, "gammaHat")
        alpha <- sapply(res, getElement, "opt.alpha")
      }
      ## report metrics
      names(nres[[ks]])[km] <- method[km]
      nres[[ks]][[km]] <- c(list(scenario=scenario[ks], method=method[km],
                                 betaHat=betaHat, betaHat_re=betaHat_re,
                                 gammaHat=gammaHat, alpha=alpha),
                            reportMetrics(betaHat, beta, gammaHat, ind.o, deltaMat),
                            re=reportMetrics(betaHat_re, beta, gammaHat, ind.o, deltaMat))
    }
  }

  
  
  
  return(nres)
}


reportMetrics <- function(betaHatMat, beta, gammaHatMat, ind.o, deltaMat){
  if(is.null(betaHatMat))
    return(NA)
  numout <- length(ind.o)
  mse.total <- mean(apply(betaHatMat, 2, function(betaHat) sum((beta - betaHat)^2)))
  mse.total.sd <- sd(apply(betaHatMat, 2, function(betaHat) sum((beta - betaHat)^2))) 
  abias.total <- mean(apply(betaHatMat, 2, function(betaHat) sum(abs(beta - betaHat))))
  abias.total.sd <- sd(apply(betaHatMat, 2, function(betaHat) sum(abs(beta - betaHat))))
  bias.indi <- apply((betaHatMat-beta), 1, mean)
  mse.indi<- apply((betaHatMat-beta)^2, 1, mean)
  se<- apply(betaHatMat, 1, sd)
  est <- apply(betaHatMat, 1, mean)
  ## variable selection metric
  true_set <- which(beta!=0)
  pnum <- length(true_set)
  betaMetricMat <- apply(betaHatMat, 2, function(betaHat){
    active_set <- which(betaHat!=0)
    msize <- length(active_set)
    common_size <- length(intersect(true_set,active_set))
    cfr <- ifelse(common_size == pnum & msize== pnum, 1, 0)  # correctly fit
    pdr <- common_size/pnum 
    fdr <- ifelse(msize==0,0,(msize - common_size)/msize)
    fpr <- (msize-common_size)/(p-pnum)
    fnr<- (pnum-common_size)/pnum
    c(msize, cfr, pdr, fdr, fpr, fnr)
  })
  betaMetric <- c(apply(betaMetricMat, 1, mean), apply(betaMetricMat, 1, sd))
  names(betaMetric) <-c("msize", "cfr", "pdr", "fdr", "fpr", "fnr", 
                        "msize.sd", "cfr.sd", "pdr.sd", "fdr.sd", "fpr.sd", "fnr.sd")
  metric_list <- list(est=est, se=se, mse.indi=mse.indi, bias.indi=bias.indi, mse.total=mse.total, mse.total.sd=mse.total.sd,
                      abias.total=abias.total, abias.total.sd=abias.total.sd,
                      select=betaMetric)
  if(!is.null(gammaHatMat)){
    n <- nrow(gammaHatMat)
    m <- ncol(gammaHatMat)
    vout <- rep(FALSE,n)
    vout[ind.o] <- TRUE
    maskprob <- swamprob <- mask.event <- mask.cens<- R <- rep(NA, m)
    maskprobPos <- maskprobNeg <- rep(NA,m)
    num.oevent <- num.ocens<- rep(NA, m)
    for(j in 1:m){
      ind.oh <- which(abs(gammaHatMat[,j])>0)
      swamprob[j] <- (length(ind.oh) - length(intersect(ind.o, ind.oh)))/(n-numout)
      R[j] <- sum((gammaHatMat[,j]!=0)==vout)/n
      if(numout!=0){
        maskprob[j] <- 1-length(intersect(ind.o, ind.oh))/numout
        maskprobPos[j] <- 1-length(intersect(ind.o[1:floor(numout/2)], ind.oh))/floor(numout/2)
        maskprobNeg[j] <- 1-length(intersect(ind.o[-(1:floor(numout/2))], ind.oh))/ceiling(numout/2)
        ind.oevent <- intersect(which(deltaMat[,j]==1), ind.o)
        ind.ocens <- intersect(which(deltaMat[,j]==0), ind.o)
        num.oevent[j] <- length(ind.oevent)
        num.ocens[j] <- length(ind.ocens)
        if(num.oevent[j]!=0)
          mask.event[j] <- 1-length(intersect(ind.oevent, ind.oh))/length(ind.oevent)
        if(num.ocens[j]!=0)
          mask.cens[j] <- 1-length(intersect(ind.ocens, ind.oh))/length(ind.ocens)
      }
      
    }
    metric_list <- c(metric_list, 
                     list(maskprob=round(mean(maskprob), digits = 2), 
                          maskprobPos=round(mean(maskprobPos), digits = 2),
                          maskprobNeg=round(mean(maskprobNeg), digits = 2),
                          mask.event=round(mean(mask.event), digits = 2), 
                          num.oevent=round(mean(num.oevent), digits = 2),
                          mask.cens=round(mean(mask.cens), digits = 2), 
                          num.ocens=round(mean(num.ocens), digits = 2),
                          swamprob=round(mean(swamprob), digits = 2), 
                          maskprob.sd=round(sd(maskprob), digits = 2), 
                          maskprobPos.sd=round(sd(maskprobPos), digits = 2),
                          maskprobNeg.sd=round(sd(maskprobNeg), digits = 2),
                          mask.event.sd=round(sd(mask.event), digits = 2), 
                          num.oevent.sd=round(sd(num.oevent), digits = 2),
                          mask.cens.sd=round(sd(mask.cens), digits = 2), 
                          num.ocens.sd=round(sd(num.ocens), digits = 2),
                          swamprob.sd=round(sd(swamprob), digits = 2)))
  }
  return(metric_list)
}
