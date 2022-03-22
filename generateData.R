## generate simulation data

generateData <- function(n, beta=c(1,2,-1), prop.out=0.1, prop.cens=0.2,
                         cens.mode=c("uniform","non-outlier only", "censored outlier"),
                         x.con=F,
                         h0=2, k.shape=1, lp.con=T, 
                         y.con=F, gamma=rep(0,n), seed=NULL){
  num.out <- n * prop.out
  num.binary <- 2
  p <- length(beta)
  ## generate x matrix
  mu.x <- rep(0, p)
  Sigma.x <- matrix(rep(0.1,p*p),p,p)+diag(0.9,p)
  X <- mvrnorm(n, mu=mu.x, Sigma = Sigma.x)
  ## generate binary variable
  index.binary <- which((1:sum(beta!=0))%%2==0)
  X[,index.binary] <- ifelse(X[,index.binary]>0,1,0)
  if(x.con){ # contamination on x direction
    X[(n-num.out+1):n, ] <- matrix(rchisq(num.out*p,df=5),nrow=num.out,ncol=p)
  }
  if(!lp.con){# no contamination in lp
    gamma=rep(0,n)
  }
  if(lp.con & all.equal(gamma, rep(0,n))& num.out>0){

    # k.gamma <-c(rep(2, floor(num.out*0.5)), rep(-2,ceiling(num.out*0.5)))
    k.gamma <-c(rep(5, floor(num.out*0.5)), rep(-5,ceiling(num.out*0.5)))
    # k.gamma <-c(rep(8, floor(num.out*0.5)), rep(-8,ceiling(num.out*0.5)))
    gamma <- c(rep(0, n-num.out), k.gamma)
  }
  ## generate failure time
  r.u <- runif(n)
  ftime <- (-log(r.u)/(h0^k.shape*exp(X %*% beta + gamma )))^(1/k.shape)
  if(y.con){ # contamination on failure time
    ftime[(n-num.out+1):n] <- ftime[(n-num.out+1):n]+10
  }
  if(cens.mode=="uniform"|num.out==0|prop.cens==0){
    obs <- getObservedTime(ftime, prop.cens)
    survtime <- obs$survtime
    delta <- obs$delta
    ncens=obs$ncens
  } else if (cens.mode=="non-outlier only"){
    obs <- getObservedTime(ftime[1:(n-num.out)], prop.cens = n/(n-num.out) * prop.cens)
    survtime <- c(obs$survtime, ftime[(n-num.out+1):n])
    delta <- c(obs$delta, rep(1, num.out))
    ncens <- obs$ncens
  } else if(cens.mode=="censored outlier"){
    num.cens.out <- round(num.out * prop.cens/2)
    obs <- getObservedTime(ftime[1:(n-num.out)], prop.cens=(n*prop.cens-num.cens.out)/(n-num.out))
    survtime <- c(obs$survtime, ftime[(n-num.out+1):n])
    delta <- c(obs$delta, rep(1, num.out))
    if(num.cens.out!=0){
      ind.cens.out <- sample(which(gamma<0), num.cens.out)
      delta[ind.cens.out] <- 0
    }
    ncens <- obs$ncens + num.cens.out
  }
  
  y <- Surv(survtime, as.matrix(delta,nrow=n))
  list(y=y, X=X, survtime=survtime, delta=delta, L=obs$L, ncens=ncens)
}


## generate and save data to file
saveDataByScenario <- function(m = 100, n, beta, prop.out=0.1, prop.cens=0.2, 
                               seed = NULL, scenario = c("NU1", "NU2", "NN1", "NN2", "NC1", "NC2"), 
                               path = "../data/", ...) {
  # browser()
  ns <- length(scenario)
  for (i in 1:ns) {
    dfile <- paste0(path, scenario[i],"-", n, "X", length(beta),
                   "(",prop.cens*100, "C", prop.out*100, "O).rda")
    x.con<- !(substr(scenario[i], 1, 1) == "N")
    cens.code <- substr(scenario[i], 2, 2)
    cens.mode <- ifelse(cens.code=="U", "uniform", ifelse(cens.code=="N", "non-outlier only", "censored outlier"))
    base.type =substr(scenario[i], 3, 3)
    if(base.type=="1"){# exponential
      h0 <- 2
      k.shape=1
    }else if (base.type=="2"){ ## weibull
      h0 <- 2
      k.shape <- 0.5 
    }
    datas <- array(list(), m)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    for (j in 1:m) {
      datas[[j]] <- generateData(n, beta = beta, x.con=x.con, prop.out=prop.out, 
                                 prop.cens=prop.cens, seed = NULL,h0=h0, k.shape=k.shape,cens.mode = cens.mode, ... )
    }
    save(datas, file = dfile)
  }
}

loadDataByScenario <- function(n, beta, prop.out=0.1, prop.cens=0.2, 
                              scenario = c("NU1", "NU2", "NN1", "NN2", "NC1", "NC2"), path = "../data/") {
  scenario <- match.arg(scenario)
  fname <- paste0(scenario,"-", n, "X", length(beta),
                  "(",prop.cens*100, "C", prop.out*100, "O)")
  dfile <- paste0(path, fname, ".rda")
  load(dfile)
  return(list(datas=datas, name=fname))
}

getObservedTime <- function(ftime, prop.cens, L.init=10, tol=0.01){
  n <- length(ftime)
  if(prop.cens==0){
    return(list(survtime=ftime, delta=rep(1,n), L=NA, necens=0))
  }

  L.low <- min(ftime)
  L.up <- max(ftime)
  if(L.low>L.init | L.up < L.init)
    L.init <- min((L.low + L.up)/2, L.low+L.init)
  L.curr <- L.init
  ncens.low <- n*(prop.cens-tol)
  ncens.up <- n * (prop.cens+tol)
  while(TRUE){
    # get current number of censoring
    ctime <- pmin(runif(n, 0, L.curr+2), L.curr)
    delta <- ifelse(ftime<=ctime, 1, 0)
    ncens.curr <- n- sum(delta)
    if(ncens.low <= ncens.curr & ncens.curr <=ncens.up) break
    if(ncens.curr< ncens.low){
      L.up <- L.curr
      L.curr <- max((L.up+L.low)/2,L.curr/2)
    }else if(ncens.up<ncens.curr){
      L.low <- L.curr
      L.curr <- min((L.low+L.up)/2, 2*L.curr)
    }
  }
  survtime<- pmin(ftime, ctime)
  return(list(survtime=survtime, delta=delta, L=L.curr, ncens=ncens.curr))
}