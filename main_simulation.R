rm(list=ls())
library(survival)
library(MASS)
library(ncvreg2)
library(gbm)
library(dplyr)
library(glmnet)
source('generateData.R')
source('simulation.R')
#################################################################
##                        Data Settings                        ##
#################################################################
seed <- 2021
set.seed(seed)
m <- 100
p <- 8
beta <- c(c(1,2,-1), rep(0, p-3))
#pnum <- 6
##-2.359632  2.898054  3.922570 -2.261744  2.915566  2.013772
#beta <- c(ifelse(runif(pnum)>=0.5,1,-1)* (2+abs(rnorm(pnum))), rep(0, p-pnum))
#method=c("n.cox","r.cox", "coxr", "prcox" )
method <- c("np.cox","rp.cox", "prcox")
scenario <- c("NC2")
n.array <- c(100, 300, 500)
prop.out.array <- c(0,0.05, 0.10, 0.20)
prop.cens.array <- c(0, 0.2, 0.4)

# 
# ##################################################################
##                   Generate Simulation data                   ##
##################################################################
#saveDataByScenario(m=2, n=100, beta = c(1,2,-1), prop.out = 0.1, prop.cens =0.2, scenario = "NC2", seed = 2021)
for(n in n.array){
  for(prop.out in prop.out.array){
    for(prop.cens in prop.cens.array){
      saveDataByScenario(m=m, n=n, beta = beta, prop.out = prop.out, prop.cens = prop.cens, scenario = scenario, seed = seed)
    }
  }
}


##################################################################
##                          Simulation                          ##
##################################################################
for(scen in scenario){
  res <- list()
  for(n in n.array){
    for(prop.out in prop.out.array){
      for(prop.cens in prop.cens.array){
        print(system.time(res <- c(res, simulation(m=m, n, beta, prop.out, prop.cens, method, scenario=as.vector(scen), seed = seed
                                                   #,alpha=1))))
                                              ,alpha= seq(0.9, 0.1, length.out=10)))))
      }
    }
  }
  file_name <- paste0("res_", scen, ".RData")
  save(res, file = file_name)
}
