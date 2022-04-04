set.seed(2021)

# basic data features
n.depVars <- 2
n.periods <- 3 #2 quarters or 6 months
n.obs <- 5000

# how many montecarlo runs used for the various calls
n.sim <- 5000

#####
# Libraries
#####
if(!require(tsDyn)) { install.packages("tsDyn") }
if(!require(mvtnorm)) { install.packages("mvtnorm") }
if(!require(quantreg)) { install.packages("quantreg") }
if(!require(vars)) { install.packages("vars") }
if(!require(zoo)) { install.packages("zoo") }
library(quantreg)
library(tsDyn)
library(mvtnorm)
library(vars)
library(zoo)

######
# Functions
#####

# this function gives three different options
# the first is a simple geometric mean
# the second is a (known) weighted mean
# and the third is a "random" weights that don't necessarily construct a proxy average
bridge.functions <- function(dat,type){
  if(type == "mean"){
    tmp <- rollmean(dat,k=3,fill = T,align="right")
    #tmp[index(tmp)%%3 != 0] <- NA
  } else if (type == "weighted"){
    tmp <- filter(dat, c(.6,.3,.1), sides = 1)
    # want to make this random
  } else if (type == "random"){
    a <- runif(3,min=0,max=.8)
    tmp <- filter(dat, a, sides = 1)
  }
  return(tmp)
}

######
# State-Space Representation of the Variables
######

rsltMat <- matrix(0,ncol= 2*(n.periods*n.depVars+1),nrow = n.obs)

# matrix of parameters
B2 <- matrix(round(runif(n.depVars+n.depVars*n.depVars*n.periods,-.5,.5),1),
             nrow = n.depVars,ncol=1+n.depVars*n.periods)

# technically at this point should test for stationarity just to be sure
# but for the seed, the first draw seems okay

#####
# Test that VAR recovers B2
#####

# this is to test that the VAR recovers the true parameters appropriately, and there wasn't some mistake
# it doesn't do the BEST, but does okay at nsim = 5000 (errors +/i 1e^-5)
# for(i in 1:n.sim){

  x.star <- data.frame(VAR.sim(B=B2, n=n.obs,lag=3, include="const"))
  names(x.star) <- c("y1","y2")

  if(length(x.star[is.na(x.star)]) == 0 & max(x.star) < Inf & min(x.star) > -Inf){
    var.test <- VAR(x.star,p=3)
    rsltMat[i,1:(n.depVars*n.periods+1)] <- coef(var.test)$y1[,1]-c(B2[1,2:dim(B2)[2]],B2[1,1])
    rsltMat[i,(n.depVars*n.periods+2):(2*(n.depVars*n.periods+1))] <- coef(var.test)$y2[,1]-c(B2[2,2:dim(B2)[2]],B2[2,1])
    
  }
 }

#dim(rsltMat)
#rsltMat <- rsltMat[rowSums(rsltMat) != 0,] # generally doing something wrong here, they're not even close to zero
#dim(rsltMat)
#colMeans(rsltMat)

#####
# BUILD STATE SPACE DATA SET
#####

x.star <- data.frame(VAR.sim(B=B2, n=n.obs,lag=3, include="const"))
names(x.star) <- c("y1","y2")
x.obs <- x.star
x.obs[,2] <- bridge.functions(x.star[,2],type ="mean")

## the response function

beta <- c(2,4)

#x <- runif(n.obs, 0, 4)
#u <- rnorm(n.obs) # homoskedastic errors
#y <- x*beta + u

#plot(x, y)
#abline(lm(y~x), lwd=2, col="blue")
#rq(y~x, seq(0.1,0.9, by=0.1))
#lapply(seq(0.1, 0.9, by=0.1), function(tt) abline(rq(y~x, tau=tt)))

u <- 1/2*(x.star[,1]+x.star[,2])*rnorm(n.obs) # heteroskedastic errors
y.obs <- as.matrix(x.star)%*%beta + u

rq(y.obs~as.matrix(x.star), seq(0.1,0.9, by=0.1))

rq(y.obs~as.matrix(x.obs), seq(0.1,0.9, by=0.1))

