# Date: 20250618
# Author: Fanjing Guo
# Title: Functions of common normalization methods

# PACKAGE ####
library(preprocessCore)

## Constant Sum Normalization (CSN) ###
CSN_normalize <- function(dat) {
  p <- ncol(dat)
  
  time.start <- Sys.time()
  
  samp.sum <- rowSums(dat)
  coef.csn <- 1 / samp.sum
  dat.norm <- diag(coef.csn) %*% dat
  
  # record running time
  time.end <- Sys.time()
  dt <- time.end - time.start
  run_time <- paste(round(dt, 2), attributes(dt)$units)
  
  norm.res <- list(dat.norm = dat.norm, 
                   norm.coef = coef.csn, 
                   RunTime = run_time)
  return(norm.res)
}


## L2-norm Normalization (L2N) ####
L2_Normalize <- function(dat) {
  N <- nrow(dat)
  
  time.start <- Sys.time()
  
  ss.all <- apply(dat, MARGIN = 1, FUN = function(x) sum(x^2))
  l2.all <- sqrt(ss.all)
  ss.m <- mean(ss.all)
  l2.m <- sqrt(ss.m)
  
  coef.l2n <- l2.m / l2.all
  dat.norm <- diag(coef.l2n) %*% dat
  
  # record running time
  time.end <- Sys.time()
  dt <- time.end - time.start
  run_time <- paste(round(dt, 2), attributes(dt)$units)
  
  norm.res <- list(dat.norm = dat.norm, 
                   norm.coef = coef.l2n, 
                   RunTime = run_time)
  
  return(norm.res)
}  


## Probabilistic Quotient Normalization (PQN) ####
PQN_normalize <- function(dat) {
  # dat: row is observation, column is feature
  N <- nrow(dat)
  p <- ncol(dat)
  
  time.start <- Sys.time()
  
  ref.samp <- rep(0, p)
  for (i in c(1: p)) {
    ref.samp[i] <- median(dat[, i], na.rm = T)
  }
  ratio_by_ref <- matrix(ref.samp, nrow = N, ncol = p, byrow = T) / dat
  
  coef.pqn <- rep(0, N)
  for (i in c(1: N)) {
    coef.pqn[i] <- median(ratio_by_ref[i, ], na.rm = T)
  }
  
  coef.pqn.mat <- coef.pqn %*% t(rep(1, p))
  dat.norm <- coef.pqn.mat * dat
  
  # record running time
  time.end <- Sys.time()
  dt <- time.end - time.start
  run_time <- paste(round(dt, 2), attributes(dt)$units)
  
  norm.res <- list(dat.norm = dat.norm, 
                   norm.coef = coef.pqn, 
                   RunTime = run_time)
  
  return(norm.res)
}


## Quantile normalization ####
Quantile_normalize <- function(dat) {
  
  time.start <- Sys.time()
  
  dat <- t(dat)
  dat.norm <- normalize.quantiles.robust(dat)
  dat.norm <- t(dat.norm)
  
  # record running time
  time.end <- Sys.time()
  dt <- time.end - time.start
  run_time <- paste(round(dt, 2), attributes(dt)$units)
  
  norm.res <- list(dat.norm = dat.norm, 
                   RunTime = run_time)
  
  return(norm.res)
}




