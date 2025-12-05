# Date: 20251102
# Author: Fanjing Guo
# Title: Functions to Generate a Simulated Dataset

# PACKAGE ####
library(mixOmics)

## Generate Simulated Data ####
GenSimuDat <- function(X, N = 100, selsamp.ratio = 0.9, 
                       indiv.sd.ratio = 0.25, grpVar.ratio = 0.3, 
                       grpEff.lamda = c(1, 3), grpEff.sd.ratio = 0.2, 
                       subgrpNum = 3, topPC.num = 10, 
                       diltF.m = 10, diltF.sd = 3, print.txt = F) {
  
  # INPUT
  # X: sample-by-feature data matrix that serves as the source of the base spectra.
  # N: (default=100) integer, the sample size of the simulated dataset.
  # selsamp.ratio: (default=0.9) a probability value, the proportion of samples drawn from X in step 1.
  # indiv.sd.ratio: (default=0.25) a numeric value between 0 and 1, the ratio of the standard deviation of individual differences to the mean feature intensity.
  # grpVar.ratio: (default=0.3) a numeric value between 0 and 1, the proportion of randomly selected features with inter-group differences.
  # grpEff.lamda: (default=c(1,3)) a numeric vector with two elements, which determines the range of the inter-group difference magnitudes. 
  # grpEff.sd.ratio (default=0.3) a numeric value between 0 and 1, which also determines the range of the inter-group difference magnitudes.
  # subgrpNum: (default=3) integer, the number of sub-groups.
  # topPC.num: (default=10) integer, the number of principal components to add intra-group heterogeneity.
  # diltF.m: (default=10) a positive real value determines the distribution of dilution coefficients.
  # diltF.sd: (default=3) a positive real value determines the distribution of dilution coefficients.
  # print.txt: (default=FALSE) logical, whether to print progress messages.
  
  N0 <- nrow(X)
  p <- ncol(X)
  
  # Step 1: Generate base spectra
  if (print.txt) {
    cat('Generate base data...')
  }
  
  Base_dat <- matrix(0, N, p)
  for (i in c(1: N)) {
    selsamp <- X[sample(c(1: N0), round(selsamp.ratio * N0)), ]
    Base_dat[i, ] <- apply(selsamp, MARGIN = 2, FUN = median)
  }
  
  if (print.txt) {
    cat('finished.\n')
  }
  
  # Step 2: Add individual differences
  if (print.txt) {
    cat('Add individual differences...')
  }
  
  ref.peak <- apply(Base_dat, MARGIN = 2, FUN = mean)
  Ei.mat <- matrix(0, nrow = N, ncol = p)
  for (i in c(1: p)) {
    indiv_diff.sd <- ref.peak[i] * indiv.sd.ratio
    Ei.mat[, i] <- rnorm(N, mean = 0, sd = indiv_diff.sd)
    
    dat.i <- Base_dat[, i]
    dat.in <- dat.i + Ei.mat[, i]
    
    negIdx <- which(dat.in < 0)
    if (length(negIdx) != 0) {
      negVal.min <- min(dat.in[negIdx])
      dm.i <- 1e-5 - negVal.min # ensure the data values >= 1e-5
      Base_dat[, i] <- dat.i + dm.i
    }
  }
  
  Indiv_dat <- Base_dat + Ei.mat
  
  era.peak.num <- sum(as.vector((abs(Ei.mat) - Base_dat) >= 0))
  era.peak.ratio <- era.peak.num / (N * p)
  
  if (print.txt) {
    cat('finished.\n')
  }
  
  # Step 3: Add Inter-group Heterogeneity
  if (print.txt) {
    cat('Add inter-group heterogeneity...')
  }
  
  Group_dat <- Indiv_dat
  Indiv_dat.grp1 <- Indiv_dat[c(1: round(N/2)), ]
  Indiv_dat.grp2 <- Indiv_dat[-c(1: round(N/2)), ]
  
  grpVar.num <- round(grpVar.ratio * p)
  
  # set 70% variables up-regulated, other 30% down-regulated
  grpVar.up.num <- round(grpVar.num * 0.7)
  grpVar.down.num <- grpVar.num - grpVar.up.num
  
  sel.DEmet <- sample(c(1: p), grpVar.num)
  lamda.mat <- matrix(0, nrow(Indiv_dat.grp2), grpVar.num)
  
  # up-regulated
  if (grpVar.up.num > 0) {
    up.min <- - grpEff.lamda[2]
    up.max <- - grpEff.lamda[1]
    lamda.up <- runif(grpVar.up.num, min = up.min, max = up.max)
    m.g1 <- colMeans(Indiv_dat.grp1[, sel.DEmet[c(1: grpVar.up.num)]]) 
    m.g2 <- colMeans(Indiv_dat.grp2[, sel.DEmet[c(1: grpVar.up.num)]]) 
    lamda.up.m <- m.g1 / exp(lamda.up) - m.g2
    lamda.up.sd <- grpEff.sd.ratio * abs(lamda.up.m)
    lamda.up.mlog <- log(lamda.up.m) - log(sqrt((lamda.up.sd / lamda.up.m)^2 + 1))
    lamda.up.sdlog <- sqrt(log((lamda.up.sd / lamda.up.m)^2 + 1))
    for (i in c(1: grpVar.up.num)) {
      lamda.mat[, i] <- rlnorm(nrow(Indiv_dat.grp2), lamda.up.mlog[i], lamda.up.sdlog[i])
    }
  }
  
  # down-regulated
  if (grpVar.down.num > 0) {
    down.min <- grpEff.lamda[1]
    down.max <- grpEff.lamda[2]
    lamda.down <- runif(grpVar.down.num, min = down.min, max = down.max)
    if (grpVar.down.num > 1) {
      m.g1 <- colMeans(Indiv_dat.grp1[, sel.DEmet[c((grpVar.up.num + 1): grpVar.num)]]) 
      m.g2 <- colMeans(Indiv_dat.grp2[, sel.DEmet[c((grpVar.up.num + 1): grpVar.num)]]) 
    }
    else {
      m.g1 <- mean(Indiv_dat.grp1[, sel.DEmet[grpVar.num]]) 
      m.g2 <- mean(Indiv_dat.grp2[, sel.DEmet[grpVar.num]]) 
    }
    
    lamda.down.m <- m.g2 - m.g1 / exp(lamda.down)
    lamda.down.sd <- grpEff.sd.ratio * abs(lamda.down.m)
    lamda.down.mlog <- log(lamda.down.m) - log(sqrt((lamda.down.sd / lamda.down.m)^2 + 1))
    lamda.down.sdlog <- sqrt(log((lamda.down.sd / lamda.down.m)^2 + 1))
    for (i in c(1: grpVar.down.num)) {
      lamda.mat[, grpVar.up.num + i] <- - rlnorm(nrow(Indiv_dat.grp2), lamda.down.mlog[i], lamda.down.sdlog)
    }
  }
  
  Indiv_dat.grp2[, sel.DEmet] <- Indiv_dat.grp2[, sel.DEmet] + lamda.mat
  Indiv_dat.grp2 <- abs(Indiv_dat.grp2)
  Group_dat[c(1: round(N/2)), ] <- Indiv_dat.grp1
  Group_dat[-c(1: round(N/2)), ] <- Indiv_dat.grp2
  
  # group
  group <- rep(0, N)
  group[c(1: round(N/2))] <- 1
  group[-c(1: round(N/2))] <- 2
  
  if (print.txt) {
    cat('finished.\n')
  }
  
  # Step 4: Add Intra-group Heterogeneity
  if (print.txt) {
    cat('Add intra-group heterogeneity...')
  }
  
  Heter_dat <- Group_dat
  
  if (subgrpNum > 1) {
    
    candi.Eff_down <- seq(-0.4, -0.1, 0.01)
    candi.Eff_up <- seq(0.1, 0.4, 0.01)
    candi.Eff <- c(candi.Eff_down, candi.Eff_up)
    
    heterEff <- matrix(0, subgrpNum - 1, topPC.num)
    for (i in c(1: topPC.num)) {
      heterEff[, i] <- exp(sample(candi.Eff, subgrpNum - 1))
    }
    
    # For group 1
    for (g in c(1: 2)) {
      dat.grpg <- Group_dat[group == g, ]
      dat.grpg.sca <- scale(dat.grpg, center = T, scale = T)
      sca.attr <- attributes(dat.grpg.sca)
      dat.grpg.m <- sca.attr$`scaled:center`
      dat.grpg.sd <- sca.attr$`scaled:scale`
      dat.grpg.pca <- pca(dat.grpg.sca, ncomp = NULL, 
                          center = F, scale = F)
      dat.grpg.load <- dat.grpg.pca$loadings$X
      dat.grpg.pcs <- dat.grpg.pca$variates$X
      
      grpg.num <- sum(group == g)
      subgrp.sampnum <- round(grpg.num / subgrpNum)
      
      for (i in c(2: subgrpNum)) {
        heterEff.i <- heterEff[i - 1, ]
        
        start.index <- (i - 1) * subgrp.sampnum + 1
        if (i != subgrpNum) {
          end.index <- i * subgrp.sampnum
        }
        else {
          end.index <- grpg.num
        }
        
        dat.grpg.pcs[c(start.index: end.index), c(1: topPC.num)] <- 
          dat.grpg.pcs[c(start.index: end.index), c(1: topPC.num)] %*% diag(heterEff.i)
      }
      
      dat.grpg.sca <- dat.grpg.pcs %*% t(dat.grpg.load)
      dat.grpg <- dat.grpg.sca %*% diag(dat.grpg.sd)
      dat.grpg <- dat.grpg + rep(1, nrow(dat.grpg)) %*% t(dat.grpg.m)
      
      Heter_dat[group == g, ] <- dat.grpg
    }
    
  }
  
  Heter_dat <- abs(Heter_dat)
  simu_init <- Heter_dat
  
  if (print.txt) {
    cat('finished.\n')
  }
  
  # Step 5: Add Dilution Effects
  if (print.txt) {
    cat('Add dilution effects...')
  }
  
  diltF.meanlog <- log(diltF.m) - log(sqrt((diltF.sd / diltF.m)^2 + 1))
  diltF.sdlog <- sqrt(log((diltF.sd / diltF.m)^2 + 1))
  
  diltF.alpha <- rlnorm(N, meanlog = diltF.meanlog, sdlog = diltF.sdlog)
  diltF.A <- diag(diltF.alpha)
  Simu_dat <- diltF.A %*% Heter_dat
  
  if (print.txt) {
    cat('finished.\n')
  }
  
  # Output results
  para <- list(Ei = Ei.mat, 
               era.peak.ratio = era.peak.ratio,
               group = group,
               DEmetIdx = sel.DEmet,
               grpEff = lamda.mat,
               diltF = diltF.alpha)
  
  Simu_result <- list(BaseDat = Base_dat, 
                      IndivDat = Indiv_dat, 
                      GrpDat = Group_dat, 
                      HeterDat = Heter_dat, 
                      SimuDat = Simu_dat, 
                      InitDat = simu_init, 
                      Para = para)
  
  # OUTPUT
  # BaseDat: the matrix of simulated base data.
  # IndivDat: the simulated data matrix after adding individual differences.
  # GrpDat: the simulated data matrix after adding individual differences and inter-group differences.
  # HeterDat: the simulated data matrix after adding individual differences, inter-group differences and intra-group heterogeneity.
  # SimuDat: the simulated data matrix after adding individual differences, inter-group differences, intra-group heterogeneity and dilution effects.
  # InitDat: the simulated data matrix before adding dilution effects.
  # Para: a list containing six elements: 
  # Ei: the matrix of the added individual differences;
  # era.peak.ratio: the ratio of peaks that are masked by individual differences;
  # group: the group label;
  # DEmetIdx: the indices of inter-group differential features;
  # grpEff: the matrix of the added inter-group differences;
  # diltF: the dilution coefficients.
  
  return(Simu_result)
}