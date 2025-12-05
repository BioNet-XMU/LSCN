# Date: 20251102
# Author: Fanjing Guo
# Title: Functions of LSCN Algorithm

# PACKAGE ####
library(philentropy)
library(mixOmics)
library(igraph)
library(fpc)
library(rlang)

## Find the contaminated metabolites between samples ####
findContamVar <- function(dat, k = NA, print.txt = T) {
  
  N <- nrow(dat)
  p <- ncol(dat)
  
  met.name <- colnames(dat)
  
  # compute the original sample correlation
  samp.corr.ori <- cor(t(dat))
  
  if (is.na(k)) {
    ContMet.num <- round(0.05 * p)
    if (ContMet.num > 10) {
      ContMet.num <- 10 # a max number of 10
    }
  }
  else {
    ContMet.num <- k
  }
  
  if (ContMet.num == 0) {
    ContVar.mat <- NULL
  }
  else {
    out_metNum.mat <- matrix(0, N, N)
    ContVar.mat <- matrix('', N, N)
    all.out_met <- c()
    newSamp.corr <- list()
    c <- 0
    
    samppair.num <- N * (N - 1) / 2
    print.sn <- round(quantile(c(1: samppair.num), seq(0.01, 1, 0.01)))
    sn <- 0
    
    if (print.txt) {
      cat('Detecting contaminated metabolites: \n')
    }
    
    for (i in c(1: (N - 1))) {
      for (j in c((i + 1): N)) {
        sn <- sn + 1
        
        # find outlier metabolites between two samples
        modes <- log2(abs(dat[i, ] / dat[j, ]))
        modes.m <- mean(modes)
        modes.sd <- sd(modes)
        modes_min_out <- modes.m - 3 * modes.sd
        modes_max_out <- modes.m + 3 * modes.sd
        
        out_met <- which(modes < modes_min_out | modes > modes_max_out)
        out_met.Num <- length(out_met)
        out_metNum.mat[i, j] <- out_met.Num
        
        # get original rank for the correlation between sample i and j
        sampi.corr.ori <- samp.corr.ori[i, -i]
        sampij.rank.ori <- order(order(sampi.corr.ori, decreasing = T))[j - 1]
        
        sampj.corr.ori <- samp.corr.ori[j, -j]
        sampji.rank.ori <- order(order(sampj.corr.ori, decreasing = T))[i]
        
        # determine contaminated variables
        if (length(out_met) > ContMet.num) {
          all.drank <- rep(0, length(out_met))
          for (k in c(1: length(out_met))) {
            out_met.occur <- which(all.out_met == out_met[k])
            if (length(out_met.occur) == 0) {
              c <- c + 1
              all.out_met <- append(all.out_met, out_met[k])
              
              dat_delout <- dat[, -out_met[k]]
              
              samp.corr <- cor(t(dat_delout))
              newSamp.corr[[c]] <- samp.corr
            }
            else {
              samp.corr <- newSamp.corr[[out_met.occur]]
            }
            
            # get new rank for the correlation between sample i and j
            sampi.corr <- samp.corr[i, -i]
            sampij.rank <- order(order(sampi.corr, decreasing = T))[j - 1]
            
            sampj.corr <- samp.corr[j, -j]
            sampji.rank <- order(order(sampj.corr, decreasing = T))[i]
            
            # decide whether the metabolite is contaminated
            drank.ij <- abs(sampij.rank - sampij.rank.ori)
            drank.ji <- abs(sampji.rank - sampji.rank.ori)
            drank.max <- max(drank.ij, drank.ji)
            all.drank[k] <- drank.max
          }
          cont_met <- out_met[order(all.drank, decreasing = T)[c(1: ContMet.num)]]
        }
        else {
          modes_med <- median(modes)
          out_met.rankIdx <- order(abs(modes - modes_med), decreasing = T)
          cont_met <- out_met.rankIdx[c(1: ContMet.num)]
        }
        cont_met <- paste0(cont_met, collapse = ',')
        ContVar.mat[i, j] <- cont_met
        
        if (print.txt) {
          print.sn_10 <- print.sn[seq(10, 100, 10)]
          if (length(which(print.sn_10 == sn)) != 0) {
            cat(names(print.sn_10[print.sn_10 == sn]))
          }
          else if (length(which(print.sn == sn)) != 0) {
            cat('-')
          }
        }
      }
    }
    
    if (print.txt) {
      cat('  End\n')
    }
    
  }
  
  return(ContVar.mat)
}


## Calculate sample correlation with sample pair specific contaminated metabolites ####
getCorrwithCont <- function(dat, doLog = T, ContVar.mat = NULL) {
  N <- nrow(dat)
  
  if (doLog) {
    if (sum(dat <= 0) != 0) {
      stop('Data contains zero or negative values!')
    }
    else {
      dat <- log2(abs(dat))
    }
  }
  
  dat.sca <- scale(dat, center = T, scale = T)
  dat.pca <- pca(dat.sca, ncomp = NULL, center = F, scale = F)
  cum.prop.var <- dat.pca$cum.var
  pc.num <- min(which(cum.prop.var >= 0.95))
  if (pc.num < 3) {
    pc.num <- 3 # the minimum number of PCs
  }
  dat.pca.load <- dat.pca$loadings$X[, c(1: pc.num)]
  
  SampCorr <- diag(1, N, N)
  if (!is.null(ContVar.mat)) {
    for (i in c(1: (N - 1))) {
      for (j in c((i + 1): N)) {
        contmet_ij <- as.numeric(unlist(strsplit(ContVar.mat[i, j], ',')))
        load.noCont <- dat.pca.load
        if (length(contmet_ij) != 0) {
          # the loadings of contaminated metabolites are set to zero
          load.noCont[contmet_ij, ] <- 0 
        }
        
        # reconstruct the data using new loadings
        pca.score <- dat.sca %*% load.noCont
        pca.score.sca <- scale(pca.score, center = T, scale = T)
        pca.score.ij <- pca.score.sca[c(i, j), ]
        
        # compute the correlation between sample i and j 
        corr.ij <- cor(pca.score.ij[1, ], pca.score.ij[2, ])
        
        SampCorr[i, j] <- corr.ij
        SampCorr[j, i] <- SampCorr[i, j]
      }
    }
  }
  else {
    pca.score <- dat.pca$variates$X[, c(1: pc.num)]
    pca.score.sca <- scale(pca.score, center = T, scale = T)
    SampCorr <- cor(t(pca.score.sca))
  }
  
  
  return(SampCorr)
}


## Get similarity between samples from sample correlation ####
getSimifromcorr <- function(corr.mat, neighbour.thresh = 0.1, 
                             Simi.index = 'WRA') {
  
  N <- nrow(corr.mat)
  samp.similar <- matrix(0, N, N)
  
  corr.mat[corr.mat < neighbour.thresh] <- 0
  
  for (i in c(1: (N - 1))) {
    
    corr.i <- corr.mat[, i]
    for (j in c((i + 1): N)) {
      
      corr.j <- corr.mat[, j]
      
      i.neibIdx <- which(corr.i != 0)
      j.neibIdx <- which(corr.j != 0)
      
      coneibIdx <- intersect(i.neibIdx, j.neibIdx)
      
      # Common Neighbors (CN)
      if (Simi.index == 'CN') {
        if (length(coneibIdx) != 0) {
          CN_ij <- sum(corr.i[coneibIdx]) + sum(corr.j[coneibIdx])
        }
        else {
          CN_ij <- 0
        }
        samp.similar[i, j] <- CN_ij
        samp.similar[j, i] <- samp.similar[i, j]
      }
      
      # Jaccard’s Coefficient (JC)
      if (Simi.index == 'JC') {
        coneib.num <- length(coneibIdx)
        allneib.num <- length(unique(c(i.neibIdx, j.neibIdx)))
        
        JC_ij <- coneib.num / allneib.num
        
        samp.similar[i, j] <- JC_ij
        samp.similar[j, i] <- samp.similar[i, j]
        
      }
      
      # Salton Index (SI)
      if (Simi.index == 'SI') {
        coneib.num <- length(coneibIdx)
        i.degree <- length(i.neibIdx)
        j.degree <- length(j.neibIdx)
        
        SI_ij <- coneib.num / sqrt(i.degree * j.degree)
        
        samp.similar[i, j] <- SI_ij
        samp.similar[j, i] <- samp.similar[i, j]
        
      }
      
      # Weighted Preferential Attachment (WPA)
      if (Simi.index == 'WPA') {
        WPA_ij <- sum(corr.i[i.neibIdx]) * sum(corr.j[j.neibIdx])
        
        samp.similar[i, j] <- WPA_ij
        samp.similar[j, i] <- samp.similar[i, j]
        
      }
      
      # Weighted Adamic–Adar index (WAA)
      if (Simi.index == 'WAA') {
        if (length(coneibIdx) != 0) {
          WAA_ij.vec <- rep(0, length(coneibIdx))
          
          for (z in c(1: length(coneibIdx))) {
            wei.ijz.sum <- corr.i[coneibIdx[z]] * corr.j[coneibIdx[z]]
            corr.z <- corr.mat[, coneibIdx[z]]
            
            corr.z <- corr.z[-coneibIdx[z]]
            corr.z[corr.z < neighbour.thresh] <- 0
            
            wei.zneib.sum <- sum(corr.z[corr.z != 0])
            
            WAA_ij.vec[z] <- wei.ijz.sum / log(1 + wei.zneib.sum)
          }
          
          WAA_ij <- sum(WAA_ij.vec)
        }
        else {
          WAA_ij <- 0
        }
        
        samp.similar[i, j] <- WAA_ij
        samp.similar[j, i] <- samp.similar[i, j]
        
      }
      
      # Weighted Resource Allocation index(WRA)
      if (Simi.index == 'WRA') {
        
        if (length(coneibIdx) != 0) {
          wei.coneib.prods <- corr.i[coneibIdx] * corr.j[coneibIdx]
          corr.coneib <- corr.mat[, coneibIdx]
          if (length(coneibIdx) > 1) {
            wei.coneib_neib.sums <- colSums(corr.coneib)
          }
          else {
            wei.coneib_neib.sums <- sum(corr.coneib)
          }
          
          WRA_ij <- sum(wei.coneib.prods / wei.coneib_neib.sums)
        }
        else {
          WRA_ij <- 0
        }
        
        samp.similar[i, j] <- WRA_ij
        samp.similar[j, i] <- samp.similar[i, j]
        
      }
      
    }
  }
  
  # Standardized
  simi.min <- min(samp.similar[upper.tri(samp.similar)])
  simi.max <- max(samp.similar[upper.tri(samp.similar)])
  samp.similar <- (samp.similar - simi.min) / (simi.max - simi.min)
  diag(samp.similar) <- 0
  
  return(samp.similar)
}


## Determine a reasonable similarity threshold ####
getSimiThresh <- function(samp.similar, print.txt = F) {
  similars <- samp.similar[upper.tri(samp.similar)]
  
  simi.quant <- 0.7
  simi.thresh <- quantile(similars, simi.quant)
  
  # get adjacency matrix of the samples
  A.mat <- 1 * (samp.similar > simi.thresh)
  # construct a graph from the adjacency matrix
  G <- graph_from_adjacency_matrix(A.mat, weighted = NULL, mode = "undirected")
  
  # detect the conectivity of the graph
  connect.g <- is_connected(G)
  
  if (!connect.g) {
    simi.quants <- seq(0.7, 0.2, -0.05)
    
    for (i in c(1: length(simi.quants))) {
      simi.thresh <- quantile(similars, simi.quants[i])
      A.mat <- 1 * (samp.similar > simi.thresh)
      G <- graph_from_adjacency_matrix(A.mat, weighted = NULL, mode = "undirected")
      
      connect.g <- is_connected(G)
      
      if (connect.g) {
        if (i < length(simi.quants)) {
          set.quant <- simi.quants[i + 1]
          set.thresh <- quantile(similars, set.quant)
        }
        else {
          set.thresh <- simi.thresh
        }
        
        break
      }
      else {
        set.thresh <- simi.thresh
      }
    }
    
  }
  else {
    set.thresh <- simi.thresh
  }
  
  simi.perc <- sum(similars > set.thresh) / length(similars) * 100
  
  if(print.txt) {
    cat('Select a similarity threshold: ', 
        round(set.thresh, 4), paste0('(', round(simi.perc, 2), '%)'), '\n')
  }
  
  return(set.thresh)
}


## PQN used between two samples in the algorithm ####
TwoSamp_PQN <- function(S, S1) {

  if (sum(S <= 0) != 0) {
    S[S <= 0] <- 1e-5
  }
  
  if (sum(S1 <= 0) != 0) {
    S1[S1 <= 0] <- 1e-5
  }
  
  Q <- S1 / S
  norm.coef <- median(Q)
  
  return(norm.coef)
}


## Automatically adjust the learning rate ####
adjust_learning_rate <- function(conv_history, ita, adjust_freq = 5, 
                                 decayRate, incRate, 
                                 max_ita = 0.95, min_ita = 0.001) {
  if (length(conv_history) >= adjust_freq) {
    recent_conv <- tail(conv_history, adjust_freq)
    improvements <- diff(recent_conv)
    mean_improve <- mean(improvements)
    
    mean_conv <- mean(recent_conv)
    min_improvement <- 10^(floor(log10(mean_conv)) - 2)
    
    if (mean_improve > - min_improvement) {
      if (any(improvements > 0) & ita * decayRate > min_ita) {
        ita_new <- ita * decayRate
      }
      else if (all(improvements < 0) & ita * incRate < max_ita) {
        ita_new <- ita * incRate
      }
      else {
        ita_new <- ita
      }
    }
    else {
      ita_new <- ita
    }
  }
  else {
    ita_new <- ita
  }
  
  if (ita_new == ita) {
    adjusted <- F
  }
  else {
    adjusted  <- T
  }
  
  ita.res <- list(New_ita = ita_new, 
                  adjusted = adjusted)
  return(ita.res)
}


## LSCN (Main) ####
LSCN <- function(X, doLog = T, k = NA, corr.thresh = NA, similar.type = 'WRA', 
                 simi.quant = NA, ita0 = 0.4, adjust_freq = 5, 
                 decayRate = 0.8, incRate = 1.1, max_ite = 100, 
                 eps = 1e-4, print.txt = F) {
  
  # INPUT
  # X: a data matrix with samples in rows and features in columns.
  # doLog: (default=TRUE) logical, whether the data should be log-transformed.
  # k: (default=NA) integer, the number of contaminated features. If NA, k = 0.05*nrow(X).
  # corr.thresh: (default=NA) positive real, the sample PCC cut-off threshold. If NA, the adaptive selection strategy will be used.
  # similar.type: (default='WRA') one of ('CN','JC','SI','WPA','WAA','WRA'). Specifies the method used to estimate sample similarity.
  # simi.quant: (default=NA) a probability value which determines the sample similarity threshold. if NA, the adaptive selection strategy will be used.
  # ita0: (default=0.4) a numeric value in the range (0,0.5], the analogous learning rate.
  # adjust_freq: (default=5) integer, the number of iterations between each dynamic adjustment of the analogous learning rate.
  # decayRate: (default=0.8) a numeric value between 0 and 1, the decay rate employed in the dynamic adjustment process of the analogous learning rate.
  # incRate: (default=1.1) a numeric value larger than 1, the increase rate employed in the dynamic adjustment process of the analogous learning rate.
  # max_ite: (default=100) integer, the maximum number of iterations.
  # eps: (default=1e-4) positive real, the convergence tolerance.
  # print.txt: (default=FALSE) logical, whether to print progress messages.
  
  N0 <- nrow(X) # sample number
  p <- ncol(X) # feature number
  
  # Initialization
  X.norm <- X
  NORM.COEF <- rep(1, N0)
  Wei.NORM.COEF <- rep(1, N0)
  
  auto_thresh  <- F
  ite <- 0
  ita <- ita0
  conv.vals <- c()
  ita.rec <- c()
  
  ita_adjusted <- F
  ita_adjust_iteration <- 0
  
  time.start <- Sys.time()
  
  # Detect contaminated features
  Cont.mat <- findContamVar(X, k = k, print.txt = print.txt)
  while(TRUE) {
    
    # Initialize PQN normalization factors
    norm.coef <- rep(1, N0)
    
    # Calculate PCC of each sample pair without contaminated features
    samp.corr <- getCorrwithCont(X.norm, doLog = doLog, ContVar.mat = Cont.mat)
    
    # Determine PCC threshold
    if (is.na(corr.thresh)) {
      
      auto_thresh <- T
      corr.threshs <- rep(0, N0)
      for (i in c(1: N0)) {
        sc.i <- samp.corr[-i, i] # the PCCs with all other samples for sample i
        corr.threshs[i] <- quantile(sc.i, 0.7) # 0.7 quantile values of the PCCs for sample i
      }
      corr.thresh <- min(corr.threshs)
      corr.thresh <- max(corr.thresh, 0)
      
      if (print.txt) {
        cat('PCC threshold: ', corr.thresh, '\n')
      }
      
    }
    
    # get sample similarity from sample PCCs
    samp.similar <- getSimifromcorr(samp.corr, neighbour.thresh = corr.thresh,
                                    Simi.index = similar.type)
    
    # Determine similarity threshold
    if (is.na(simi.quant)) {
      simi.thresh <- getSimiThresh(samp.similar, print.txt = print.txt)
    }
    else {
      simi.thresh <- quantile(samp.similar[upper.tri(samp.similar)], simi.quant)
      samp.nw <- 1 * (samp.similar > simi.thresh)
      samp.neibNum <- colSums(samp.nw)
      
      while(any(samp.neibNum == 0) & simi.quant > 0.2) {
        simi.quant <- simi.quant - 0.05
        
        simi.thresh <- quantile(samp.similar[upper.tri(samp.similar)], simi.quant)
        samp.nw <- 1 * (samp.similar > simi.thresh)
        samp.neibNum <- colSums(samp.nw)
      }
      
      if(print.txt) {
        simi.perc <- sum(samp.nw[upper.tri(samp.nw)]) / (N0 * (N0-1) / 2) * 100
        cat('Select a similarity threshold: ', 
            round(simi.thresh, 4), paste0('(', round(simi.perc, 2), '%)'), '\n')
      }
    }
    
    # Normalization 
    for (i in c(1: N0)) {
      
      i.sampdat <- as.matrix(t(X.norm[i, ])) # data of sample i
      sampIdx <- c(1: N0)[-i]
      simi.sampi <- samp.similar[, i]
      simi.sampi <- simi.sampi[-i] # the sample similarity with all other samples for sample i
      
      # determine the neighbor set for sample i according to sample similarity
      near.sampIdx <- sampIdx[simi.sampi > simi.thresh]
      k <- length(near.sampIdx) # number of neighbors
      
      if (k != 0) {
        near.sampdat <- X.norm[near.sampIdx, ]
        near.simi <- samp.similar[near.sampIdx, i]
        
        # compute weight
        w <- near.simi / sum(near.simi)
        
        # Calculate PQN normalization factor for sample i to each neighbor
        normCoef <- rep(0, k)
        if (k == 1) {
          near.sampdat <- as.matrix(t(near.sampdat))
        }
        for (j in c(1: k)) {
          near.samp.j <- as.matrix(t(near.sampdat[j, ])) # data of neighbor j
          normCoef[j] <- TwoSamp_PQN(i.sampdat, near.samp.j)
        }
        
        # calculate weighted normalization factor (WNF)
        wei.normCoef <- sum(w * normCoef)
        
        # using learning rate
        norm.coef[i] <- 1 + ita * (wei.normCoef - 1)
        
        # record WNF
        Wei.NORM.COEF[i] <- wei.normCoef
      }
      
    }
    
    # Normalize all samples
    X.norm <- diag(norm.coef) %*% X.norm
    
    # Update normalization factor
    NORM.COEF <- NORM.COEF * norm.coef
    
    if (auto_thresh) {
      corr.thresh <- NA
    }

    conv <- median(abs(Wei.NORM.COEF - 1))
    conv.vals <- append(conv.vals, conv)
    
    ita.rec <- append(ita.rec, ita)
    
    ite <- ite + 1
    
    if (!ita_adjusted) {
      ita.res <- adjust_learning_rate(conv_history = conv.vals, adjust_freq = adjust_freq,
                                      ita = ita, decayRate = decayRate,
                                      incRate = incRate)
      ita_adjusted <- ita.res$adjusted
      ita <- ita.res$New_ita
    }
    
    if (ita_adjusted & ita_adjust_iteration < adjust_freq - 1) {
      ita_adjust_iteration <- ita_adjust_iteration + 1
    }
    else {
      ita_adjusted <- F
      ita_adjust_iteration <- 0
    }
    
    
    if (print.txt) {
      cat(paste0(round(conv, 5), ' (iteration: ', ite, ')'), '\n')
    }
    
    if (conv < eps) {
      break
    }
    else if (ite > 10 & ite < max_ite) {
      conv_last10 <- conv.vals[c((ite - 10): ite)]
      conv_deccor <- cor(c((ite - 10): ite), conv_last10, method = 'spearman')
      conv_last10.rsd <- sd(conv_last10) / mean(conv_last10)
      if (abs(conv_deccor) < 0.1 & conv_last10.rsd < 0.1) {
        break
      }
    }
    else if (ite == max_ite) {
      break
    }
  }
  
  # record running time
  time.end <- Sys.time()
  dt <- time.end - time.start
  run_time <- paste(round(dt, 2), attributes(dt)$units)
  
  Norm.Result <- list(dat.norm = X.norm,
                      NORM.COEF = NORM.COEF, 
                      conv.vals = conv.vals, 
                      ita_history = ita.rec, 
                      RunTime = run_time)
  
  # OUTPUT
  # dat.norm: the final normalized data matrix
  # NORM.COEF: the normalization factors applied to the original data.
  # conv.vals: the recorded convergence criterion value in each iteration.
  # ita_history: the used analogous learning rate in each iteration.
  # RunTime: run time.
  
  return(Norm.Result)
  
}
