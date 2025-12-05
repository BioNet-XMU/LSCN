# Date: 20250618
# Author: Fanjing Guo
# Title: Functions used in Evaluation

# PACKAGE ####
library(ropls)
library(networkD3)
library(Rgraphviz)
library(igraph)
library(RBGL)
library(MASS)
library(car)
library(factoextra)
library(cluster)
library(e1071)
library(pROC)
library(ggsci)
library(colorspace)
library(ggstar)
library(ggraph)
library(ggplot2)
library(ggExtra)
library(ggforce)
library(caret)
library(MLmetrics)
library(progress)
library(VennDiagram)
library(nlme)
library(DMwR2)
library(energy)
library(ROCR)


## Pareto Scaling ####
Pareto_scale <- function(x) {
  N <- nrow(x)
  p <- ncol(x)
  x <- scale(x, center = T, scale = F)
  x.sd <- apply(x, MARGIN = 2, FUN = sd)
  x.sd.mat <- matrix(x.sd, nrow = N, ncol = p, byrow = T)
  x.par <- x / sqrt(x.sd.mat)
  
  return(x.par)
}


## Get PCA scores ####
getPCscore <- function(dat, scale.type = 'auto', delmet = NULL, 
                       pc.num = NA, VarPerc = 0.95, result = 'score') {
  
  if (!is.null(delmet)) {
    dat <- dat[, -delmet]
  }
  else {
    dat.sd <- apply(dat, MARGIN = 2, FUN = sd)
    sd.zeroIdx <- which(dat.sd == 0)
    if (length(sd.zeroIdx) != 0) {
      dat <- dat[, -sd.zeroIdx]
    }
  }
  
  if (scale.type == 'pareto') {
    dat.sca <- Pareto_scale(dat)
  }
  
  if(scale.type == 'auto') {
    dat.sca <- scale(dat, center = T, scale = T)
  }
  
  if(scale.type == 'min-max') {
    dat.min <- apply(dat, MARGIN = 2, FUN = min)
    dat.max <- apply(dat, MARGIN = 2, FUN = max)
    
    dat.sca <- dat
    for (i in c(1: ncol(dat))) {
      dat.sca[, i] <- (dat[, i] - dat.min[i]) / (dat.max[i] - dat.min[i])
    }
  }
  
  if (scale.type == 'none') {
    dat.sca <- dat
  }
  
  dat.pca <- pca(dat.sca, ncomp = NULL, center = F, scale = F)
  
  if (is.null(pc.num)) {
    pc.num <- ncol(dat.pca$variates$X)
  }
  else if (is.na(pc.num)) {
    cumVar <- dat.pca$cum.var
    pc.num <- min(which(cumVar >= VarPerc))
  }
  
  pca.varX <- dat.pca$variates$X[, c(1: pc.num)]
  
  if (result == 'score') {
    pca.res <- pca.varX
  }
  else if (result == 'pca') {
    pca.load <- dat.pca$loadings$X[, c(1: pc.num)]
    expl.var <- dat.pca$prop_expl_var$X
    cum.var <- dat.pca$cum.var
    
    pca.res <- list(score = pca.varX, 
                    loading = pca.load, 
                    explVar = expl.var, 
                    cumVar = cum.var)
  }
  
  return(pca.res)
}


## Compute the DHS Index ####
NetworkHeter <- function(A) {
  N <- nrow(A)
  d <- colSums(A)
  
  meanD <- sum(d) / N
  meanD2 <- sum(d^2) / N
  
  nw.heter <- meanD2 / (meanD^2)
  
  return(nw.heter)
}

Heter.Similar <- function(dat1, dat2, group = NULL, doLog = F) {
  
  N <- nrow(dat1)
  
  # dat1 is the reference or ground truth
  if (is.null(group)) {
    group <- rep(1, N)
  }
  
  grp.f <- as.factor(group)
  grp.lev <- levels(grp.f)
  grp.num <- nlevels(grp.f)
  
  if (doLog) {
    dat1 <- log(dat1)
    dat2 <- log(dat2)
  }
  
  rem.quant <- seq(0.1, 0.9, 0.1)
  heter.simi <- matrix(0, length(rem.quant), grp.num)
  for (i in c(1: grp.num)) {
    dat1.grpi <- dat1[group == grp.lev[i], ]
    dat2.grpi <- dat2[group == grp.lev[i], ]
    N.grpi <- sum(group == grp.lev[i])
    
    dat1.grpi.sd <- apply(dat1.grpi, MARGIN = 2, FUN = sd)
    del1.sd0 <- which(dat1.grpi.sd == 0)
    
    dat2.grpi.sd <- apply(dat2.grpi, MARGIN = 2, FUN = sd)
    del2.sd0 <- which(dat2.grpi.sd == 0)
    
    del.sd0 <- unique(c(del1.sd0, del2.sd0))
    if (length(del.sd0) != 0) {
      dat1.grpi <- dat1.grpi[, -del.sd0]
      dat2.grpi <- dat2.grpi[, -del.sd0]
    }
    
    score1.grpi <- getPCscore(dat1.grpi)
    score2.grpi <- getPCscore(dat2.grpi)
    
    eu.dist1 <- as.matrix(dist(score1.grpi))
    eu.dist2 <- as.matrix(dist(score2.grpi))
    
    hsimi <- rep(0, length(rem.quant))
    dCorr <- rep(0, length(rem.quant))
    for (j in c(1: length(rem.quant))) {
      all.dist1 <- eu.dist1[upper.tri(eu.dist1)]
      dist.thresh1 <- quantile(all.dist1, rem.quant[j])
      A1 <- 1 * (eu.dist1 <= dist.thresh1)
      diag(A1) <- 0
      
      all.dist2 <- eu.dist2[upper.tri(eu.dist2)]
      dist.thresh2 <- quantile(all.dist2, rem.quant[j])
      A2 <- 1 * (eu.dist2 <= dist.thresh2)
      diag(A2) <- 0
      
      # Correlation of the degrees 
      d1 <- colSums(A1)
      d2 <- colSums(A2)
      if (sd(d1) != 0 & sd(d2) != 0) {
        dCorr <- cor(d1, d2, method = 'pearson')
        dCorr <- (dCorr + 1) / 2
      }
      else {
        dCorr <- NA
      }
      
      # Similarity of degree heterogeneity
      h1 <- NetworkHeter(A1)
      h2 <- NetworkHeter(A2)
      hsimi <- exp(- abs(h2 - h1) / h1)
      
      # Similarity of Overall Heterogeneity
      heter.simi[j, i] <- 2 * dCorr * hsimi / (dCorr + hsimi)
    }
  }

  heter.simi <- median(rowMeans(heter.simi, na.rm = T))
  
  return(heter.simi)
}


## Compute the LSS Index ####
struct.similar <- function(dat1, dat2, group = NULL, doLog = F) {
  
  N <- nrow(dat1)
  
  if (is.null(group)) {
    group <- rep(1, N)
  }
  
  grp.f <- as.factor(group)
  grp.lev <- levels(grp.f)
  grp.num <- nlevels(grp.f)
  
  if (doLog) {
    dat1 <- log(dat1)
    dat2 <- log(dat2)
  }
  
  kr <- seq(0.1, 0.3, 0.1)
  Ms <- matrix(0, length(kr), grp.num)
  for (i in c(1: grp.num)) {
    dat1.grpi <- dat1[group == grp.lev[i], ]
    dat2.grpi <- dat2[group == grp.lev[i], ]
    N.grpi <- sum(group == grp.lev[i])
    
    dat1.grpi.sd <- apply(dat1.grpi, MARGIN = 2, FUN = sd)
    del1.sd0 <- which(dat1.grpi.sd == 0)
    
    dat2.grpi.sd <- apply(dat2.grpi, MARGIN = 2, FUN = sd)
    del2.sd0 <- which(dat2.grpi.sd == 0)
    
    del.sd0 <- unique(c(del1.sd0, del2.sd0))
    if (length(del.sd0) == 0) {
      del.sd0 <- NULL
    }
    
    score1.grpi <- getPCscore(dat1.grpi, delmet = del.sd0)
    score2.grpi <- getPCscore(dat2.grpi, delmet = del.sd0)
    
    eu.dist1 <- as.matrix(dist(score1.grpi))
    eu.dist2 <- as.matrix(dist(score2.grpi))
    
    for (ii in c(1: length(kr))) {
      k <- ceiling(kr[ii] * N.grpi)
      
      sum.rkU <- rep(0, N.grpi)
      sum.rkV <- rep(0, N.grpi)
      for (j in c(1: N.grpi)) {
        sampj1.dist <- eu.dist1[, j]
        sampj2.dist <- eu.dist2[, j]
        
        j_ <- c(1: N.grpi)[-j]
        rankIdx1 <- j_[order(sampj1.dist[j_])]
        rankIdx2 <- j_[order(sampj2.dist[j_])]
        
        rankIdx1.k <- rankIdx1[c(1: k)]
        rankIdx2.k <- rankIdx2[c(1: k)]
        
        # sample in knn of dat but not in knn of dat.ori
        Uk_j <- setdiff(rankIdx2.k, rankIdx1.k)
        if (length(Uk_j) != 0) {
          drk <- rep(0, length(Uk_j))
          for (u in c(1: length(Uk_j))) {
            drk[u] <- which(rankIdx1 == Uk_j[u]) - k
          }
          sum.rkU[j] <- sum(drk)
        }
        else {
          sum.rkU[j] <- 0
        }
        
        # sample in knn of dat.ori but not in knn of dat
        Vk_j <- setdiff(rankIdx1.k, rankIdx2.k)
        if (length(Vk_j) != 0) {
          drk <- rep(0, length(Vk_j))
          for (v in c(1: length(Vk_j))) {
            drk[v] <- which(rankIdx2 == Vk_j[v]) - k
          }
          sum.rkV[j] <- sum(drk)
        }
        else {
          sum.rkV[j] <- 0
        }
      }
      
      sum.rkU <- sum(sum.rkU)
      sum.rkV <- sum(sum.rkV)
      
      C <- 2 / (N.grpi * k * (2 * N.grpi - 3 * k - 1))
      M1.k <- 1 - C * sum.rkU
      M2.k <- 1 - C * sum.rkV
      
      M.k <- 2 * M1.k * M2.k / (M1.k + M2.k)
      
      Ms[ii, i] <- M.k
    }
  }
  
  M <- median(rowMeans(Ms))
  
  return(M)
}


## Improved k-means Clustering ####
# construct a coclass probability network 
CocProbNet <- function(dat, k = 5, repT = 1000, thresh = NULL, 
                       pc.num = NA, dolog = F, reduc_dim = T, 
                       weight = F) {
  
  N <- nrow(dat)
  if (dolog) {
    dat <- log2(dat)
  }
  
  dat.sd <- apply(dat, MARGIN = 2, FUN = sd)
  del.sd0 <- which(dat.sd == 0)
  if (length(del.sd0) != 0) {
    dat <- dat[, -del.sd0]
  }
  dat <- scale(dat, center = T, scale = T)
  
  if (reduc_dim) {
    dat <- getPCscore(dat, scale.type = 'none', pc.num = pc.num)
  }
  
  cocfreq.mat <- matrix(0, N, N)
  for (ii in c(1: length(k))) {
    ki <- k[ii]
    for (t in c(1: repT)) {
      km <- kmeans(dat, centers = ki, nstart = 10)
      dat.clust <- km$cluster
      
      for (i in c(1: (N - 1))) {
        i.clustId <- dat.clust[i]
        js.Idx <- c((i + 1): N)
        js.clustId <- dat.clust[js.Idx]
        
        coc.sampIdx <- js.Idx[js.clustId == i.clustId]
        cocfreq.mat[i, coc.sampIdx] <- cocfreq.mat[i, coc.sampIdx] + 1
        cocfreq.mat[coc.sampIdx, i] <- cocfreq.mat[coc.sampIdx, i] + 1
      }
    }
  }
  
  cocProb <- cocfreq.mat / (length(k) * repT)
  
  if (!is.null(thresh)) {
    cocNet <- 1 * (cocProb >= thresh)
    
    if (weight) {
      samp.simi <- as.matrix(distance(dat, method = 'cosine', 
                                      mute.message = T))
      samp.simi <- (samp.simi + 1) / 2
      cocNet <- cocNet * samp.simi
    }
  }
  else {
    if (weight) {
      samp.simi <- as.matrix(distance(dat, method = 'cosine', 
                                      mute.message = T))
      samp.simi <- (samp.simi + 1) / 2
      cocNet <- cocProb * samp.simi
    }
    else {
      cocNet <- cocProb
    }
  }
  
  return(cocNet)
  
}

ImproKmCluster <- function(dat, k = 3, thresh = NULL, doLog = F, 
                           reduce_dimension = T, weighted = F, 
                           optimal = T, varyk = 1, sig = c(0.5, 0.5), repT = 1000) {
  
  if (!is.na(varyk)) {
    if (k - varyk < 2) {
      ks <- c(2 : (k + varyk))
    }
    else {
      ks <- c((k - varyk) : (k + varyk))
    }
  }
  else {
    ks <- k
  }
  
  if (optimal) {
    thresh <- NULL
    A <- CocProbNet(dat, k = ks, thresh = thresh, weight = weighted, 
                    dolog = doLog, reduc_dim = reduce_dimension, repT = repT)
    max.prob <- max(A)
    
    threshs <- seq(0.5, max.prob, 0.1)
    modulas <- rep(0, length(threshs))
    clustNums <- rep(0, length(threshs))
    cmty.list <- list()
    for (i in c(1: length(threshs))) {
      thresh <- threshs[i]
      A1 <- 1 * (A >= thresh)
      g <- graph_from_adjacency_matrix(A1, mode = 'undirected', 
                                       weighted = TRUE, 
                                       diag = FALSE)
      cmty.i <- cluster_louvain(g, weights = NULL)
      modulas[i] <- modularity(cmty.i)
      clustNums[i] <- max(membership(cmty.i))
      cmty.list[[i]] <- cmty.i
    }
    
    if (any(duplicated(modulas))) {
      modulas.uni <- unique(modulas)
      modulas.uni.rank <- rank(-modulas.uni)
      modulas.rank <- rep(0, length(modulas))
      for (j in c(1: length(modulas.uni))) {
        modulas.rank[modulas == modulas.uni[j]] <- modulas.uni.rank[j]
      }
      modulas.rank <- exp(modulas.rank)
    }
    else {
      modulas.rank <- exp(rank(-modulas))
    }
    
    clustNums.diff <- abs(clustNums - mean(ks))
    if (any(duplicated(clustNums.diff))) {
      clustNums.diff.uni <- unique(clustNums.diff)
      clustNums.diff.uni.rank <- rank(clustNums.diff.uni)
      clustNum.diff.rank <- rep(0, length(clustNums.diff))
      for (j in c(1: length(clustNums.diff.uni))) {
        clustNum.diff.rank[clustNums.diff == clustNums.diff.uni[j]] <- clustNums.diff.uni.rank[j]
      }
      clustNum.diff.rank <- exp(clustNum.diff.rank)
    }
    else {
      clustNum.diff.rank <- exp(rank(clustNums.diff))
    }
    
    clustPerf.rank <- sig[1] * modulas.rank + sig[2] * clustNum.diff.rank
    cmty.idx <- which(clustPerf.rank == min(clustPerf.rank))
    
    cmty.idx <- cmty.idx[length(cmty.idx)]
    cmty <- cmty.list[[cmty.idx]]
  }
  else {
    A <- CocProbNet(dat, k = k, thresh = thresh, weight = weighted, 
                    reduc_dim = reduce_dimension, repT = repT)
    g <- graph_from_adjacency_matrix(A, mode = 'undirected', 
                                     weighted = TRUE, 
                                     diag = FALSE)
    cmty <- cluster_louvain(g, weights = NULL)
  }
  
  clusId <- membership(cmty)
  
  return(clusId)
}


## Metrics about clustering ####
getVmeasure <- function(cluster, group, beta = 1) {
  N <- length(cluster)
  
  clust.f <- as.factor(cluster)
  cnum <- nlevels(clust.f)
  clust.lev <- levels(clust.f)
  
  grp.f <- as.factor(group)
  gnum <- nlevels(grp.f)
  grp.lev <- levels(grp.f)
  
  Nc_vec <- rep(0, gnum)
  Hc_vec <- rep(0, gnum)
  for (c in c(1: gnum)) {
    Nc_vec[c] <- sum(group == grp.lev[c])
    tmp <- Nc_vec[c] / N
    Hc_vec[c] <- - tmp * log(tmp)
  }
  HC <- sum(Hc_vec)
  
  Nk_vec <- rep(0, cnum)
  Hk_vec <- rep(0, cnum)
  for (k in c(1: cnum)) {
    Nk_vec[k] <- sum(cluster == clust.lev[k])
    tmp <- Nk_vec[k] / N
    Hk_vec[k] <- - tmp * log(tmp)
  }
  HK <- sum(Hk_vec)
  
  Hck_vec <- c()
  Hkc_vec <- c()
  for (c in c(1: gnum)) {
    Nc <- Nc_vec[c]
    for (k in c(1: cnum)) {
      Nk <- Nk_vec[k]

      Nck <- sum(group == grp.lev[c] & cluster == clust.lev[k])
      
      if (Nck == 0) {
        Hck <- 0
        Hkc <- 0
      }
      else {
        tmp1 <- Nck / N
        tmp2 <- Nck / Nk
        tmp3 <- Nck / Nc
        
        Hck <- - tmp1 * log(tmp2)
        Hkc <- - tmp1 * log(tmp3)
      }
      
      Hck_vec <- append(Hck_vec, Hck)
      Hkc_vec <- append(Hkc_vec, Hkc)
    }
  }
  
  H_CK <- sum(Hck_vec)
  H_KC <- sum(Hkc_vec)
  
  # Homogeneity
  if (H_CK == 0) {
    h <- 1
  }
  else {
    h <- 1 - H_CK / HC
  }
  
  # Completeness
  if (H_KC == 0) {
    c <- 1
  }
  else {
    c <- 1 - H_KC / HK
  }
  
  # V-measure
  if (h == 0 & c == 0) {
    v <- 0
  }
  else {
    v <- (1 + beta) * h * c / (beta * h + c)
  }
  
  Vmeasure <- list(Homogeneity = h, 
                   Completeness = c, 
                   V_measure = v)
  
  return(Vmeasure)
  
}


## PCA Visualization ####
PCA_scatter_withQC <- function(data, group, xlim, ylim, 
                               x_label = NULL, y_label = NULL, 
                               legend.pos = 'inside', 
                               legend.inspos = c(0.85, 0.85), 
                               title = NULL, title.cex = 3, fon = 'sans') {
  
  dat <- as.data.frame(data)
  
  grp.f <- as.factor(group[group != 'QC'])
  s.lev <- levels(grp.f)
  group <- factor(group, 
                  levels = c('QC', s.lev))
  
  set_col <- c("black", "blue", "red", "green", 
               "magenta", "yellow", "cyan", "gray")
  
  # Plot
  set_shape <- c(16, 17, 18, 10, 9, 8, 13, 11)
  
  x.len <- xlim[2] - xlim[1]
  y.len <- ylim[2] - ylim[1]
  
  dx <- x.len / 8
  dy <- y.len / 18
  
  pMain <- ggplot(data = dat, aes(x = dat[, 1], y = dat[, 2])) + 
    geom_point(aes(color = group, shape = group), size = 5) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_hline(yintercept = 0, linetype = 2) + 
    annotate(geom = 'text', x = xlim[1] + dx, y = ylim[2] - dy, 
             label = title, size = 8, color = 'black', family = 'sans', fontface = 'bold') + 
    scale_color_manual(values = set_col) + 
    scale_shape_manual(values = set_shape) + 
    theme_bw() + 
    theme(legend.position = legend.pos, 
          legend.position.inside = legend.inspos, 
          axis.title.x = element_text(size = 25, family = fon, face = 'bold'), 
          axis.title.y = element_text(size = 25, family = fon, face = 'bold'), 
          axis.text.x = element_text(size = 22, family = fon, face = "bold", color = 'black'), 
          axis.text.y = element_text(size = 22, family = fon, face = "bold", color = 'black'), 
          axis.ticks = element_line(color = 'black', linewidth = 1), 
          axis.ticks.length = unit(2, 'mm'), 
          panel.border = element_rect(linewidth = 2), 
          panel.grid = element_blank(), 
          plot.title = element_text(size = 25, family = fon, face = 'bold', hjust = 0.5)) + 
    xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2]) + 
    xlab(x_label) + ylab(y_label)
  
  
  if (legend.pos != 'none') {
    pMain <- pMain + theme(legend.title = element_text(size = 18, family = fon, face = 'bold'), 
                           legend.text = element_text(size = 18, family = fon, face = 'bold'), 
                           legend.background = element_rect(color = 'black', linewidth = 0.8)) + 
      labs(color = 'group', shape = 'group', title = NULL)
  }
  else {
    pMain <- pMain + labs(title = NULL)
  }
  
  pMain
  
}


## get DE metabolites (Multiple groups) ####
# get VIP values for metabolites
getVIPforMet <- function(dat, group, CvI = 5, doLog = F) {
  
  if(doLog) {
    dat <- log2(1 + dat)
  }
  
  dat.plsda <- opls(dat, group, orthoI = 0, predI = 3, crossvalI = CvI, 
                    fig.pdfC = 'none', info.txtC = 'none')
  met.VIP <- getVipVn(dat.plsda)
  
  return(met.VIP)
}

getDEmet <- function(dat, group, doLog = T, p.level = 0.05, 
                      VIP.level = 1, adjustP = F) {
  p <- ncol(dat)
  
  grp.f <- as.factor(group)
  grp.lev <- levels(grp.f)
  grp.num <- nlevels(grp.f)
  
  if (doLog) {
    dat1 <- log2(dat)
  }
  else {
    dat1 <- dat
  }
  
  # anova, p-value
  aov.Pval <- rep(0, p)
  aov.Fval <- rep(0, p)
  for (i in c(1: p)) {
    var_zero <- F
    for (g in c(1: grp.num)) {
      dat.i.grpg <- dat1[group == grp.lev[g], i]
      sd.i.grpg <- sd(dat.i.grpg)
      
      if (sd.i.grpg == 0) {
        var_zero <- T
        aov.Pval[i] <- NA
        aov.Fval[i] <- NA
        break
      }
    }
    
    if (!var_zero) {
      meti.dat <- data.frame(grp = grp.f, 
                             val = dat1[, i])
      
      levT <- leveneTest(val ~ grp, data = meti.dat)
      levTPval <- levT$`Pr(>F)`[1]
      
      if (levTPval < 0.05){
        varEqual <- F
      } else{
        varEqual <- T
      }
      
      aov.i <- oneway.test(val ~ grp, data = meti.dat, 
                           var.equal = varEqual)
      aov.Pval[i] <- aov.i$p.value
      aov.Fval[i] <- aov.i$statistic
    }
  }
  
  if (adjustP) {
    aov.Pval <- p.adjust(aov.Pval, method = 'fdr')
  }
  
  aov.DEmet <- which(aov.Pval < p.level)
  
  # plsda, VIP
  dat.sd <- apply(dat, MARGIN = 2, FUN = sd)
  del.sd0 <- which(dat.sd == 0)
  
  if (length(del.sd0) != 0) {
    met.VIP <- rep(NA, p)
    met.VIP[-del.sd0] <- getVIPforMet(dat[, -del.sd0], 
                                      grp.f)
  } 
  else {
    met.VIP <- getVIPforMet(dat, grp.f)
  }
  
  
  VIP.DEmet <- which(met.VIP > VIP.level)
  
  DEmetIdx <- intersect(aov.DEmet, VIP.DEmet)
  
  DEmet.res <- list(DEmetIdx = DEmetIdx, 
                    Pval = aov.Pval, 
                    Fval = aov.Fval, 
                    VIP = met.VIP)
  return(DEmet.res)
}

getdiffPCC <- function(dat1, dat2, cor.method = 'pearson', 
                       conf_alpha = 0.05, topk = NA) {
  
  N <- nrow(dat1)
  p <- ncol(dat1)
  
  # compute PCC
  sd1 <- apply(dat1, MARGIN = 2, FUN = sd)
  sd2 <- apply(dat2, MARGIN = 2, FUN = sd)
  sd0_idx1 <- which(sd1 == 0)
  sd0_idx2 <- which(sd2 == 0)
  sd0_idx <- union(sd0_idx1, sd0_idx2)
  if (length(sd0_idx) != 0) {
    pcc1 <- diag(1, p, p)
    pcc1[-sd0_idx, -sd0_idx] <- cor(dat1[, -sd0_idx], method = cor.method)
    
    pcc2 <- diag(1, p, p)
    pcc2[-sd0_idx, -sd0_idx] <- cor(dat2[, -sd0_idx], method = cor.method)
  }
  else {
    pcc1 <- cor(dat1, method = cor.method)
    pcc2 <- cor(dat2, method = cor.method)
  }
  
  pcc1 <- round(pcc1, 5)
  pcc2 <- round(pcc2, 5)
  diag(pcc1) <- 0
  diag(pcc2) <- 0
  
  # process correlation = 1
  if (sum(pcc1 == 1) != 0 | sum(pcc2 == 1) != 0) {
    pcc1[pcc1 == 1] <- 0.99999
    pcc2[pcc2 == 1] <- 0.99999
  }
  
  # Fisher-z transformation
  fz1 <- 1 / 2 * log((matrix(1, p, p) + pcc1) / (matrix(1, p, p) - pcc1))
  fz2 <- 1 / 2 * log((matrix(1, p, p) + pcc2) / (matrix(1, p, p) - pcc2))
  
  # Z test
  SE <- sqrt(2 / (N - 3))
  dfz <- abs(fz1 - fz2)
  Z <- dfz / SE
  # Z to P
  P <- 2 * pnorm(Z, lower.tail = F)
  P <- matrix(p.adjust(P, method = "BH"),
              nrow(P), ncol(P))
  
  tot.pcc.num <- p * (p - 1) / 2
  diff.pcc.num <- sum(P[upper.tri(P)] < conf_alpha)
  diff.pcc.ratio <- diff.pcc.num / tot.pcc.num
  
  if (is.na(topk)) {
    A <- 1 * (P < conf_alpha)
  }
  else {
    all.pval <- P[upper.tri(P)]
    pval.thresh <- sort(all.pval)[topk]
    A <- 1 * (P <= pval.thresh)
  }
  diag(A) <- 0
  
  diffPCC.res <- list(diff.pcc.num = diff.pcc.num, 
                      diff.ratio = diff.pcc.ratio, 
                      P.Mat = P, 
                      A = A)
  return(diffPCC.res)
}


## Classification Model ####
SRsamp <- function(N, n, group, group1 = NULL) {
  grp.f <- as.factor(group)
  grp.num <- nlevels(grp.f)
  grp.lev <- levels(grp.f)
  
  g.list <- list()
  for (i in c(1: grp.num)) {
    grpi.Idx <- which(group == grp.lev[i])
    N.grpi <- length(grpi.Idx)
    
    if (!is.null(group1)) {
      grp1.i <- group1[group == grp.lev[i]]
      grp1.i.f <- as.factor(grp1.i)
      grp1.i.lev <- levels(grp1.i.f)
      grp1.i.num <- nlevels(grp1.i.f)
    }
    else {
      grp1.i.num <- 0
    }
    
    subIdx.list <- list()
    if (grp1.i.num > 1) {
      grpi.list <- list()
      for (j in c(1: grp1.i.num)) {
        grpij.Idx <- grpi.Idx[grp1.i == grp1.i.lev[j]]
        N.grpij <- length(grpij.Idx)
        
        # disturb sample order
        grpij.Idx.disturb <- sample(grpij.Idx, N.grpij)
        subNum.ij <- round(N.grpij / n)
        sub_subIdx.list <- list()
        for (k in c(1: n)) {
          startIdx <- (k - 1) * subNum.ij + 1
          if (k < n) {
            endIdx <- k * subNum.ij
          }
          else {
            endIdx <- N.grpij
          }
          
          sub_subIdx.list[[k]] <- grpij.Idx.disturb[c(startIdx: endIdx)]
        }
        grpi.list[[j]] <- sub_subIdx.list
      }
      
      for (k in c(1: n)) {
        s.k <- c()
        for (j in c(1: grp1.i.num)) {
          sub_subIdxlist.j <- grpi.list[[j]]
          s.k <- append(s.k, sub_subIdxlist.j[[k]])
        }
        subIdx.list[[k]] <- s.k
      }
    }
    else {
      # disturb sample order
      grpi.Idx.disturb <- sample(grpi.Idx, N.grpi)
      subNum.i <- round(N.grpi / n)
      for (k in c(1: n)) {
        startIdx <- (k - 1) * subNum.i + 1
        if (k < n) {
          endIdx <- k * subNum.i
        }
        else {
          endIdx <- N.grpi
        }
        
        subIdx.list[[k]] <- grpi.Idx.disturb[c(startIdx: endIdx)]
      }
    }
    
    g.list[[i]] <- subIdx.list
  }
  
  S.list <- list()
  for (i in c(1: n)) {
    s.i <- c()
    for (g in c(1: grp.num)) {
      subIdxlist.g <- g.list[[g]]
      s.i <- append(s.i, subIdxlist.g[[i]])
    }
    S.list[[i]] <- s.i
  }
  
  return(S.list)
}

Safe_ClassTrain <- function(dat, group, method = method, tr_contrl, 
                            preProc, metric) {
  repeat {
    RF.mod <- tryCatch(
      expr = {
        model <- train(x = dat, y = group, method = method, 
                       trControl = tr_contrl, metric = metric, 
                       preProcess = preProc, ntree = 200, 
                       importance = T)
        return(model)
      }, 
      warning = function(w) {
        return(model)
      }
    )
    
    if (!is.null(RF.mod)) {
      break
    }
  }
  
  return(RF.mod)
}

CVforClas_Pred <- function(train_dat, test_dat, tr.group, ts.group, 
                           method = 'rf', doLog = F, resamp_method = 'cv', 
                           fold = 10, preproc = NULL, metric = NA, 
                           repeat_times = 100, show.progress = F) {
  if (doLog) {
    train_dat <- log2(train_dat)
    test_dat <- log2(test_dat)
  }
  
  trgrp.f <- as.factor(tr.group)
  trgrp.lev <- levels(trgrp.f)
  trgrp.num <- nlevels(trgrp.f)
  
  tsgrp.f <- as.factor(ts.group)
  tsgrp.lev <- levels(tsgrp.f)
  tsgrp.num <- nlevels(tsgrp.f)
  
  if (trgrp.num == 2) {
    sF <- twoClassSummary
    
    if (is.na(metric)) {
      metric <- 'ROC'
    }
  }
  else if (trgrp.num >= 3) {
    sF <- multiClassSummary
    
    if (is.na(metric)) {
      metric <- 'prAUC'
    }
  }
  
  tc <- trainControl(method = resamp_method, number = fold, 
                     summaryFunction = sF, 
                     classProbs = T, savePredictions = 'final')
  
  Preds <- matrix(0, nrow(test_dat), tsgrp.num)
  
  if (show.progress) {
    pb <- progress_bar$new(total = repeat_times, width = 60, 
                           format = "Classification models run [:bar] :percent in :elapsed", 
                           clear = F)
  }
  

  for (t in c(1: repeat_times)) {
    
    rf.model <- Safe_ClassTrain(dat = train_dat, group = trgrp.f, 
                                method = method, tr_contrl = tc, 
                                preProc = preproc, metric = metric)
    
    pred <- predict(rf.model, test_dat, type = 'prob')
    Preds <- Preds + pred
    
    if (show.progress) {
      pb$tick()
    }
    
  }
  
  meanPred <- Preds / repeat_times

  return(meanPred)
}

getClasMetric <- function(Pred.prob, group.GT, prob.lev = 0.5, ident_grp.idx = NA) {
  
  grp.f <- as.factor(group.GT)
  grp.lev <- levels(grp.f)
  grp.num <- nlevels(grp.f)
  
  metric.names <- c('Precision', 'Recall', 'F1-score', 
                    'Accuracy', 'AUROC')
  
  if (is.matrix(Pred.prob)) {
    if (ncol(Pred.prob) == grp.num) {
      EvalMetrics <- matrix(0, grp.num, length(metric.names))
      for (i in c(1: grp.num)) {
        pred.prob.i <- Pred.prob[, i]
        pred.i <- 1 * (pred.prob.i > prob.lev)
        GT.i <- 1 * (group.GT == grp.lev[i])
        
        TP.num <- sum(pred.i == 1 & GT.i == 1)
        TN.num <- sum(pred.i == 0 & GT.i == 0)
        FP.num <- sum(pred.i == 1 & GT.i == 0)
        FN.num <- sum(pred.i == 0 & GT.i == 1)
        
        if (TP.num == 0) {
          EvalMetrics[i, 1] <- 0
          EvalMetrics[i, 2] <- 0
        }
        else {
          EvalMetrics[i, 1] <- TP.num / (TP.num + FP.num)
          EvalMetrics[i, 2] <- TP.num / (TP.num + FN.num)
        }
        
        if (EvalMetrics[i, 1] == 0 & EvalMetrics[i, 2] == 0) {
          EvalMetrics[i, 3] <- 0
        }
        else {
          EvalMetrics[i, 3] <- 2 * EvalMetrics[i, 1] * EvalMetrics[i, 2] / (EvalMetrics[i, 1] + EvalMetrics[i, 2])
        }
        
        EvalMetrics[i, 4] <- (TP.num + TN.num) / length(GT.i)
        
        pred.res.i <- prediction(pred.prob.i, GT.i)
        EvalMetrics[i, 5] <- performance(pred.res.i, measure = 'auc')@y.values[[1]]
        
      }
      colnames(EvalMetrics) <- metric.names
    }
    else {
      stop('Not provide predictions of multiple groups!')
    }
  }
  else if (is.vector(Pred.prob) & is.numeric(ident_grp.idx)) {
    EvalMetrics <- rep(0, length(metric.names))
    pred <- 1 * (Pred.prob > prob.lev)
    GT <- 1 * (group.GT == grp.lev[ident_grp.idx])
    
    TP.num <- sum(pred == 1 & GT == 1)
    TN.num <- sum(pred == 0 & GT == 0)
    FP.num <- sum(pred == 1 & GT == 0)
    FN.num <- sum(pred == 0 & GT == 1)
    
    if (TP.num == 0) {
      EvalMetrics[1] <- 0
      EvalMetrics[2] <- 0
    }
    else {
      EvalMetrics[1] <- TP.num / (TP.num + FP.num)
      EvalMetrics[2] <- TP.num / (TP.num + FN.num)
    }
    
    if (EvalMetrics[1] == 0 & EvalMetrics[2] == 0) {
      EvalMetrics[3] <- 0
    }
    else {
      EvalMetrics[3] <- 2 * EvalMetrics[1] * EvalMetrics[2] / (EvalMetrics[1] + EvalMetrics[2])
    }
    
    EvalMetrics[4] <- (TP.num + TN.num) / length(GT)
    
    pred.res <- prediction(Pred.prob, GT)
    EvalMetrics[5] <- performance(pred.res, measure = 'auc')@y.values[[1]]
  }
  else {
    stop('Error!')
  }
  
  
  return(EvalMetrics)
}


## Cumulative Frequency of RSDs ####
getCumFreq <- function(rsd.vals, rsds = seq(0, 1, 0.01)) {
  tot.num <- length(rsd.vals)
  rsd.num <- length(rsds)
  
  cumFreq <- rep(0, rsd.num)
  for (i in c(1: rsd.num)) {
    freq.i <- sum(rsd.vals < rsds[i]) / tot.num
    cumFreq[i] <- freq.i
  }
  
  return(cumFreq)
}

getAUCforRSD <- function(dat.list, group = NULL, 
                         rsds = seq(0.01, 1, 0.01)) {
  dat.num <- length(dat.list)
  dat.name <- names(dat.list)
  
  N <- nrow(dat.list[[1]])
  if (is.null(group)) {
    group <- rep(1, N)
  }
  
  grp.f <- as.factor(group)
  grp.num <- nlevels(grp.f)
  grp.lev <- levels(grp.f)
  
  cumFreq.mat <- matrix(0, length(rsds), dat.num)
  rsdAUCs <- rep(0, dat.num)
  for (i in c(1: dat.num)) {
    dat.i <- dat.list[[i]]
    
    cumFreq.i <- rep(0, length(rsds))
    for (g in c(1: grp.num)) {
      dati.g <- dat.i[group == grp.lev[g], ]
      m.ig <- colMeans(dati.g)
      sd.ig <- apply(dati.g, MARGIN = 2, FUN = sd)
      rsd.ig <- sd.ig / m.ig
      cum.freq <- getCumFreq(rsd.ig, rsds = rsds)
      cumFreq.i <- cumFreq.i + cum.freq
    }
    cumFreq.i <- cumFreq.i / grp.num
    cumFreq.mat[, i] <- cumFreq.i
    
    rsd.num <- length(rsds)
    sum.area <- 0
    for (j in c(1: (rsd.num - 1))) {
      area.j <- (cumFreq.i[j] + cumFreq.i[j + 1]) * (rsds[j + 1] - rsds[j]) / 2
      sum.area <- sum.area + area.j
    }
    rsdAUCs[i] <- sum.area
  }
  
  colnames(cumFreq.mat) <- dat.name
  rsdAUCs <- as.data.frame(t(rsdAUCs))
  colnames(rsdAUCs) <- dat.name
  
  RSD.res <- list(cumFreqs = cumFreq.mat, 
                  rsdAUCs = rsdAUCs)
  return(RSD.res)
}

