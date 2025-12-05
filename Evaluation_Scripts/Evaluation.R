# Date: 20250618
# Author: Fanjing Guo
# Title: Evaluation

# FUNCTION ####
library(here)
source(here("LSCN_Function", "Funcs_LSCN_Algorithm.R"))
source(here("Evaluation_Scripts", "Funcs_Norm_Methods.R"))
source(here("Evaluation_Scripts", "Funcs_Evaluation.R"))
source(here("Simulation_Code", "Funcs_simulation.R"))

# Get data ####
## Data used for generating simulated data ####
dat <- read.csv(here("Data", "Diabet_data_proc_with_Info_new1.csv"), header = T)
group <- dat$Group
order <- dat$Order
dat <- dat[, -c(1: 3)]
Metabo <- colnames(dat)
dat <- as.matrix(dat)
N <- nrow(dat)
p <- ncol(dat)

dat.healthy <- dat[group == "Control Group", ]


## Real-world data ####
#### Dataset I. ST003484 ####
# Untargeted LC-MS (Not containing QCs)
dat <- read.csv(here("Data", "dat_ST003484_pos.csv"), header = T)
dat.title <- "ST003484"
N <- nrow(dat)
SampName <- dat$Samples
group <- dat$Factors
group <- gsub(' ', '', group)
group.info <- strsplit(group, '\\|')
g1.info <- sapply(group.info, "[", 1)
g1.info <- strsplit(g1.info, ':')
g1 <- sapply(g1.info, "[", 2)
g2.info <- sapply(group.info, "[", 2)
g2.info <- strsplit(g2.info, ':')
g2 <- as.numeric(sapply(g2.info, "[", 2))
g3.info <- sapply(group.info, "[", 3)
g3.info <- strsplit(g3.info, ':')
g3 <- sapply(g3.info, "[", 2)
g3.f <- as.factor(g3)
g3.lev <- levels(g3.f)
g3.num <- nlevels(g3.f)
for (i in c(1: g3.num)) {
  g3[g3 == g3.lev[i]] <- LETTERS[i]
}
groups <- data.frame(Condition = g1, 
                     Time = g2, 
                     Treatment = g3)
dat <- as.matrix(dat[, -c(1: 2)])
library(DMwR2)
dat[dat == 0] <- NA
na.num <- colSums(is.na(dat))
delmetIdx <- which(na.num / N > 0.5)
if (length(delmetIdx) != 0) {
  dat <- dat[, -delmetIdx]
}
dat <- as.data.frame(t(dat))
dat <- as.matrix(t(knnImputation(dat)))
p <- ncol(dat)
datIdx <- 'I'
containQC <- F

kegg_id <- read.csv(here("Data", "kegg_id_I.csv"), header = T)
MetName <- kegg_id$metabolite_name
kegg_id <- kegg_id$KEGG.ID


### Dataset II. MTBLS290 ####
# Untargeted LC-MS (Containing QCs)
dat <- read.csv(here("Data", "dat_MTBLS290.csv"), header = T)
dat.title <- 'MTBLS290'
g1 <- dat$Group
dat <- dat[g1 != 'BK', ]
g1 <- g1[g1 != 'BK']
N <- nrow(dat)
SampName <- dat$Sample
p.time <- strsplit(SampName, 't')
p.time <- sapply(p.time, 
                 FUN = function(x) x[2])
g2 <- p.time
# delete T1 and T7
T1_7.idx <- which(g2 == "1" | g2 == "7")
g2 <- as.numeric(g2[-T1_7.idx])
dat <- dat[-T1_7.idx, ]
N <- nrow(dat)
SampName <- SampName[-T1_7.idx]
g1 <- g1[-T1_7.idx]
sampname.split <- strsplit(SampName, 't')
groups <- data.frame(Disease = g1, 
                     TimePoint = g2)
dat <- as.matrix(dat[, -c(1: 2)])
dat[dat == 0] <- NA
library(DMwR2)
na.num <- colSums(is.na(dat))
delmetIdx <- which(na.num / N > 0.5)
if (length(delmetIdx) != 0) {
  dat <- dat[, -delmetIdx]
}
dat <- as.data.frame(t(dat))
dat <- as.matrix(t(knnImputation(dat)))
rownames(dat) <- SampName
MetName <- colnames(dat) # too long
p <- ncol(dat)
MetName <- paste0('f', c(1: p))

datIdx <- 'II'

containQC <- T
SampType <- g1

#### Dataset III. ST002178 ####
# Targeted LC-MS (Containing QCs)
dat <- read.csv(here("Data", "dat_ST002178_pos.csv"), header = T)
dat.title <- 'ST002178'
N <- nrow(dat)
SampName <- dat$Samples
group.info <- dat$Factors
group.info <- gsub(' ', '', group.info)
group.info <- strsplit(group.info, '\\|')
g1.info <- sapply(group.info, "[", 1)
g1.info <- strsplit(g1.info, ':')
g1 <- sapply(g1.info, "[", 2)
g1[grep('qualitycontrol', g1)] <- 'QC'
g2.info <- sapply(group.info, "[", 2)
g2.info <- strsplit(g2.info, ':')
g2 <- sapply(g2.info, "[", 2)
g2[g1 == 'QC'] <- NA
g2 <- gsub('-', '_', g2)
groups <- data.frame(Age = g1, 
                     Treatment = g2)
dat <- as.matrix(dat[, -c(1: 2)])
MetName <- colnames(dat)
p <- ncol(dat)

datIdx <- 'III'

containQC <- T
SampType <- g1

kegg_id <- read.csv(here("Data", "kegg_id_III.csv"), header = T)
MetName <- kegg_id$metabolite_name
kegg_id <- kegg_id$KEGG.ID

filtered_met <- read.csv(here("Data", "filtered_met_III.csv"), header = T)
filtered_met <- filtered_met$metabolite_name

filtered_met_kegg_id <- rep('', length(filtered_met))
for (i in c(1: length(filtered_met))) {
  filtered_met_kegg_id[i] <- kegg_id[MetName == filtered_met[i]]
}
kegg_id <- filtered_met_kegg_id
MetName <- filtered_met

### Dataset IV. Tea DMSO ####
# NMR (Not containing QCs)
dat <- read.csv(here("Data", "teadata_DMSO.csv"), header = T)
dat.title <- "teadata_DMSO"
groups <- data.frame(Variety = dat$Variety, 
                     Site = dat$Site, 
                     Year = dat$Year)
SampName <- dat$Name
dat <- as.matrix(dat[, -c(1: 4)])
N <- nrow(dat)
rownames(dat) <- paste0('s', c(1: N))
MetName <- colnames(dat)
delIdx <- which(MetName == 'uk.51052')
dat <- dat[, -delIdx]
p <- ncol(dat)
MetName <- colnames(dat)

datIdx <- 'IV'

containQC <- F

## Normalization ####
# LSCN
if (dat.title == 'teadata_DMSO') {
  sq <- 0.9
} else {
  sq <- NA
}
dat.nor <- LSCN(dat, doLog = F, ita0 = 0.4, max_ite = 50, 
                simi.quant = sq, print.txt = T)
dat.LSCN.real <- dat.nor$dat.norm
coef.LSCN <- dat.nor$NORM.COEF
rownames(dat.LSCN.real) <- SampName
colnames(dat.LSCN.real) <- MetName

# CSN
dat.nor <- CSN_normalize(dat)
dat.CSN.real <- dat.nor$dat.norm
coef.CSN <- dat.nor$norm.coef
rownames(dat.LSCN.real) <- SampName
colnames(dat.LSCN.real) <- MetName

# L2N
dat.nor <- L2_Normalize(dat)
dat.L2N.real <- dat.nor$dat.norm
coef.L2N <- dat.nor$norm.coef
rownames(dat.L2N.real) <- SampName
colnames(dat.L2N.real) <- MetName

# PQN
dat.nor <- PQN_normalize(dat)
dat.PQN.real <- dat.nor$dat.norm
coef.PQN <- dat.nor$norm.coef
rownames(dat.PQN.real) <- SampName
colnames(dat.PQN.real) <- MetName

# Quantile
dat.nor <- Quantile_normalize(dat)
dat.QT.real <- dat.nor$dat.norm
coef.QT <- dat.QT.real / dat
rownames(dat.QT.real) <- SampName
colnames(dat.QT.real) <- MetName


med.bef <- median(dat)
med.LSCN <- median(dat.LSCN.real)
med.CSN <- median(dat.CSN.real)
med.L2N <- median(dat.L2N.real)
med.PQN <- median(dat.PQN.real)
med.QT <- median(dat.QT.real)

dat.LSCN.real <- dat.LSCN.real * (med.bef / med.LSCN)
dat.CSN.real <- dat.CSN.real * (med.bef / med.CSN)
dat.L2N.real <- dat.L2N.real * (med.bef / med.L2N)
dat.PQN.real <- dat.PQN.real * (med.bef / med.PQN)
dat.QT.real <- dat.QT.real * (med.bef / med.QT)



## Get subject samples ####
if (containQC) {
  sampIdx.nQC <- which(SampType != 'QC')
  N <- length(sampIdx.nQC)
  dat.list <- list(Unnorm = dat[sampIdx.nQC, ], 
                   CSN = dat.CSN.real[sampIdx.nQC, ], 
                   L2N = dat.L2N.real[sampIdx.nQC, ], 
                   PQN = dat.PQN.real[sampIdx.nQC, ], 
                   QT = dat.QT.real[sampIdx.nQC, ], 
                   LSCN = dat.LSCN.real[sampIdx.nQC, ])
  
  groups <- groups[sampIdx.nQC, ]
} else {
  sampIdx.nQC <- c(1: N)
  dat.list <- list(Unnorm = dat, 
                   CSN = dat.CSN.real, 
                   L2N = dat.L2N.real, 
                   PQN = dat.PQN.real, 
                   QT = dat.QT.real, 
                   LSCN = dat.LSCN.real)
}

## Determine the major group and minor group ####
Group_levels <- rep('', ncol(groups))
for (i in c(1: ncol(groups))) {
  group.i <- groups[, i]
  grpi.fac <- as.factor(group.i)
  grpi.lev <- levels(grpi.fac)
  Group_levels[i] <- paste0(grpi.lev, collapse = ',')
}

Fvals <- matrix(0, ncol(groups), length(dat.list))
for (i in c(1: ncol(groups))) {
  for (j in c(1: length(dat.list))) {
    dat.j <- dat.list[[j]]
    dat.j.log <- log2(dat.j)
    sd.j <- apply(dat.j.log, MARGIN = 2, FUN = sd)
    if (any(sd.j == 0)) {
      dat.j.log <- dat.j.log[, sd.j != 0]
    }
    dat.pca <- pca(dat.j.log, center = T, scale = T)
    dat.pcs <- dat.pca$variates$X
    group <- groups[, i]
    if (sum(is.na(group)) != 0) {
      dat.pcs <- dat.pcs[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    grpi.fac <- as.factor(group)
    grpi.lev <- levels(grpi.fac)
    grpi.num <- nlevels(grpi.fac)
    
    var_zero <- F
    for (g in c(1: grpi.num)) {
      dat.j.grpg <- dat.pcs[group == grpi.lev[g], 1]
      sd.j.grpg <- sd(dat.j.grpg)
      
      if (sd.j.grpg == 0) {
        var_zero <- T
        Fvals[i, j] <- NA
        break
      }
    }
    
    if (!var_zero) {
      dat1 <- data.frame(grp = grpi.fac, 
                         val = dat.pcs[, 1])
      
      levT <- leveneTest(val ~ grp, data = dat1)
      levTPval <- levT$`Pr(>F)`[1]
      
      if (levTPval < 0.05){
        varEqual <- F
      } else{
        varEqual <- T
      }
      
      aov.res <- oneway.test(val ~ grp, data = dat1, 
                             var.equal = varEqual)
      Fvals[i, j] <- aov.res$statistic
    }
  }
}
med_Fval <- apply(Fvals, MARGIN = 1, FUN = median)
F_rank <- apply(- Fvals, MARGIN = 2, FUN = rank)
med_Fr <- apply(F_rank, MARGIN = 1, FUN = median)

grp.df <- data.frame(ID = c(1: ncol(groups)), 
                     Group = colnames(groups), 
                     Levels = Group_levels, 
                     medFval = med_Fval, 
                     medRank = med_Fr)
print(grp.df)

if (dat.title != 'teadata_DMSO') {
  majgrpId <- 1
  majGroup <- groups[, majgrpId]
  
  mingrpId <- 2
  minGroup <- groups[, mingrpId]
} else {
  majgrpId <- 1
  majGroup <- groups[, majgrpId]
  
  mingrpId <- 3
  minGroup <- groups[, mingrpId]
}


# Eval 1 - Influence of Normalization on Heterogeneity ####

## a. (Simulated) Intar-group LSS & Heterogeneity Similarity ####
### I. Individual differences ####
library(mixOmics)
rep.time <- 100
n.ratio <- seq(0.05, 0.45, 0.1)

HeterSI.simu <- matrix(0, rep.time, length(n.ratio))
HeterSI.LSCN <- matrix(0, rep.time, length(n.ratio))
HeterSI.CSN <- matrix(0, rep.time, length(n.ratio))
HeterSI.L2N <- matrix(0, rep.time, length(n.ratio))
HeterSI.PQN <- matrix(0, rep.time, length(n.ratio))
HeterSI.QT <- matrix(0, rep.time, length(n.ratio))

structG.si.simu <- matrix(0, rep.time, length(n.ratio))
structG.si.LSCN <- matrix(0, rep.time, length(n.ratio))
structG.si.CSN <- matrix(0, rep.time, length(n.ratio))
structG.si.L2N <- matrix(0, rep.time, length(n.ratio))
structG.si.PQN <- matrix(0, rep.time, length(n.ratio))
structG.si.QT <- matrix(0, rep.time, length(n.ratio))

for (i in c(1: length(n.ratio))) {
  
  heters.si.simu <- rep(0, rep.time)
  heters.si.LSCN <- rep(0, rep.time)
  heters.si.CSN <- rep(0, rep.time)
  heters.si.L2N <- rep(0, rep.time)
  heters.si.PQN <- rep(0, rep.time)
  heters.si.QT <- rep(0, rep.time)
  
  structSIG.simu <- rep(0, rep.time)
  structSIG.LSCN <- rep(0, rep.time)
  structSIG.CSN <- rep(0, rep.time)
  structSIG.L2N <- rep(0, rep.time)
  structSIG.PQN <- rep(0, rep.time)
  structSIG.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, 
                               indiv.sd.ratio = n.ratio[i])
    dat_simu <- dat_simu.res$SimuDat
    grp_simu <- dat_simu.res$Para$group
    dat.GT <- dat_simu.res$InitDat
    
    # Unnorm
    heters.si.simu[k] <- Heter.Similar(dat.GT, dat_simu, group = grp_simu)
    structSIG.simu[k] <- struct.similar(dat.GT, dat_simu, group = grp_simu)
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    dat.nor.CSN <- dat.nor$dat.norm
    heters.si.CSN[k] <- Heter.Similar(dat.GT, dat.nor.CSN, group = grp_simu)
    structSIG.CSN[k] <- struct.similar(dat.GT, dat.nor.CSN, group = grp_simu)
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    dat.nor.L2N <- dat.nor$dat.norm
    heters.si.L2N[k] <- Heter.Similar(dat.GT, dat.nor.L2N, group = grp_simu)
    structSIG.L2N[k] <- struct.similar(dat.GT, dat.nor.L2N, group = grp_simu)
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    dat.nor.PQN <- dat.nor$dat.norm
    heters.si.PQN[k] <- Heter.Similar(dat.GT, dat.nor.PQN, group = grp_simu)
    structSIG.PQN[k] <- struct.similar(dat.GT, dat.nor.PQN, group = grp_simu)
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    heters.si.QT[k] <- Heter.Similar(dat.GT, dat.nor.QT, group = grp_simu)
    structSIG.QT[k] <- struct.similar(dat.GT, dat.nor.QT, group = grp_simu)
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, max_ite = 30, 
                    print.txt = F)
    dat.nor.LSCN <- dat.nor$dat.norm
    heters.si.LSCN[k] <- Heter.Similar(dat.GT, dat.nor.LSCN, group = grp_simu)
    structSIG.LSCN[k] <- struct.similar(dat.GT, dat.nor.LSCN, group = grp_simu)
  }
  
  HeterSI.simu[, i] <- heters.si.simu
  HeterSI.LSCN[, i] <- heters.si.LSCN
  HeterSI.CSN[, i] <- heters.si.CSN
  HeterSI.L2N[, i] <- heters.si.L2N
  HeterSI.PQN[, i] <- heters.si.PQN
  HeterSI.QT[, i] <- heters.si.QT
  
  structG.si.simu[, i] <- structSIG.simu
  structG.si.LSCN[, i] <- structSIG.LSCN
  structG.si.CSN[, i] <- structSIG.CSN
  structG.si.L2N[, i] <- structSIG.L2N
  structG.si.PQN[, i] <- structSIG.PQN
  structG.si.QT[, i] <- structSIG.QT
  
  cat(i, '\n')
}

# Degree Heterogeneity Similarity (DHS)
vals <- c(as.vector(HeterSI.LSCN), as.vector(HeterSI.CSN), 
          as.vector(HeterSI.L2N), as.vector(HeterSI.PQN), 
          as.vector(HeterSI.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(n.ratio)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

nr <- rep(rep(n.ratio, each = rep.time), 5)

dat.heterSI <- data.frame(dattype = dat.type, 
                          NR = nr, 
                          heterSI = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)

ggplot() + 
  geom_smooth(data = dat.heterSI, aes(x = NR, y = heterSI,
                                      color = dattype, fill = dattype, 
                                      linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = n.ratio) + 
  theme_bw() + 
  xlab(NULL) + ylab('The DHS Index') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

# Local Sample Structure Similarity (LSS)
vals <- c(as.vector(structG.si.LSCN), as.vector(structG.si.CSN), 
          as.vector(structG.si.L2N), as.vector(structG.si.PQN), 
          as.vector(structG.si.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(n.ratio)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

nr <- rep(rep(n.ratio, each = rep.time), 5)

dat.structSI <- data.frame(dattype = dat.type, 
                           NR = nr, 
                           structSI = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.structSI, aes(x = NR, y = structSI,
                                       color = dattype, fill = dattype, 
                                       linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = n.ratio) + 
  theme_bw() + 
  xlab('Individual Differences (r)') + 
  ylab('The LSS Index') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


### II. Intra-group Heterogeneity ####
rep.time <- 100
sg.num <- seq(2, 6, 1)

HeterSI.simu <- matrix(0, rep.time, length(sg.num))
HeterSI.LSCN <- matrix(0, rep.time, length(sg.num))
HeterSI.CSN <- matrix(0, rep.time, length(sg.num))
HeterSI.L2N <- matrix(0, rep.time, length(sg.num))
HeterSI.PQN <- matrix(0, rep.time, length(sg.num))
HeterSI.QT <- matrix(0, rep.time, length(sg.num))

structG.si.simu <- matrix(0, rep.time, length(sg.num))
structG.si.LSCN <- matrix(0, rep.time, length(sg.num))
structG.si.CSN <- matrix(0, rep.time, length(sg.num))
structG.si.L2N <- matrix(0, rep.time, length(sg.num))
structG.si.PQN <- matrix(0, rep.time, length(sg.num))
structG.si.QT <- matrix(0, rep.time, length(sg.num))

for (i in c(1: length(sg.num))) {
  
  heters.si.simu <- rep(0, rep.time)
  heters.si.LSCN <- rep(0, rep.time)
  heters.si.CSN <- rep(0, rep.time)
  heters.si.L2N <- rep(0, rep.time)
  heters.si.PQN <- rep(0, rep.time)
  heters.si.QT <- rep(0, rep.time)
  
  structSIG.simu <- rep(0, rep.time)
  structSIG.LSCN <- rep(0, rep.time)
  structSIG.CSN <- rep(0, rep.time)
  structSIG.L2N <- rep(0, rep.time)
  structSIG.PQN <- rep(0, rep.time)
  structSIG.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, 
                               subgrpNum = sg.num[i])
    dat_simu <- dat_simu.res$SimuDat
    grp_simu <- dat_simu.res$Para$group
    dat.GT <- dat_simu.res$InitDat
    
    # Unnorm
    heters.si.simu[k] <- Heter.Similar(dat.GT, dat_simu, group = grp_simu)
    structSIG.simu[k] <- struct.similar(dat.GT, dat_simu, group = grp_simu)
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    dat.nor.CSN <- dat.nor$dat.norm
    heters.si.CSN[k] <- Heter.Similar(dat.GT, dat.nor.CSN, group = grp_simu)
    structSIG.CSN[k] <- struct.similar(dat.GT, dat.nor.CSN, group = grp_simu)
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    dat.nor.L2N <- dat.nor$dat.norm
    heters.si.L2N[k] <- Heter.Similar(dat.GT, dat.nor.L2N, group = grp_simu)
    structSIG.L2N[k] <- struct.similar(dat.GT, dat.nor.L2N, group = grp_simu)
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    dat.nor.PQN <- dat.nor$dat.norm
    heters.si.PQN[k] <- Heter.Similar(dat.GT, dat.nor.PQN, group = grp_simu)
    structSIG.PQN[k] <- struct.similar(dat.GT, dat.nor.PQN, group = grp_simu)
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    heters.si.QT[k] <- Heter.Similar(dat.GT, dat.nor.QT, group = grp_simu)
    structSIG.QT[k] <- struct.similar(dat.GT, dat.nor.QT, group = grp_simu)
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, max_ite = 30, 
                    print.txt = F)
    dat.nor.LSCN <- dat.nor$dat.norm
    heters.si.LSCN[k] <- Heter.Similar(dat.GT, dat.nor.LSCN, group = grp_simu)
    structSIG.LSCN[k] <- struct.similar(dat.GT, dat.nor.LSCN, group = grp_simu)
  }
  
  HeterSI.simu[, i] <- heters.si.simu
  HeterSI.LSCN[, i] <- heters.si.LSCN
  HeterSI.CSN[, i] <- heters.si.CSN
  HeterSI.L2N[, i] <- heters.si.L2N
  HeterSI.PQN[, i] <- heters.si.PQN
  HeterSI.QT[, i] <- heters.si.QT
  
  structG.si.simu[, i] <- structSIG.simu
  structG.si.LSCN[, i] <- structSIG.LSCN
  structG.si.CSN[, i] <- structSIG.CSN
  structG.si.L2N[, i] <- structSIG.L2N
  structG.si.PQN[, i] <- structSIG.PQN
  structG.si.QT[, i] <- structSIG.QT
  
  cat(i, '\n')
}

# Degree Heterogeneity Similarity (DHS)
vals <- c(as.vector(HeterSI.LSCN), as.vector(HeterSI.CSN), 
          as.vector(HeterSI.L2N), as.vector(HeterSI.PQN), 
          as.vector(HeterSI.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(sg.num)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

sgn <- rep(rep(sg.num, each = rep.time), 5)

dat.heterSI <- data.frame(dattype = dat.type, 
                          SGN = sgn, 
                          heterSI = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)

ggplot() + 
  geom_smooth(data = dat.heterSI, aes(x = SGN, y = heterSI,
                                      color = dattype, fill = dattype, 
                                      linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = sg.num) + 
  theme_bw() + 
  xlab(NULL) + ylab('The DHS Index') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'top', 
        legend.text = element_text(size = 18, family = 'sans', face = 'bold'), 
        legend.title = element_blank(), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

# Local Sample Structure Similarity (LSS)
vals <- c(as.vector(structG.si.LSCN), as.vector(structG.si.CSN), 
          as.vector(structG.si.L2N), as.vector(structG.si.PQN), 
          as.vector(structG.si.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(sg.num)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

sgn <- rep(rep(sg.num, each = rep.time), 5)

dat.structSI <- data.frame(dattype = dat.type, 
                           SGN = sgn, 
                           structSI = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.structSI, aes(x = SGN, y = structSI,
                                       color = dattype, fill = dattype, 
                                       linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = sg.num) + 
  theme_bw() + 
  xlab('Intra-group Heterogeneity (n)') + 
  ylab('The LSS Index') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


### III. Dilution Effects ####
rep.time <- 100
sv.m <- seq(10, 90, 20)

HeterSI.simu <- matrix(0, rep.time, length(sv.m))
HeterSI.LSCN <- matrix(0, rep.time, length(sv.m))
HeterSI.CSN <- matrix(0, rep.time, length(sv.m))
HeterSI.L2N <- matrix(0, rep.time, length(sv.m))
HeterSI.PQN <- matrix(0, rep.time, length(sv.m))
HeterSI.QT <- matrix(0, rep.time, length(sv.m))

structG.si.simu <- matrix(0, rep.time, length(sv.m))
structG.si.LSCN <- matrix(0, rep.time, length(sv.m))
structG.si.CSN <- matrix(0, rep.time, length(sv.m))
structG.si.L2N <- matrix(0, rep.time, length(sv.m))
structG.si.PQN <- matrix(0, rep.time, length(sv.m))
structG.si.QT <- matrix(0, rep.time, length(sv.m))

for (i in c(1: length(sv.m))) {
  
  heters.si.simu <- rep(0, rep.time)
  heters.si.LSCN <- rep(0, rep.time)
  heters.si.CSN <- rep(0, rep.time)
  heters.si.L2N <- rep(0, rep.time)
  heters.si.PQN <- rep(0, rep.time)
  heters.si.QT <- rep(0, rep.time)
  
  structSIG.simu <- rep(0, rep.time)
  structSIG.LSCN <- rep(0, rep.time)
  structSIG.CSN <- rep(0, rep.time)
  structSIG.L2N <- rep(0, rep.time)
  structSIG.PQN <- rep(0, rep.time)
  structSIG.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, 
                               diltF.m = sv.m[i])
    dat_simu <- dat_simu.res$SimuDat
    grp_simu <- dat_simu.res$Para$group
    dat.GT <- dat_simu.res$InitDat
    
    # Unnorm
    heters.si.simu[k] <- Heter.Similar(dat_simu, dat.GT, group = grp_simu)
    structSIG.simu[k] <- struct.similar(dat.GT, dat_simu, group = grp_simu)
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    dat.nor.CSN <- dat.nor$dat.norm
    heters.si.CSN[k] <- Heter.Similar(dat.GT, dat.nor.CSN, group = grp_simu)
    structSIG.CSN[k] <- struct.similar(dat.GT, dat.nor.CSN, group = grp_simu)
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    dat.nor.L2N <- dat.nor$dat.norm
    heters.si.L2N[k] <- Heter.Similar(dat.GT, dat.nor.L2N, group = grp_simu)
    structSIG.L2N[k] <- struct.similar(dat.GT, dat.nor.L2N, group = grp_simu)
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    dat.nor.PQN <- dat.nor$dat.norm
    heters.si.PQN[k] <- Heter.Similar(dat.GT, dat.nor.PQN, group = grp_simu)
    structSIG.PQN[k] <- struct.similar(dat.GT, dat.nor.PQN, group = grp_simu)
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    heters.si.QT[k] <- Heter.Similar(dat.GT, dat.nor.QT, group = grp_simu)
    structSIG.QT[k] <- struct.similar(dat.GT, dat.nor.QT, group = grp_simu)
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, max_ite = 30, 
                    print.txt = F)
    dat.nor.LSCN <- dat.nor$dat.norm
    heters.si.LSCN[k] <- Heter.Similar(dat.GT, dat.nor.LSCN, group = grp_simu)
    structSIG.LSCN[k] <- struct.similar(dat.GT, dat.nor.LSCN, group = grp_simu)
  }
  
  HeterSI.simu[, i] <- heters.si.simu
  HeterSI.LSCN[, i] <- heters.si.LSCN
  HeterSI.CSN[, i] <- heters.si.CSN
  HeterSI.L2N[, i] <- heters.si.L2N
  HeterSI.PQN[, i] <- heters.si.PQN
  HeterSI.QT[, i] <- heters.si.QT
  
  structG.si.simu[, i] <- structSIG.simu
  structG.si.LSCN[, i] <- structSIG.LSCN
  structG.si.CSN[, i] <- structSIG.CSN
  structG.si.L2N[, i] <- structSIG.L2N
  structG.si.PQN[, i] <- structSIG.PQN
  structG.si.QT[, i] <- structSIG.QT
  
  cat(i, '\n')
}

# Degree Heterogeneity Similarity (DHS)
vals <- c(as.vector(HeterSI.LSCN), as.vector(HeterSI.CSN), 
          as.vector(HeterSI.L2N), as.vector(HeterSI.PQN), 
          as.vector(HeterSI.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(sv.m)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

svms <- rep(rep(sv.m, each = rep.time), 5)

dat.heterSI <- data.frame(dattype = dat.type, 
                          SVM = svms, 
                          heterSI = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)

ggplot() + 
  geom_smooth(data = dat.heterSI, aes(x = SVM, y = heterSI,
                                      color = dattype, fill = dattype, 
                                      linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = sv.m) + 
  theme_bw() + 
  xlab(NULL) + ylab('The DHS Index') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

# Local Sample Structure Similarity (LSS)
vals <- c(as.vector(structG.si.LSCN), as.vector(structG.si.CSN), 
          as.vector(structG.si.L2N), as.vector(structG.si.PQN), 
          as.vector(structG.si.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(sv.m)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

svms <- rep(rep(sv.m, each = rep.time), 5)

dat.structSI <- data.frame(dattype = dat.type, 
                           SVM = svms, 
                           structSI = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.structSI, aes(x = SVM, y = structSI,
                                       color = dattype, fill = dattype, 
                                       linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = sv.m) + 
  theme_bw() + 
  xlab('Dilution Effect (m)') + 
  ylab('The LSS Index') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


## b. (MTBLS290) Observation of QC samples ####
### I. Scatter plot of QCs ####
# Unnorm
pca.bef <- getPCscore(dat, result = 'pca')
pcs.bef <- pca.bef$score[, c(1: 2)]
expl.bef <- pca.bef$explVar[c(1: 2)]

# LSCN
pca.LSCN <- getPCscore(dat.LSCN.real, result = 'pca')
pcs.LSCN <- pca.LSCN$score[, c(1: 2)]
pcs.LSCN <- - pcs.LSCN
expl.LSCN <- pca.LSCN$explVar[c(1: 2)]

# CSN
pca.CSN <- getPCscore(dat.CSN.real, result = 'pca')
pcs.CSN <- pca.CSN$score[, c(1: 2)]
expl.CSN <- pca.CSN$explVar[c(1: 2)]

# L2N
pca.L2N <- getPCscore(dat.L2N.real, result = 'pca')
pcs.L2N <- pca.L2N$score[, c(1: 2)]
pcs.L2N <- - pcs.L2N
expl.L2N <- pca.L2N$explVar[c(1: 2)]

# PQN
pca.PQN <- getPCscore(dat.PQN.real, result = 'pca')
pcs.PQN <- pca.PQN$score[, c(1: 2)]
pcs.PQN <- - pcs.PQN
expl.PQN <- pca.PQN$explVar[c(1: 2)]

# QT
pca.QT <- getPCscore(dat.QT.real, result = 'pca')
pcs.QT <- pca.QT$score[, c(1: 2)]
expl.QT <- pca.QT$explVar[c(1: 2)]

# plot
all.score <- data.frame(s1 = c(pcs.bef[, 1], pcs.LSCN[, 1], 
                               pcs.CSN[, 1], pcs.L2N[, 1], 
                               pcs.PQN[, 1], pcs.QT[, 1]), 
                        s2 = c(pcs.bef[, 2], pcs.LSCN[, 2], 
                               pcs.CSN[, 2], pcs.L2N[, 2], 
                               pcs.PQN[, 2], pcs.QT[, 2]))

x.range <- c(min(all.score$s1) - 1, max(all.score$s1) + 1)
y.range <- c(min(all.score$s2) - 1, max(all.score$s2) + 1)


expl.bef1 <- format(round(expl.bef[1] * 100, 1), nsmall = 1)
expl.bef2 <- format(round(expl.bef[2] * 100, 1), nsmall = 1)
PCA_scatter_withQC(pcs.bef, group = SampType,
                   xlim = x.range, ylim = y.range, 
                   x_label = paste0('PC1: ', expl.bef1, '%'), 
                   y_label = paste0('PC2: ', expl.bef2, '%'), 
                   title = 'Unnorm')

expl.LSCN1 <- format(round(expl.LSCN[1] * 100, 1), nsmall = 1)
expl.LSCN2 <- format(round(expl.LSCN[2] * 100, 1), nsmall = 1)
PCA_scatter_withQC(pcs.LSCN, group = SampType,
                   xlim = x.range, ylim = y.range, 
                   x_label = paste0('PC1: ', expl.LSCN1, '%'), 
                   y_label = paste0('PC2: ', expl.LSCN2, '%'), 
                   title = 'LSCN')

expl.CSN1 <- format(round(expl.CSN[1] * 100, 1), nsmall = 1)
expl.CSN2 <- format(round(expl.CSN[2] * 100, 1), nsmall = 1)
PCA_scatter_withQC(pcs.CSN, group = SampType,
                   xlim = x.range, ylim = y.range, 
                   x_label = paste0('PC1: ', expl.CSN1, '%'), 
                   y_label = paste0('PC2: ', expl.CSN2, '%'), 
                   title = 'CSN')

expl.L2N1 <- format(round(expl.L2N[1] * 100, 1), nsmall = 1)
expl.L2N2 <- format(round(expl.L2N[2] * 100, 1), nsmall = 1)
PCA_scatter_withQC(pcs.L2N, group = SampType,
                   xlim = x.range, ylim = y.range, 
                   x_label = paste0('PC1: ', expl.L2N1, '%'), 
                   y_label = paste0('PC2: ', expl.L2N2, '%'), 
                   title = 'L2-norm')

expl.PQN1 <- format(round(expl.PQN[1] * 100, 1), nsmall = 1)
expl.PQN2 <- format(round(expl.PQN[2] * 100, 1), nsmall = 1)
PCA_scatter_withQC(pcs.PQN, group = SampType, 
                   xlim = x.range, ylim = y.range, 
                   x_label = paste0('PC1: ', expl.PQN1, '%'), 
                   y_label = paste0('PC2: ', expl.PQN2, '%'), 
                   title = 'PQN')

expl.QT1 <- format(round(expl.QT[1] * 100, 1), nsmall = 1)
expl.QT2 <- format(round(expl.QT[2] * 100, 1), nsmall = 1)
PCA_scatter_withQC(pcs.QT, group = SampType, 
                   xlim = x.range, ylim = y.range, 
                   x_label = paste0('PC1: ', expl.QT1, '%'), 
                   y_label = paste0('PC2: ', expl.QT2, '%'), 
                   title = 'Quantile')


### II. Euclidean distances of QCs ####
# QC samples in original space
dat.sca <- scale(dat, center = T, scale = T)
dat.LSCN.sca <- scale(dat.LSCN.real, center = T, scale = T)
dat.CSN.sca <- scale(dat.CSN.real, center = T, scale = T)
dat.L2N.sca <- scale(dat.L2N.real, center = T, scale = T)
dat.PQN.sca <- scale(dat.PQN.real, center = T, scale = T)
dat.QT.sca <- scale(dat.QT.real, center = T, scale = T)
QCdat.list <- list(Unnorm = dat.sca[-sampIdx.nQC, ], 
                   CSN = dat.CSN.sca[-sampIdx.nQC, ], 
                   L2N = dat.L2N.sca[-sampIdx.nQC, ], 
                   PQN = dat.PQN.sca[-sampIdx.nQC, ], 
                   QT = dat.QT.sca[-sampIdx.nQC, ], 
                   LSCN = dat.LSCN.sca[-sampIdx.nQC, ])

dat.num <- length(QCdat.list)
dat.name <- names(QCdat.list)
N.QC <- nrow(dat[-sampIdx.nQC, ])
tot.QCdist.num <- N.QC * (N.QC - 1) / 2

QCdist <- matrix(0, tot.QCdist.num, dat.num)
for (i in c(1: dat.num)) {
  QCdat.i <- QCdat.list[[i]]
  QCdist.i.mat <- as.matrix(dist(QCdat.i))
  QCdist[, i] <- QCdist.i.mat[upper.tri(QCdist.i.mat)]
}
apply(QCdist, MARGIN = 2, FUN = median)

## c. (MTBLS290) Network construction ####
library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)

set_col <- c("#0066FF", "#3F59D8", "#7F4CB2", "#BF3F8C", "#FF3366")
set_col1 <- c("#0066FF", "#7F4CB2", "#FF3366")
set_col2 <- c("#3F59D8", "#BF3F8C")

if (p < 100) {
  layout_method <- 'kk'
} else {
  layout_method <- 'fr'
}

pz <- 4
k_nn <- 3 # 0617

dat.num <- length(dat.list)
dat.name <- names(dat.list)

minGroup_1 <- minGroup[!is.na(minGroup)]
modls <- matrix(0, 3, dat.num)
for (i in c(1: dat.num)) {
  
  # draw network
  # all samples
  dat.i <- dat.list[[i]][!is.na(minGroup), ]
  A.i <- getSampNetwork(dat.i, k = k_nn, method = 'euclidean', 
                        doLog = F, reduce_dim = T)
  graph.i <- graph_from_adjacency_matrix(A.i, mode = 'undirected', 
                                         weighted = T, diag = FALSE)
  
  minGroup_1 <- minGroup[!is.na(minGroup)]
  p1 <- ggraph(graph = graph.i, layout = 'kk') +
    geom_edge_link() +
    geom_node_point(aes(color = as.factor(minGroup_1)),
                    size = pz, stroke = 2) + 
    scale_color_manual(values = set_col) +
    theme_graph() +
    theme(legend.position = 'none', 
          legend.text = element_text(family = 'sans', size = 18),
          legend.title = element_text(family = 'sans', size = 20),
          plot.title = element_text(size = 22, family = 'sans',
                                    face = 'bold', hjust = 0.5)) +
    labs(color = 'Disease State',
         title = dat.name[i])
  
  modls[1, i] <- modularity(graph.i, minGroup_1)
  
  # sample 2, 4, 6
  samp.idx <- which(minGroup_1 == 2 | minGroup_1 == 4 | minGroup_1 == 6)
  A.i <- getSampNetwork(dat.i[samp.idx, ], 
                        k = k_nn, method = 'euclidean', 
                        doLog = F, reduce_dim = T)
  graph.i <- graph_from_adjacency_matrix(A.i, mode = 'undirected', 
                                         weighted = T, diag = FALSE)
  
  mingrp1 <- minGroup_1[samp.idx]
  p2 <- ggraph(graph = graph.i, layout = 'kk') +
    geom_edge_link() +
    geom_node_point(aes(color = as.factor(mingrp1)),
                    size = pz, stroke = 2) + 
    scale_color_manual(values = set_col1) +
    theme_graph() +
    theme(legend.position = 'none', 
          plot.title = element_text(size = 22, family = 'sans',
                                    face = 'bold', hjust = 0.5)) +
    labs(title = '2,4,6')
  
  modls[2, i] <- modularity(graph.i, mingrp1)
  
  # sample 3, 5
  samp.idx <- which(minGroup_1 == 3 | minGroup_1 == 5)
  A.i <- getSampNetwork(dat.i[samp.idx, ], 
                        k = k_nn, method = 'euclidean', 
                        doLog = F, reduce_dim = T)
  graph.i <- graph_from_adjacency_matrix(A.i, mode = 'undirected', 
                                         weighted = T, diag = FALSE)
  
  mingrp2 <- minGroup_1[samp.idx]
  p3 <- ggraph(graph = graph.i, layout = 'kk') +
    geom_edge_link() +
    geom_node_point(aes(color = as.factor(mingrp2)),
                    size = pz, stroke = 2) + 
    scale_color_manual(values = set_col2) +
    theme_graph() +
    theme(legend.position = 'none', 
          plot.title = element_text(size = 22, family = 'sans',
                                    face = 'bold', hjust = 0.5)) +
    labs(title = '3,5')
  
  modls[3, i] <- modularity(graph.i, mingrp2)
  
  combined_plot <- p1 + (p2 / p3)
  
  print(combined_plot)
  
}


## d. (MTBLS290) ANOVA of patient samples ####
dat.num <- length(dat.list)
time_point <- as.factor(minGroup[majGroup == 'Patient'])
residual.all <- matrix(0, p, dat.num)

for (i in c(1: dat.num)) {
  dat.i <- dat.list[[i]]
  dat.i.p <- dat.i[majGroup == 'Patient', ]
  dat.i.p <- log2(dat.i.p)
  N.i.p <- nrow(dat.i.p)
  
  for (j in c(1: p)) {
    anova_result <- aov(dat.i.p[, j] ~ time_point)
    resi.r <- summary(anova_result)[[1]][2, 'Sum Sq'] / (N.i.p - 1) / var(dat.i.p[, j])
    residual.all[j, i] <- resi.r
  }
}

colMeans(residual.all)

# plot
library(ggbreak)
vals <- colMeans(residual.all)
dat.type <- factor(dat.name, 
                   levels = dat.name)
dat.res <- data.frame(datty = dat.type, 
                      residual = vals)

set_col <- c('gray', "#FF66FF", "#FFCC00", "#0066FF", "#00CC66", "#FF3366")
ggplot(data = dat.res, aes(x = datty, y = residual, 
                           fill = datty)) + 
  geom_bar(stat = 'identity', color = 'black', linewidth = 1, 
           position = position_dodge(width = 0.5)) + 
  geom_text(aes(label = round(residual, 3)), vjust = -0.5, size = 6) + 
  scale_fill_manual(values = set_col) + 
  scale_y_break(c(0.05, round(min(vals) - 0.0001, 4)), 
                ticklabels = seq(round(min(vals) - 0.001, 2), round(max(vals) + 0.001, 2), 0.01),
                scales = 10) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, round(max(vals) + 0.002, 3)),
                     breaks = c(0, seq(round(min(vals) - 0.001, 2), round(max(vals) + 0.001, 2), 0.01))) + 
  theme_bw() + 
  xlab(NULL) + 
  ylab('Residual Ratio') + 
  theme(legend.position = 'none', 
        axis.title = element_text(size = 20), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y.left = element_text(size = 18), 
        axis.text.x = element_text(size = 18), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


## e. (MTBLS290) Intra-group Clustering ####
library(igraph)
library(ggraph)
library(tidygraph)

if (dat.title == 'teadata_DMSO') {
  sig  <- c(1, 0)
} else {
  sig <- c(0.5, 0.5)
}

majgrp.f <- as.factor(majGroup)
majgrp.lev <- levels(majgrp.f)
majgrp.num <- nlevels(majgrp.f)

dat.num <- length(dat.list)
dat.name <- names(dat.list)

Vmeasure <- list()
Vmeasure.sd <- list()
for (i in c(1: dat.num)) {
  
  vm.i <- matrix(0, majgrp.num, 3)
  vm.sd.i <- matrix(0, majgrp.num, 3)
  for (g in c(1: majgrp.num)) {
    mingrp.g <- minGroup[majGroup == majgrp.lev[g]]
    mingrp.g.f <- as.factor(mingrp.g)
    mingrp.g.lev <- levels(mingrp.g.f)
    mingrp.g.num <- nlevels(mingrp.g.f)
    
    if (mingrp.g.num > 1) {
      
      # clustering
      rep.time <- 100
      vms <- matrix(0, rep.time, 3)
      for (t in c(1: rep.time)) {
        clusId.i.g <- ImproKmCluster(dat.list[[i]][majGroup == majgrp.lev[g], ], 
                                     k = mingrp.g.num, varyk = 1, optimal = T, 
                                     sig = sig)
        vms[t, ] <- unlist(getVmeasure(clusId.i.g, mingrp.g))
      }

      vm.i[g, ] <- colMeans(vms)
      vm.sd.i[g, ] <- apply(vms, MARGIN = 2, FUN = sd)
      
    } else {
      vm.i[g, ] <- NA
      vm.sd.i[g, ] <- NA
    }
    
  }
  
  Vmeasure[[i]] <- vm.i
  Vmeasure.sd[[i]] <- vm.sd.i
}


# Eval 2 - Preservation of Biological Effect ####

## a. (Simulated) Effects of Differential Features ####
### I. Individual Differences ####
rep.time <- 100
n.ratio <- seq(0.05, 0.45, 0.1)

DEmetFval.LSCN <- matrix(0, rep.time, length(n.ratio))
DEmetFval.CSN <- matrix(0, rep.time, length(n.ratio))
DEmetFval.L2N <- matrix(0, rep.time, length(n.ratio))
DEmetFval.PQN <- matrix(0, rep.time, length(n.ratio))
DEmetFval.QT <- matrix(0, rep.time, length(n.ratio))

DEmetVIP.LSCN <- matrix(0, rep.time, length(n.ratio))
DEmetVIP.CSN <- matrix(0, rep.time, length(n.ratio))
DEmetVIP.L2N <- matrix(0, rep.time, length(n.ratio))
DEmetVIP.PQN <- matrix(0, rep.time, length(n.ratio))
DEmetVIP.QT <- matrix(0, rep.time, length(n.ratio))

for (i in c(1: length(n.ratio))) {
  
  mfval.LSCN <- rep(0, rep.time)
  mfval.CSN <- rep(0, rep.time)
  mfval.L2N <- rep(0, rep.time)
  mfval.PQN <- rep(0, rep.time)
  mfval.QT <- rep(0, rep.time)

  mVIP.LSCN <- rep(0, rep.time)
  mVIP.CSN <- rep(0, rep.time)
  mVIP.L2N <- rep(0, rep.time)
  mVIP.PQN <- rep(0, rep.time)
  mVIP.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, grpEff.lamda = c(0.3, 1), 
                               diltF.m = 30, diltF.sd = 10, 
                               indiv.sd.ratio = n.ratio[i])
    dat_simu <- dat_simu.res$SimuDat
    grp_simu <- dat_simu.res$Para$group
    dat.GT <- dat_simu.res$InitDat
    DEmet.GT <- dat_simu.res$Para$DEmetIdx
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    dat.nor.CSN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.CSN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.CSN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.CSN[k] <- mean(DE.VIP, na.rm = T)
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    dat.nor.L2N <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.L2N, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.L2N[k] <- mean(DE.Fval, na.rm = T)
    mVIP.L2N[k] <- mean(DE.VIP, na.rm = T)
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    dat.nor.PQN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.PQN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.PQN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.PQN[k] <- mean(DE.VIP, na.rm = T)
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.QT, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.QT[k] <- mean(DE.Fval, na.rm = T)
    mVIP.QT[k] <- mean(DE.VIP, na.rm = T)
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, 
                    max_ite = 30, print.txt = F)
    dat.nor.LSCN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.LSCN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.LSCN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.LSCN[k] <- mean(DE.VIP, na.rm = T)
    
  }
  
  DEmetFval.LSCN[, i] <- mfval.LSCN
  DEmetFval.CSN[, i] <- mfval.CSN
  DEmetFval.L2N[, i] <- mfval.L2N
  DEmetFval.PQN[, i] <- mfval.PQN
  DEmetFval.QT[, i] <- mfval.QT
  
  DEmetVIP.LSCN[, i] <- mVIP.LSCN
  DEmetVIP.CSN[, i] <- mVIP.CSN
  DEmetVIP.L2N[, i] <- mVIP.L2N
  DEmetVIP.PQN[, i] <- mVIP.PQN
  DEmetVIP.QT[, i] <- mVIP.QT
  
  cat(i, '\n')
}

# F-value
vals <- c(as.vector(DEmetFval.LSCN), as.vector(DEmetFval.CSN), 
          as.vector(DEmetFval.L2N), as.vector(DEmetFval.PQN), 
          as.vector(DEmetFval.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(n.ratio)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

nr <- rep(rep(n.ratio, each = rep.time), 5)

dat.DEFval <- data.frame(dattype = dat.type, 
                         NR = nr, 
                         Fval = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.DEFval, aes(x = NR, y = Fval,
                                     color = dattype, fill = dattype, 
                                     linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = grpVar.ratio) + 
  theme_bw() + 
  xlab('Individual Differences (r)') + 
  ylab('Average F-value of Potential DE Features') + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

# VIP scores
vals <- c(as.vector(DEmetVIP.LSCN), as.vector(DEmetVIP.CSN), 
          as.vector(DEmetVIP.L2N), as.vector(DEmetVIP.PQN), 
          as.vector(DEmetVIP.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(n.ratio)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

nr <- rep(rep(n.ratio, each = rep.time), 5)

dat.DEVIP <- data.frame(dattype = dat.type, 
                        NR = nr, 
                        VIP = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.DEVIP, aes(x = NR, y = VIP,
                                    color = dattype, fill = dattype, 
                                    linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = grpVar.ratio) + 
  theme_bw() + 
  xlab('Individual Differences (r)') + 
  ylab('Average VIP of Potential DE Features') + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


### II. Inter-group Heterogeneity ####
rep.time <- 100
grpVar.ratio <- seq(0.05, 0.45, 0.1)

DEmetFval.LSCN <- matrix(0, rep.time, length(grpVar.ratio))
DEmetFval.CSN <- matrix(0, rep.time, length(grpVar.ratio))
DEmetFval.L2N <- matrix(0, rep.time, length(grpVar.ratio))
DEmetFval.PQN <- matrix(0, rep.time, length(grpVar.ratio))
DEmetFval.QT <- matrix(0, rep.time, length(grpVar.ratio))

DEmetVIP.LSCN <- matrix(0, rep.time, length(grpVar.ratio))
DEmetVIP.CSN <- matrix(0, rep.time, length(grpVar.ratio))
DEmetVIP.L2N <- matrix(0, rep.time, length(grpVar.ratio))
DEmetVIP.PQN <- matrix(0, rep.time, length(grpVar.ratio))
DEmetVIP.QT <- matrix(0, rep.time, length(grpVar.ratio))

for (i in c(1: length(grpVar.ratio))) {
  
  mfval.LSCN <- rep(0, rep.time)
  mfval.CSN <- rep(0, rep.time)
  mfval.L2N <- rep(0, rep.time)
  mfval.PQN <- rep(0, rep.time)
  mfval.QT <- rep(0, rep.time)

  mVIP.LSCN <- rep(0, rep.time)
  mVIP.CSN <- rep(0, rep.time)
  mVIP.L2N <- rep(0, rep.time)
  mVIP.PQN <- rep(0, rep.time)
  mVIP.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, grpEff.lamda = c(0.3, 1), 
                               diltF.m = 30, grpVar.ratio = grpVar.ratio[i])
    dat_simu <- dat_simu.res$SimuDat
    grp_simu <- dat_simu.res$Para$group
    DEmet.GT <- dat_simu.res$Para$DEmetIdx
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    dat.nor.CSN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.CSN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.CSN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.CSN[k] <- mean(DE.VIP, na.rm = T)
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    dat.nor.L2N <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.L2N, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.L2N[k] <- mean(DE.Fval, na.rm = T)
    mVIP.L2N[k] <- mean(DE.VIP, na.rm = T)
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    dat.nor.PQN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.PQN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.PQN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.PQN[k] <- mean(DE.VIP, na.rm = T)
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.QT, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.QT[k] <- mean(DE.Fval, na.rm = T)
    mVIP.QT[k] <- mean(DE.VIP, na.rm = T)
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, 
                    max_ite = 30, print.txt = F)
    dat.nor.LSCN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.LSCN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.LSCN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.LSCN[k] <- mean(DE.VIP, na.rm = T)
    
  }
  
  DEmetFval.LSCN[, i] <- mfval.LSCN
  DEmetFval.CSN[, i] <- mfval.CSN
  DEmetFval.L2N[, i] <- mfval.L2N
  DEmetFval.PQN[, i] <- mfval.PQN
  DEmetFval.QT[, i] <- mfval.QT

  DEmetVIP.LSCN[, i] <- mVIP.LSCN
  DEmetVIP.CSN[, i] <- mVIP.CSN
  DEmetVIP.L2N[, i] <- mVIP.L2N
  DEmetVIP.PQN[, i] <- mVIP.PQN
  DEmetVIP.QT[, i] <- mVIP.QT
  
  cat(i, '\n')
}

# F-value
vals <- c(as.vector(DEmetFval.LSCN), as.vector(DEmetFval.CSN), 
          as.vector(DEmetFval.L2N), as.vector(DEmetFval.PQN), 
          as.vector(DEmetFval.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(grpVar.ratio)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

var.ratio <- rep(rep(grpVar.ratio, each = rep.time), 5)

dat.DEFval <- data.frame(dattype = dat.type, 
                         grpVar.ratio = var.ratio, 
                         Fval = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.DEFval, aes(x = grpVar.ratio, y = Fval,
                                     color = dattype, fill = dattype, 
                                     linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = grpVar.ratio) + 
  theme_bw() + 
  xlab('Inter-group Heterogeneity (k)') + 
  ylab('Average F-value of Potential DE Features') + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

# VIP scores
vals <- c(as.vector(DEmetVIP.LSCN), as.vector(DEmetVIP.CSN), 
          as.vector(DEmetVIP.L2N), as.vector(DEmetVIP.PQN), 
          as.vector(DEmetVIP.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(grpVar.ratio)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

var.ratio <- rep(rep(grpVar.ratio, each = rep.time), 5)

dat.DEVIP <- data.frame(dattype = dat.type, 
                        grpVar.ratio = var.ratio, 
                        VIP = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.DEVIP, aes(x = grpVar.ratio, y = VIP,
                                    color = dattype, fill = dattype, 
                                    linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = grpVar.ratio) + 
  theme_bw() + 
  xlab('Inter-group Heterogeneity') + 
  ylab('Average VIP of Potential DE Features') + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


### III. Intra-group Heterogeneity ####
rep.time <- 100
sg.num <- seq(2, 6, 1)

DEmetFval.LSCN <- matrix(0, rep.time, length(sg.num))
DEmetFval.CSN <- matrix(0, rep.time, length(sg.num))
DEmetFval.L2N <- matrix(0, rep.time, length(sg.num))
DEmetFval.PQN <- matrix(0, rep.time, length(sg.num))
DEmetFval.QT <- matrix(0, rep.time, length(sg.num))

DEmetVIP.LSCN <- matrix(0, rep.time, length(sg.num))
DEmetVIP.CSN <- matrix(0, rep.time, length(sg.num))
DEmetVIP.L2N <- matrix(0, rep.time, length(sg.num))
DEmetVIP.PQN <- matrix(0, rep.time, length(sg.num))
DEmetVIP.QT <- matrix(0, rep.time, length(sg.num))

for (i in c(1: length(sg.num))) {

  mfval.LSCN <- rep(0, rep.time)
  mfval.CSN <- rep(0, rep.time)
  mfval.L2N <- rep(0, rep.time)
  mfval.PQN <- rep(0, rep.time)
  mfval.QT <- rep(0, rep.time)

  mVIP.LSCN <- rep(0, rep.time)
  mVIP.CSN <- rep(0, rep.time)
  mVIP.L2N <- rep(0, rep.time)
  mVIP.PQN <- rep(0, rep.time)
  mVIP.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, grpEff.lamda = c(0.3, 1), 
                               diltF.m = 30, subgrpNum = sg.num[i])
    dat_simu <- dat_simu.res$SimuDat
    grp_simu <- dat_simu.res$Para$group
    DEmet.GT <- dat_simu.res$Para$DEmetIdx
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    dat.nor.CSN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.CSN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.CSN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.CSN[k] <- mean(DE.VIP, na.rm = T)
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    dat.nor.L2N <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.L2N, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.L2N[k] <- mean(DE.Fval, na.rm = T)
    mVIP.L2N[k] <- mean(DE.VIP, na.rm = T)
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    dat.nor.PQN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.PQN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.PQN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.PQN[k] <- mean(DE.VIP, na.rm = T)
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.QT, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.QT[k] <- mean(DE.Fval, na.rm = T)
    mVIP.QT[k] <- mean(DE.VIP, na.rm = T)
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, 
                    max_ite = 30, print.txt = F)
    dat.nor.LSCN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.LSCN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.LSCN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.LSCN[k] <- mean(DE.VIP, na.rm = T)
    
  }

  DEmetFval.LSCN[, i] <- mfval.LSCN
  DEmetFval.CSN[, i] <- mfval.CSN
  DEmetFval.L2N[, i] <- mfval.L2N
  DEmetFval.PQN[, i] <- mfval.PQN
  DEmetFval.QT[, i] <- mfval.QT

  DEmetVIP.LSCN[, i] <- mVIP.LSCN
  DEmetVIP.CSN[, i] <- mVIP.CSN
  DEmetVIP.L2N[, i] <- mVIP.L2N
  DEmetVIP.PQN[, i] <- mVIP.PQN
  DEmetVIP.QT[, i] <- mVIP.QT
  
  cat(i, '\n')
}

# F-value
vals <- c(as.vector(DEmetFval.LSCN), as.vector(DEmetFval.CSN), 
          as.vector(DEmetFval.L2N), as.vector(DEmetFval.PQN), 
          as.vector(DEmetFval.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(sg.num)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

sgn <- rep(rep(sg.num, each = rep.time), 5)

dat.DEFval <- data.frame(dattype = dat.type, 
                         SGN = sgn, 
                         Fval = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.DEFval, aes(x = SGN, y = Fval,
                                     color = dattype, fill = dattype, 
                                     linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = sg.num) + 
  theme_bw() + 
  xlab('Intra-group Heterogeneity (n)') + 
  ylab('Average F-value of Potential DE Features') + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

# VIP scores
vals <- c(as.vector(DEmetVIP.LSCN), as.vector(DEmetVIP.CSN), 
          as.vector(DEmetVIP.L2N), as.vector(DEmetVIP.PQN), 
          as.vector(DEmetVIP.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(sg.num)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

sgn <- rep(rep(sg.num, each = rep.time), 5)

dat.DEVIP <- data.frame(dattype = dat.type, 
                        SGN = sgn, 
                        VIP = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.DEVIP, aes(x = SGN, y = VIP,
                                    color = dattype, fill = dattype, 
                                    linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = sg.num) + 
  theme_bw() + 
  xlab('Intra-group Heterogeneity') + 
  ylab('Average VIP of Potential DE Features') + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


### IV. Dilution Effects ####
rep.time <- 100
sv.m <- seq(10, 90, 20)

DEmetFval.LSCN <- matrix(0, rep.time, length(sv.m))
DEmetFval.CSN <- matrix(0, rep.time, length(sv.m))
DEmetFval.L2N <- matrix(0, rep.time, length(sv.m))
DEmetFval.PQN <- matrix(0, rep.time, length(sv.m))
DEmetFval.QT <- matrix(0, rep.time, length(sv.m))

DEmetVIP.LSCN <- matrix(0, rep.time, length(sv.m))
DEmetVIP.CSN <- matrix(0, rep.time, length(sv.m))
DEmetVIP.L2N <- matrix(0, rep.time, length(sv.m))
DEmetVIP.PQN <- matrix(0, rep.time, length(sv.m))
DEmetVIP.QT <- matrix(0, rep.time, length(sv.m))

for (i in c(1: length(sv.m))) {

  mfval.LSCN <- rep(0, rep.time)
  mfval.CSN <- rep(0, rep.time)
  mfval.L2N <- rep(0, rep.time)
  mfval.PQN <- rep(0, rep.time)
  mfval.QT <- rep(0, rep.time)

  mVIP.LSCN <- rep(0, rep.time)
  mVIP.CSN <- rep(0, rep.time)
  mVIP.L2N <- rep(0, rep.time)
  mVIP.PQN <- rep(0, rep.time)
  mVIP.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, grpEff.lamda = c(0.3, 1), 
                               diltF.m = sv.m[i])
    dat_simu <- dat_simu.res$SimuDat
    grp_simu <- dat_simu.res$Para$group
    DEmet.GT <- dat_simu.res$Para$DEmetIdx
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    dat.nor.CSN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.CSN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.CSN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.CSN[k] <- mean(DE.VIP, na.rm = T)
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    dat.nor.L2N <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.L2N, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.L2N[k] <- mean(DE.Fval, na.rm = T)
    mVIP.L2N[k] <- mean(DE.VIP, na.rm = T)
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    dat.nor.PQN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.PQN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.PQN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.PQN[k] <- mean(DE.VIP, na.rm = T)
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.QT, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.QT[k] <- mean(DE.Fval, na.rm = T)
    mVIP.QT[k] <- mean(DE.VIP, na.rm = T)
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, 
                    max_ite = 30, print.txt = F)
    dat.nor.LSCN <- dat.nor$dat.norm
    DEana.res <- getDEmet(dat.nor.LSCN, group = grp_simu, adjustP = T)
    DE.Fval <- DEana.res$Fval[DEmet.GT]
    DE.VIP <- DEana.res$VIP[DEmet.GT]
    mfval.LSCN[k] <- mean(DE.Fval, na.rm = T)
    mVIP.LSCN[k] <- mean(DE.VIP, na.rm = T)
  }

  DEmetFval.LSCN[, i] <- mfval.LSCN
  DEmetFval.CSN[, i] <- mfval.CSN
  DEmetFval.L2N[, i] <- mfval.L2N
  DEmetFval.PQN[, i] <- mfval.PQN
  DEmetFval.QT[, i] <- mfval.QT

  DEmetVIP.LSCN[, i] <- mVIP.LSCN
  DEmetVIP.CSN[, i] <- mVIP.CSN
  DEmetVIP.L2N[, i] <- mVIP.L2N
  DEmetVIP.PQN[, i] <- mVIP.PQN
  DEmetVIP.QT[, i] <- mVIP.QT
  
  cat(i, '\n')
}

# F-value
vals <- c(as.vector(DEmetFval.LSCN), as.vector(DEmetFval.CSN), 
          as.vector(DEmetFval.L2N), as.vector(DEmetFval.PQN), 
          as.vector(DEmetFval.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(sv.m)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

svms <- rep(rep(sv.m, each = rep.time), 5)

dat.DEFval <- data.frame(dattype = dat.type, 
                         SVM = svms, 
                         Fval = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.DEFval, aes(x = SVM, y = Fval,
                                     color = dattype, fill = dattype, 
                                     linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = grpVar.ratio) + 
  theme_bw() + 
  xlab('Dilution Effects (m)') + 
  ylab('Average F-value of Potential DE Features') + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

# VIP score
vals <- c(as.vector(DEmetVIP.LSCN), as.vector(DEmetVIP.CSN), 
          as.vector(DEmetVIP.L2N), as.vector(DEmetVIP.PQN), 
          as.vector(DEmetVIP.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(sv.m)),
                   levels = c('CSN', 'L2-norm', 'PQN', 
                              'Quantile', 'LSCN'))

svms <- rep(rep(sv.m, each = rep.time), 5)

dat.DEVIP <- data.frame(dattype = dat.type, 
                        SVM = svms, 
                        VIP = vals)

set_col <- c('magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 4), 1)


ggplot() + 
  geom_smooth(data = dat.DEVIP, aes(x = SVM, y = VIP,
                                    color = dattype, fill = dattype, 
                                    linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = grpVar.ratio) + 
  theme_bw() + 
  xlab('Dilution Effect (m)') + 
  ylab('Average VIP of Potential DE Features') + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


## b. (Real) Identification of differential features and correlations ####
conf_lev <- 0.05

# Unnorm
DEana.res <- getDEmet(dat.list$Unnorm, group = majGroup, 
                       p.level = conf_lev, adjustP = T)
DEmetIdx.bef <- DEana.res$DEmetIdx
Pval.bef <- DEana.res$Pval
nlogP.bef <- -log10(Pval.bef)
VIP.bef <- DEana.res$VIP

# LSCN
DEana.res <- getDEmet(dat.list$LSCN, group = majGroup, 
                       p.level = conf_lev, adjustP = T)
DEmetIdx.LSCN <- DEana.res$DEmetIdx
Pval.LSCN <- DEana.res$Pval
nlogP.LSCN <- -log10(Pval.LSCN)
VIP.LSCN <- DEana.res$VIP

# CSN
DEana.res <- getDEmet(dat.list$CSN, group = majGroup, 
                       p.level = conf_lev, adjustP = T)
DEmetIdx.CSN <- DEana.res$DEmetIdx
Pval.CSN <- DEana.res$Pval
nlogP.CSN <- -log10(Pval.CSN)
VIP.CSN <- DEana.res$VIP

# L2N
DEana.res <- getDEmet(dat.list$L2N, group = majGroup, 
                       p.level = conf_lev, adjustP = T)
DEmetIdx.L2N <- DEana.res$DEmetIdx
Pval.L2N <- DEana.res$Pval
nlogP.L2N <- -log10(Pval.L2N)
VIP.L2N <- DEana.res$VIP

# PQN
DEana.res <- getDEmet(dat.list$PQN, group = majGroup, 
                       p.level = conf_lev, adjustP = T)
DEmetIdx.PQN <- DEana.res$DEmetIdx
Pval.PQN <- DEana.res$Pval
nlogP.PQN <- -log10(Pval.PQN)
VIP.PQN <- DEana.res$VIP

# QT
DEana.res <- getDEmet(dat.list$QT, group = majGroup, 
                       p.level = conf_lev, adjustP = T)
DEmetIdx.QT <- DEana.res$DEmetIdx
Pval.QT <- DEana.res$Pval
nlogP.QT <- -log10(Pval.QT)
VIP.QT <- DEana.res$VIP


DE_Pval.list <- list(CSN = Pval.CSN[DEmetIdx.CSN], 
                     L2N = Pval.L2N[DEmetIdx.L2N], 
                     PQN = Pval.PQN[DEmetIdx.PQN], 
                     QT = Pval.QT[DEmetIdx.QT], 
                     LSCN = Pval.LSCN[DEmetIdx.LSCN])

### I. Number of differential features ####
DEmetNum.df <- data.frame(CSN = length(DEmetIdx.CSN), 
                          L2N = length(DEmetIdx.L2N), 
                          PQN = length(DEmetIdx.PQN), 
                          QT = length(DEmetIdx.QT), 
                          LSCN = length(DEmetIdx.LSCN))

write.table(DEmetNum.df, here("Output", paste0(dat.title, '_DEmetNum.txt')), 
            row.names = F)

### II. P-values and VIPs of DE features ####
DEmetIdx.list <- list(CSN = DEmetIdx.CSN, 
                      L2N = DEmetIdx.L2N, 
                      PQN = DEmetIdx.PQN, 
                      QT = DEmetIdx.QT, 
                      LSCN = DEmetIdx.LSCN)

dat.num <- length(DEmetIdx.list)
commDE.mat <- matrix(0, dat.num, p)
for (i in c(1: dat.num)) {
  commDE.mat[i, DEmetIdx.list[[i]]] <- 1
}

commDE.freq <- colSums(commDE.mat)
commDE.freq.ordIdx <- order(commDE.freq, decreasing = T)
commDE.freq.sort <- commDE.freq[commDE.freq.ordIdx]
commDE.freq.ordIdx <- commDE.freq.ordIdx[commDE.freq.sort > 0]
commDE.freq.sort <- commDE.freq.sort[commDE.freq.sort > 0]
allcommDE.Idx <- commDE.freq.ordIdx[commDE.freq.sort == dat.num]
allcommDE.nlogP.LSCN <- nlogP.LSCN[allcommDE.Idx]
allcommDE.Idx.sort <- allcommDE.Idx[order(allcommDE.nlogP.LSCN, decreasing = T)]
commDE.freq.ordIdx[commDE.freq.sort == dat.num] <- allcommDE.Idx.sort

plot_name <- MetName
if (dat.title == 'ST003484') {
  plot_name[commDE.freq.ordIdx] <- c('CA', 'IP', 'Lal', 'aCD', 'Spe', 
                                     'Lly', 'Las', 'Lpr', 'IMP', 
                                     'I3a', 'Nic', 'dTMP')
} else if (dat.title == 'MTBLS290') {
  fullname <- read.csv(here("Data", "metname_II.csv"), header = T)
  DEmet_fullname <- fullname$metabolite_name[commDE.freq.ordIdx]
  write.csv(DEmet_fullname, here("Output", 'DEmet_fullname_II.csv'))
} else if (dat.title == 'ST002178') {
  name1 <- plot_name[commDE.freq.ordIdx]
  cleaned_name <- gsub("[^A-Za-z]", "", name1)
  result_name <- substr(cleaned_name, 1, 3)
  result_name[which(duplicated(result_name))] <- c('Ami1', 'Ami2', 'His1', 'Asp1', 
                                                   'Ami3', 'CGl1', 'Met1', 'Glu1', 
                                                   'Tri1', 'CGl2', 'CGl3')
  plot_name[commDE.freq.ordIdx] <- result_name
} 

commDEmet.name <- plot_name[commDE.freq.ordIdx]

met.name <- factor(c(plot_name[DEmetIdx.LSCN], 
                     plot_name[DEmetIdx.CSN], 
                     plot_name[DEmetIdx.L2N], 
                     plot_name[DEmetIdx.PQN], 
                     plot_name[DEmetIdx.QT]), 
                   levels = commDEmet.name)

commDE.nlogP <- c(nlogP.LSCN[DEmetIdx.LSCN], 
                  nlogP.CSN[DEmetIdx.CSN], 
                  nlogP.L2N[DEmetIdx.L2N], 
                  nlogP.PQN[DEmetIdx.PQN], 
                  nlogP.QT[DEmetIdx.QT])

commDE.VIP <- c(VIP.LSCN[DEmetIdx.LSCN], 
                VIP.CSN[DEmetIdx.CSN], 
                VIP.L2N[DEmetIdx.L2N], 
                VIP.PQN[DEmetIdx.PQN], 
                VIP.QT[DEmetIdx.QT])

dat.type <- factor(c(rep('LSCN', length(DEmetIdx.LSCN)), 
                     rep('CSN', length(DEmetIdx.CSN)), 
                     rep('L2-norm', length(DEmetIdx.L2N)), 
                     rep('PQN', length(DEmetIdx.PQN)), 
                     rep('Quantile', length(DEmetIdx.QT))),
                   levels = c('CSN', 'L2-norm', 
                              'PQN', 'Quantile', 'LSCN'))

dat.commDE <- data.frame(datty = dat.type,
                         DEmet = met.name,
                         nlogP = commDE.nlogP,
                         VIP = commDE.VIP)

if (length(commDEmet.name) < 30) {
  metname.size <- 12
} else if (length(commDEmet.name) >= 30 & length(commDEmet.name) < 60) {
  metname.size <- 9
} else if (length(commDEmet.name) >= 60 & length(commDEmet.name) < 90) {
  metname.size <- 6
}

if (dat.title == 'ST003484') {
  epd1 <- 0.5
  epd2 <- 0.5
  alp <- 1
}else {
  epd1 <- 2
  epd2 <- 3
  alp <- 0.7
}

# plot
p1 <- ggplot(data = dat.commDE, aes(x = DEmet, y = datty)) +
  geom_point(aes(fill = nlogP, size = VIP), 
             shape = 21, color = 'black', alpha = alp) + 
  scale_fill_gradientn(limits = c(0, max(dat.commDE$nlogP)), 
                       values = seq(0, 1, 0.1),
                       colors = c('#6699cc', '#FFFF99', '#CC3333')) + 
  scale_x_discrete(limits = commDEmet.name, expand = expansion(add = c(epd1, epd2))) + 
  theme_bw() + 
  xlab(NULL) + ylab(NULL)

if (dat.title == 'ST003484') {
  
  p1 <- p1 + theme(panel.grid = element_blank(),
                   axis.text.y = element_text(size = 20, face = 'bold'), 
                   axis.text.x = element_text(size = 15, face = 'bold'), 
                   legend.text = element_text(family = 'sans', size = 15), 
                   legend.title = element_text(family = 'sans', size = 18), 
                   panel.border = element_blank(), 
                   axis.line.x.bottom = element_line(linewidth = 1), 
                   axis.line.y.left = element_line(linewidth = 1), 
                   plot.margin = margin(t = 30, unit = 'pt'))
  
} else if (dat.title == 'MTBLS290' | dat.title == 'teadata_DMSO') {
  p1 <- p1 + theme(panel.grid = element_blank(),
                   axis.text.y = element_text(size = 20, face = 'bold'), 
                   axis.text.x = element_blank(), 
                   axis.ticks.x = element_blank(), 
                   legend.text = element_text(family = 'sans', size = 15), 
                   legend.title = element_text(family = 'sans', size = 18), 
                   panel.border = element_blank(), 
                   axis.line.x.bottom = element_line(linewidth = 1), 
                   axis.line.y.left = element_line(linewidth = 1), 
                   plot.margin = margin(t = 30, unit = 'pt'))
} else if (dat.title == 'ST002178') {
  p1 <- p1 + theme(panel.grid = element_blank(),
                   axis.text.y = element_text(size = 15, face = 'bold'), 
                   axis.text.x = element_text(size = metname.size, face = 'bold', 
                                              angle = 30, hjust = 0.5, vjust = 0.5, 
                                              margin = margin(t = -3)), 
                   legend.text = element_text(family = 'sans', size = 10), 
                   legend.title = element_text(family = 'sans', size = 12), 
                   panel.border = element_blank(), 
                   axis.line.x.bottom = element_line(linewidth = 1), 
                   axis.line.y.left = element_line(linewidth = 1), 
                   plot.margin = margin(t = 30, unit = 'pt'))
}

print(p1)


### III. Differential PCC network ####
library(ggraph)
library(igraph)

majgrp.f <- as.factor(majGroup)
majgrp.lev <- levels(majgrp.f)
majgrp.num <- nlevels(majgrp.f)

dat.num <- length(dat.list)
dat.name <- names(dat.list)
tot.pcc.num <- p * (p - 1) / 2

if (dat.title == 'MTBLS290') {
  MetName <- read.csv(here("Data", "metname_II.csv"), header = T)
  MetName <- MetName$metabolite_name
  kegg_id <- rep(NA, p)
} else if (dat.title == 'teadata_DMSO') {
  kegg_id <- rep(NA, p)
}

conf_lev <- 0.05
topk <- round(0.01 * tot.pcc.num)
mink <- 10


set_col1 <- c("#0066FF", "#FF3366")
set_col2 <- c("#FF66FF", "#FFCC00", "#0066FF", "#00CC66", "#FF3366")

diffPCC.Amat <- list()
for (i1 in c(1: (majgrp.num - 1))) {
  for (i2 in c((i1 + 1): majgrp.num)) {
    
    idendiffPCC.mat <- matrix(0, dat.num - 1, tot.pcc.num)
    diffP.conncompNum <- rep(0, dat.num - 1)
    
    cap.txt <- paste(majgrp.lev[i1], 'vs', majgrp.lev[i2])
    
    A.gi_12 <- list()
    clust.measure <- matrix(0, dat.num - 1, 3)
    for (i in c(1: (dat.num-1))) {
      dat.i <- dat.list[[i + 1]]
      name.i <- dat.name[i + 1]
      
      dat.i.gi1 <- dat.i[majGroup == majgrp.lev[i1], ]
      dat.i.gi2 <- dat.i[majGroup == majgrp.lev[i2], ]
      
      diffPCC.res <- getdiffPCC(dat.i.gi1, dat.i.gi2, 
                                cor.method = 'spearman', 
                                conf_alpha = conf_lev)
      P.i.gi_12 <- diffPCC.res$P.Mat
      A.i.gi_12 <- diffPCC.res$A
      edge.num <- sum(A.i.gi_12[upper.tri(A.i.gi_12)]) 
      if (edge.num >= topk | edge.num < mink) {
        all.pval <- P.i.gi_12[upper.tri(P.i.gi_12)]
        pval.thresh <- sort(all.pval)[topk]
        A.i.gi_12 <- 1 * (P.i.gi_12 <= pval.thresh)
      } 
      A.gi_12[[i]] <- A.i.gi_12
      idendiffPCC.mat[i, ] <- A.i.gi_12[upper.tri(A.i.gi_12)]
      
      # plot
      deg <- colSums(A.i.gi_12)
      DEPCC_metIdx.i <- which(deg != 0)
      A.i.gi_12.sub <- A.i.gi_12[DEPCC_metIdx.i, DEPCC_metIdx.i]
      graph.diffP <- graph_from_adjacency_matrix(A.i.gi_12.sub, mode = 'undirected', 
                                                 weighted = NULL, diag = FALSE)
      
      is_DEmet <- rep(0, p)
      is_DEmet[DEmetIdx.list[[i]]] <- 1
      is_DEmet <- as.factor(is_DEmet)
      is_DEmet.sub <- is_DEmet[DEPCC_metIdx.i]
      deg.sub <- deg[DEPCC_metIdx.i]
      
      # 1. Generate a layout for the entire network (including all sub-networks)
      if (p < 100) {
        global_layout <- layout_with_kk(graph.diffP)  # kk
      } else {
        global_layout <- layout_with_fr(graph.diffP)  # fr
      }
      
      # 2. Find the node index of the largest connected subgraph
      comp <- components(graph.diffP)
      max_comp_id <- which(comp$csize == max(comp$csize))
      max_nodes <- which(comp$membership == max_comp_id)
      
      large_comp_id <- which(comp$csize == max(comp$csize) | comp$csize / sum(comp$csize) > 0.3)
      if (length(large_comp_id) == 1) {
        lc_nodes <- which(comp$membership == large_comp_id)
      } else {
        lc_nodes <- c()
        for (ii in c(1: length(large_comp_id))) {
          lc_nodes <- append(lc_nodes, which(comp$membership == large_comp_id[ii]))
        }
      }
      
      # The names of all metabolites in the differential correlation network
      DEPCC_metName <- MetName[DEPCC_metIdx.i]
      write.csv(DEPCC_metName, here("Output", paste0('DPmN_', name.i, '_', 
                                                     datIdx, '_', majgrp.lev[i1], 
                                                     '_vs_', majgrp.lev[i2], '.csv')))
      # The metabolite keggID in the maximum connected subgraph
      DEPCC_LCmetIdx.i <- DEPCC_metIdx.i[lc_nodes]
      DPLCmet_name <- MetName[DEPCC_LCmetIdx.i]
      DPLCmet_keggID <- kegg_id[DEPCC_LCmetIdx.i]
      DPLCmet <- data.frame(name = DPLCmet_name, 
                            keggID = DPLCmet_keggID)
      write.csv(DPLCmet, here("Output", paste0('DPLCmet_', name.i, '_', 
                                               datIdx, '_', majgrp.lev[i1], 
                                               '_vs_', majgrp.lev[i2], '.csv')))
      
      # 3. The metabolite keggID in the maximum connected subgraph
      center_x <- mean(global_layout[max_nodes, 1])
      center_y <- mean(global_layout[max_nodes, 2])
      
      # 4. Translate the other sub-networks (move closer to the center)
      adjusted_layout <- global_layout
      for (ii in 1:length(comp$csize)) {
        if (ii != max_comp_id) {
          # Obtain the node index of the current sub-network
          current_nodes <- which(comp$membership == ii)
          
          # Obtain the node index of the current sub-network
          current_center_x <- mean(global_layout[current_nodes, 1])
          current_center_y <- mean(global_layout[current_nodes, 2])
          
          # Calculate the translation vector (moving towards the maximum subnetwork)
          delta_x <- (center_x - current_center_x) * 0.4  # 0.7是靠拢强度
          delta_y <- (center_y - current_center_y) * 0.4
          
          # Applying
          adjusted_layout[current_nodes, 1] <- adjusted_layout[current_nodes, 1] + delta_x
          adjusted_layout[current_nodes, 2] <- adjusted_layout[current_nodes, 2] + delta_y
        }
      }
      
      p1 <- ggraph(graph = graph.diffP, layout = adjusted_layout) + 
        geom_edge_link() + 
        geom_node_point(aes(size = deg.sub, fill = is_DEmet.sub), 
                        color = "white", shape = 21, stroke = 1) + 
        scale_fill_manual(values = set_col1) + 
        scale_size_continuous(range = c(3, 7)) + 
        theme_graph() + 
        theme(legend.position = 'none', 
              plot.title = element_text(size = 30, family = 'sans', 
                                        face = 'bold', hjust = 0.5), 
              plot.caption = element_text(size = 25, family = 'sans', 
                                          hjust = 0.5)) + 
        labs(title = name.i, caption = cap.txt)
      
      print(p1)
      
      comp <- components(graph.diffP)
      diffP.conncompNum[i] <- sum(comp$csize > 1)
    }
    
    diffPCC.Amat[[cap.txt]] <- A.gi_12
    
    # Number of connected components
    diffP.conncompNum.df <- as.data.frame(t(diffP.conncompNum))
    colnames(diffP.conncompNum.df) <- dat.name[-1]
    write.table(diffP.conncompNum.df, here("Output", paste0(dat.title, '_', majgrp.lev[i1], 
                                                            '_vs_', majgrp.lev[i2], '_connCompNum.txt')), 
                row.names = F)
    
    # Number of differential PCCs
    DEPCCNum.df <- as.data.frame(t(rowSums(idendiffPCC.mat)))
    colnames(DEPCCNum.df) <- dat.name[-1]
    write.table(DEPCCNum.df, here("Output", paste0(dat.title, '_', majgrp.lev[i1], 
                                                   '_vs_', majgrp.lev[i2], '_DEPCCNum.txt')), 
                row.names = F)
    
    # Overlapping (Venn diagram)
    DEPCCIdx.list <- list(CSN = which(idendiffPCC.mat[1, ] == 1), 
                          L2N = which(idendiffPCC.mat[2, ] == 1), 
                          PQN = which(idendiffPCC.mat[3, ] == 1), 
                          QT = which(idendiffPCC.mat[4, ] == 1), 
                          LSCN = which(idendiffPCC.mat[5, ] == 1))
    
    # Venn plot
    plot.name <- paste(dat.title, '_PCCVenn_', 
                       cap.txt, '.png')
    venn.plot <- venn.diagram(x = DEPCCIdx.list, category.names = names(DEPCCIdx.list), 
                              filename = plot.name, output = T, 
                              col = set_col2, fill = set_col2, cat.col = set_col2, 
                              alpha = 0.2, cat.dist = 0.08, cat.fontface = 'bold', 
                              cat.cex = 1.5, cex = 1.5)
  }
}


## c. (Real) Classification ####
library(ROCR)
dat.num <- length(dat.list)
dat.name <- names(dat.list)

fold <- 7

if (dat.title == 'ST002178') {
  preproc <- c('zv', 'pca')
} else {
  preproc <- c('zv', 'corr', 'pca')
}

prob_lev <- 0.5

group <- majGroup
group1 <- minGroup

if (dat.title == 'ST003484') {
  ident_grp <- 'Attached'
} else if (dat.title == 'MTBLS290') {
  ident_grp <- 'Patient'
} else if (dat.title == 'ST002178') {
  ident_grp <- 'Young'
} else {
  ident_grp <- 'YZ'
}

# initialize
N <- nrow(dat.list[[1]])
majgrp.f <- as.factor(majGroup)
majgrp.num <- nlevels(majgrp.f)
majgrp.lev <- levels(majgrp.f)

ident_grp.idx <- which(majgrp.lev == ident_grp)

rep.time <- 100
dat.num <- length(dat.list)
Eval.metrics1 <- matrix(0, 5, dat.num)
pred.res <- list()
for (i in c(1: dat.num)) {
  pred.res[[i]] <- matrix(NA, N, rep.time)
}

n <- 3
for (k in c(1: rep.time)) {
  test.sampIdx.list <- SRsamp(N, n, 
                              group = group, 
                              group1 = group1)
  
  for (j in c(1: n)) {
    test.sampIdx <- test.sampIdx.list[[j]]
    train.sampIdx <- c(1: N)[- test.sampIdx]
    
    train.group <- group[train.sampIdx]
    test.group <- group[test.sampIdx]
    
    for (i in c(1: dat.num)) {
      dat.i <- dat.list[[i]]
      
      dat.tr <- dat.i[train.sampIdx, ]
      dat.ts <- dat.i[test.sampIdx, ]
      
      rf.pred <- CVforClas_Pred(train_dat = dat.tr, 
                                test_dat = dat.ts, 
                                tr.group = train.group, 
                                ts.group = test.group, 
                                method = 'svmRadial', 
                                resamp_method = 'cv', 
                                fold = fold, 
                                preproc = preproc, 
                                repeat_times = 1, 
                                show.progress = T)
      
      pred.i <- rf.pred
      pred.res[[i]][test.sampIdx, k] <- pred.i[, ident_grp.idx]
    }
  }
  cat(k, '\n')
}


pred.probs <- sapply(pred.res, 
                     FUN = function(x) 
                       apply(x, MARGIN = 1, 
                             FUN = function(x) median(x, na.rm = T)))
for (i in c(1: dat.num)) {
  pred.i <- pred.probs[, i]
  eval.metric.i <- getClasMetric(pred.i, group,
                                 prob.lev = prob_lev, 
                                 ident_grp.idx = ident_grp.idx)
  Eval.metrics1[, i] <- eval.metric.i
}


# plot
perf.df <- data.frame(Precision = Eval.metrics1[1, -1], 
                      Recall = Eval.metrics1[2, -1], 
                      F1_score = Eval.metrics1[3, -1], 
                      Accuracy = Eval.metrics1[4, -1], 
                      AUROC = Eval.metrics1[5, -1])

perf.df <- data.frame(Method = factor(dat.name[-1], 
                                      levels = dat.name[-1]), 
                      (perf.df - 0.8) / 0.2)
library(ggradar)
ggradar(perf.df, grid.min = 0, grid.mid = 0.5, grid.max = 1, 
        gridline.min.colour = "grey", 
        gridline.mid.colour = "black", 
        gridline.max.colour = "black", 
        values.radar = c('0.8', '0.9', '1'), background.circle.colour = "white", 
        group.colours = c('#FF66FF', '#FFCC00', '#0066FF', '#00CC66', '#FF3366'))


# Eval 3 - Correction of Dilution Effects ####

## a. Cumulative percentage curves of RSDs ####
### I. Simulated datasets ####
rep.time <- 100
rsds <- seq(0.01, 1, 0.01)

cumFreq.simu <- matrix(0, rep.time, length(rsds))
cumFreq.CSN <- matrix(0, rep.time, length(rsds))
cumFreq.L2N <- matrix(0, rep.time, length(rsds))
cumFreq.PQN <- matrix(0, rep.time, length(rsds))
cumFreq.QT <- matrix(0, rep.time, length(rsds))
cumFreq.LSCN <- matrix(0, rep.time, length(rsds))

for (k in c(1: rep.time)) {
  dat_simu.res <- GenSimuDat(dat.healthy)
  dat_simu <- dat_simu.res$SimuDat
  grp_simu <- dat_simu.res$Para$group
  grp_simu.f <- as.factor(grp_simu)
  grp_simu.num <- nlevels(grp_simu.f)
  grp_simu.lev <- levels(grp_simu.f)
  
  # Unnorm
  for (g in c(1: grp_simu.num)) {
    dat_simu.gi <- dat_simu[grp_simu == grp_simu.lev[g], ]
    m.simu <- colMeans(dat_simu.gi)
    sd.simu <- apply(dat_simu.gi, MARGIN = 2, FUN = sd)
    rsd.simu <- sd.simu / m.simu
    cum.freq <- getCumFreq(rsd.simu, rsds = rsds)
    cumFreq.simu[k, ] <- cumFreq.simu[k, ] + cum.freq
  }
  cumFreq.simu[k, ] <- cumFreq.simu[k, ] / grp_simu.num
  
  # CSN
  dat.nor <- CSN_normalize(dat_simu)
  dat.nor.CSN <- dat.nor$dat.norm
  for (g in c(1: grp_simu.num)) {
    dat_CSN.gi <- dat.nor.CSN[grp_simu == grp_simu.lev[g], ]
    m.CSN <- colMeans(dat_CSN.gi)
    sd.CSN <- apply(dat_CSN.gi, MARGIN = 2, FUN = sd)
    rsd.CSN <- sd.CSN / m.CSN
    cum.freq <- getCumFreq(rsd.CSN, rsds = rsds)
    cumFreq.CSN[k, ] <- cumFreq.CSN[k, ] + cum.freq
  }
  cumFreq.CSN[k, ] <- cumFreq.CSN[k, ] / grp_simu.num
  
  # L2N
  dat.nor <- L2_Normalize(dat_simu)
  dat.nor.L2N <- dat.nor$dat.norm
  for (g in c(1: grp_simu.num)) {
    dat_L2N.gi <- dat.nor.L2N[grp_simu == grp_simu.lev[g], ]
    m.L2N <- colMeans(dat_L2N.gi)
    sd.L2N <- apply(dat_L2N.gi, MARGIN = 2, FUN = sd)
    rsd.L2N <- sd.L2N / m.L2N
    cum.freq <- getCumFreq(rsd.L2N, rsds = rsds)
    cumFreq.L2N[k, ] <- cumFreq.L2N[k, ] + cum.freq
  }
  cumFreq.L2N[k, ] <- cumFreq.L2N[k, ] / grp_simu.num
  
  # PQN
  dat.nor <- PQN_normalize(dat_simu)
  dat.nor.PQN <- dat.nor$dat.norm
  for (g in c(1: grp_simu.num)) {
    dat_PQN.gi <- dat.nor.PQN[grp_simu == grp_simu.lev[g], ]
    m.PQN <- colMeans(dat_PQN.gi)
    sd.PQN <- apply(dat_PQN.gi, MARGIN = 2, FUN = sd)
    rsd.PQN <- sd.PQN / m.PQN
    cum.freq <- getCumFreq(rsd.PQN, rsds = rsds)
    cumFreq.PQN[k, ] <- cumFreq.PQN[k, ] + cum.freq
  }
  cumFreq.PQN[k, ] <- cumFreq.PQN[k, ] / grp_simu.num
  
  # quantile
  dat.nor <- Quantile_normalize(dat_simu)
  dat.nor.QT <- dat.nor$dat.norm
  for (g in c(1: grp_simu.num)) {
    dat_QT.gi <- dat.nor.QT[grp_simu == grp_simu.lev[g], ]
    m.QT <- colMeans(dat_QT.gi)
    sd.QT <- apply(dat_QT.gi, MARGIN = 2, FUN = sd)
    rsd.QT <- sd.QT / m.QT
    cum.freq <- getCumFreq(rsd.QT, rsds = rsds)
    cumFreq.QT[k, ] <- cumFreq.QT[k, ] + cum.freq
  }
  cumFreq.QT[k, ] <- cumFreq.QT[k, ] / grp_simu.num
  
  # LSCN
  dat.nor <- LSCN(dat_simu, ita0 = 0.4, 
                  max_ite = 30, print.txt = F)
  dat.nor.LSCN <- dat.nor$dat.norm
  for (g in c(1: grp_simu.num)) {
    dat_LSCN.gi <- dat.nor.LSCN[grp_simu == grp_simu.lev[g], ]
    m.LSCN <- colMeans(dat_LSCN.gi)
    sd.LSCN <- apply(dat_LSCN.gi, MARGIN = 2, FUN = sd)
    rsd.LSCN <- sd.LSCN / m.LSCN
    cum.freq <- getCumFreq(rsd.LSCN, rsds = rsds)
    cumFreq.LSCN[k, ] <- cumFreq.LSCN[k, ] + cum.freq
  }
  cumFreq.LSCN[k, ] <- cumFreq.LSCN[k, ] / grp_simu.num
  
  cat(paste0(k, '--'))
  
}

# plot 
library(ggplot2)

vals <- c(as.vector(cumFreq.simu), as.vector(cumFreq.LSCN), 
          as.vector(cumFreq.CSN), as.vector(cumFreq.L2N), 
          as.vector(cumFreq.PQN), as.vector(cumFreq.QT))

dat.type <- factor(rep(c('Unnorm', 'LSCN', 'CSN', 
                         'L2N', 'PQN', 'QT'),
                       each = rep.time * length(rsds)),
                   levels = c('Unnorm', 'CSN', 'L2N', 
                              'PQN', 'QT', 'LSCN'))

rsds_ <- rep(rep(rsds, each = rep.time), 6)


# plot
dat.cumFreq <- data.frame(dataset = dat.type, 
                          rsd = rsds_, 
                          cumfreq = vals)

set_col <- c('black', 'magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 5), 1)
set_alp <- c(rep(1, 5), 0.5)

arrow_style <- arrow(length = unit(0.2, 'cm'))
arr_dat <- data.frame(xstart = c(0.4, 0.2), 
                      xend = c(0.31, 0.29), 
                      ystart = c(mean(cumFreq.simu[, 30]), 
                                 mean(cumFreq.LSCN[, 30])), 
                      yend = c(mean(cumFreq.simu[, 30]), 
                               mean(cumFreq.LSCN[, 30])), 
                      datty = c('Unnorm', 'LSCN'))

ggplot() + stat_summary(data = dat.cumFreq, 
                        aes(x = rsd, y = cumfreq, fill = dataset), 
                        fun.data = 'mean_cl_boot', geom = 'ribbon', 
                        alpha = 0.1) + 
  stat_summary(data = dat.cumFreq, 
               aes(x = rsd, y = cumfreq, 
                   color = dataset, linetype = dataset, alpha = dataset), 
               fun = 'mean', geom = 'line', linewidth = 1) + 
  geom_vline(xintercept = 0.3, color = 'cyan', linetype = 2, 
             linewidth = 0.8) + 
  geom_segment(data = arr_dat, 
               aes(x = xstart, y = ystart, 
                   xend = xend, yend = yend, 
                   color = datty), arrow = arrow_style, 
               linewidth = 1) + 
  annotate('text', x = arr_dat$xstart[1] + 0.05, 
           y = arr_dat$ystart[1], label = round(arr_dat$ystart[1], 2), 
           color = 'black', size = 6) + 
  annotate('text', x = arr_dat$xstart[2] - 0.05, 
           y = arr_dat$ystart[2], label = round(arr_dat$ystart[2], 2), 
           color = 'red', size = 6) + 
  scale_color_manual(values = set_col) + 
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_alpha_manual(values = set_alp) + 
  theme_bw() + 
  xlab('RSD') + 
  ylab('Cumulative Frequency') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'inside', 
        legend.position.inside = c(0.85, 0.2), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 18, family = 'sans', face = 'bold'), 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


### II. Real biological datasets ####
rsds <- seq(0.01, 1, 0.01)
dat.name <- names(dat.list)
rsd.res <- getAUCforRSD(dat.list, group = majGroup, 
                        rsds = rsds)
cumFreqs <- rsd.res$cumFreqs
rsdAUCs <- rsd.res$rsdAUCs

vals <- as.vector(cumFreqs)
dat.type <- factor(rep(dat.name, each = length(rsds)),
                   levels = c('Unnorm', 'CSN', 'L2N', 
                              'PQN', 'QT', 'LSCN'))
rsds_ <- rep(rsds, 6)

# plot
dat.cumFreq <- data.frame(dataset = dat.type, 
                          rsd = rsds_, 
                          cumfreq = vals)

set_col <- c('black', 'magenta', 'orange', 'blue', 'green', 'red')
set_line <- c(rep(2, 5), 1)

max_val <- max(dat.cumFreq$cumfreq)

p <- ggplot(data = dat.cumFreq, aes(x = rsd, y = cumfreq)) + 
  geom_line(aes(color = dataset, linetype = dataset), 
            linewidth = 1) + 
  geom_vline(xintercept = 0.3, color = 'cyan', linetype = 2,
             linewidth = 0.8) + 
  annotate('text', x = 0.9,
           y = 0.2 * max_val,
           label = paste0('AUC:', '\n',
                          'Unnorm: ', round(rsdAUCs$Unnorm, 2), '\n',
                          'CSN: ', round(rsdAUCs$CSN, 2), '\n',
                          'L2N: ', round(rsdAUCs$L2N, 2), '\n',
                          'PQN: ', round(rsdAUCs$PQN, 2), '\n',
                          'QT: ', round(rsdAUCs$QT, 2), '\n',
                          'LSCN: ', round(rsdAUCs$LSCN, 2), '\n'),
           color = 'black', size = 6) + 
  scale_color_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  theme_bw() + 
  xlim(0, 1) + ylim(0, 1) + 
  xlab('RSD') + ylab('Cumulative Frequency') + 
  theme(text = element_text(family = 'sans', face = 'bold'), 
        plot.title = element_text(size = 30, hjust = 0.1, 
                                  vjust = -8), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank()) + 
  labs(title = datIdx)

if (datIdx == 'I') {
  p <- p + theme(axis.title.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.x = element_blank(), 
                 axis.title.y = element_text(size = 25, family = 'sans', 
                                             face = "bold"), 
                 axis.text.y = element_text(size = 22, family = 'sans', 
                                            face = "bold", color = 'black'), 
                 legend.position = 'none', 
                 panel.border = element_rect(linewidth = 2), 
                 panel.grid = element_blank())
} else if (datIdx == 'II') {
  p <- p + theme(axis.title = element_blank(), 
                 axis.text = element_blank(), 
                 axis.ticks = element_blank(), 
                 legend.position = 'none', 
                 panel.border = element_rect(linewidth = 2), 
                 panel.grid = element_blank())
} else if (datIdx == 'III') {
  p <- p + theme(axis.title = element_text(size = 25, family = 'sans', 
                                             face = "bold"), 
                 axis.text = element_text(size = 22, family = 'sans', 
                                            face = "bold", color = 'black'), 
                 legend.position = 'none', 
                 panel.border = element_rect(linewidth = 2), 
                 panel.grid = element_blank())
} else {
  p <- p + theme(axis.title.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.title.x = element_text(size = 25, family = 'sans', 
                                             face = "bold"), 
                 axis.text.x = element_text(size = 22, family = 'sans', 
                                            face = "bold", color = 'black'), 
                 legend.position = 'right', 
                 legend.text = element_text(size = 18, family = 'sans', 
                                            face = 'bold'), 
                 legend.title = element_blank(), 
                 panel.border = element_rect(linewidth = 2), 
                 panel.grid = element_blank())
}


## b. (Simulated) Normalization Accuracy ####
### I. Individual Differences ####
rep.time <- 20
n.ratio <- seq(0.05, 0.45, 0.1)

dilt.PCCs.LSCN <- matrix(0, rep.time, length(n.ratio))
dilt.PCCs.CSN <- matrix(0, rep.time, length(n.ratio))
dilt.PCCs.L2N <- matrix(0, rep.time, length(n.ratio))
dilt.PCCs.PQN <- matrix(0, rep.time, length(n.ratio))
dilt.PCCs.QT <- matrix(0, rep.time, length(n.ratio))

for (i in c(1: length(n.ratio))) {
  
  PCCs.LSCN <- rep(0, rep.time)
  PCCs.CSN <- rep(0, rep.time)
  PCCs.L2N <- rep(0, rep.time)
  PCCs.PQN <- rep(0, rep.time)
  PCCs.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, 
                               indiv.sd.ratio = n.ratio[i])
    dat_simu <- dat_simu.res$SimuDat
    diltF.pres <- 1 / dat_simu.res$Para$diltF
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, 
                    max_ite = 30, print.txt = F)
    coef.LSCN <- dat.nor$NORM.COEF
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    coef.CSN <- dat.nor$norm.coef
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    coef.L2N <- dat.nor$norm.coef
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    coef.PQN <- dat.nor$norm.coef
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    coef.QT <- dat.nor.QT / dat_simu
    pccs <- rep(0, ncol(coef.QT))
    for (j in c(1: ncol(coef.QT))) {
      pccs[j] <- cor(diltF.pres, coef.QT[, j])
    }
    
    PCCs.LSCN[k] <- cor(diltF.pres, coef.LSCN)
    PCCs.CSN[k] <- cor(diltF.pres, coef.CSN)
    PCCs.L2N[k] <- cor(diltF.pres, coef.L2N)
    PCCs.PQN[k] <- cor(diltF.pres, coef.PQN)
    PCCs.QT[k] <- mean(pccs)

  }
  
  
  dilt.PCCs.LSCN[, i] <- PCCs.LSCN
  dilt.PCCs.CSN[, i] <- PCCs.CSN
  dilt.PCCs.L2N[, i] <- PCCs.L2N
  dilt.PCCs.PQN[, i] <- PCCs.PQN
  dilt.PCCs.QT[, i] <- PCCs.QT
  
  cat(i, '\n')
}

# plot
vals <- c(as.vector(dilt.PCCs.LSCN), as.vector(dilt.PCCs.CSN), 
          as.vector(dilt.PCCs.L2N), as.vector(dilt.PCCs.PQN), 
          as.vector(dilt.PCCs.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm', 
                         'PQN', 'Quantile'),
                       each = rep.time * length(n.ratio)),
                   levels = c('LSCN', 'CSN', 'L2-norm', 
                              'PQN', 'Quantile'))

nr <- rep(rep(n.ratio, each = rep.time), 5)

dat.dfPCC <- data.frame(dattype = dat.type, 
                        NR = nr, 
                        dfPCCs = vals)

set_col <- c('red', 'magenta', 'orange', 'blue', 'green')
set_line <- c(1, rep(2, 4))


ggplot() + 
  geom_smooth(data = dat.dfPCC, aes(x = NR, y = dfPCCs,
                                    color = dattype, fill = dattype, 
                                    linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = n.ratio) + 
  theme_bw() + 
  xlab('Individual Differences (r)') + 
  ylab('PCC of Dilution Factors') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

### II. Inter-group Heterogeneity ####
rep.time <- 20
grpVar.ratio <- seq(0.05, 0.45, 0.1)

dilt.PCCs.LSCN <- matrix(0, rep.time, length(grpVar.ratio))
dilt.PCCs.CSN <- matrix(0, rep.time, length(grpVar.ratio))
dilt.PCCs.L2N <- matrix(0, rep.time, length(grpVar.ratio))
dilt.PCCs.PQN <- matrix(0, rep.time, length(grpVar.ratio))
dilt.PCCs.QT <- matrix(0, rep.time, length(grpVar.ratio))

for (i in c(1: length(grpVar.ratio))) {

  PCCs.LSCN <- rep(0, rep.time)
  PCCs.CSN <- rep(0, rep.time)
  PCCs.L2N <- rep(0, rep.time)
  PCCs.PQN <- rep(0, rep.time)
  PCCs.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, 
                               grpVar.ratio = grpVar.ratio[i])
    dat_simu <- dat_simu.res$SimuDat
    diltF.pres <- 1 / dat_simu.res$Para$diltF
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, 
                    max_ite = 30, print.txt = F)
    coef.LSCN <- dat.nor$NORM.COEF
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    coef.CSN <- dat.nor$norm.coef
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    coef.L2N <- dat.nor$norm.coef
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    coef.PQN <- dat.nor$norm.coef
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    coef.QT <- dat.nor.QT / dat_simu
    pccs <- rep(0, ncol(coef.QT))
    for (j in c(1: ncol(coef.QT))) {
      pccs[j] <- cor(diltF.pres, coef.QT[, j])
    }
    
    PCCs.LSCN[k] <- cor(diltF.pres, coef.LSCN)
    PCCs.CSN[k] <- cor(diltF.pres, coef.CSN)
    PCCs.L2N[k] <- cor(diltF.pres, coef.L2N)
    PCCs.PQN[k] <- cor(diltF.pres, coef.PQN)
    PCCs.QT[k] <- mean(pccs)
    
    cat(k, ' ')
  }
  
  dilt.PCCs.LSCN[, i] <- PCCs.LSCN
  dilt.PCCs.CSN[, i] <- PCCs.CSN
  dilt.PCCs.L2N[, i] <- PCCs.L2N
  dilt.PCCs.PQN[, i] <- PCCs.PQN
  dilt.PCCs.QT[, i] <- PCCs.QT
  
  cat('\n', i, '\n')
}

# plot
vals <- c(as.vector(dilt.PCCs.LSCN), as.vector(dilt.PCCs.CSN),
          as.vector(dilt.PCCs.L2N), as.vector(dilt.PCCs.PQN),
          as.vector(dilt.PCCs.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm',
                         'PQN', 'Quantile'),
                       each = rep.time * length(grpVar.ratio)),
                   levels = c('LSCN', 'CSN', 'L2-norm', 
                              'PQN', 'Quantile'))

var.ratio <- rep(rep(grpVar.ratio, each = rep.time), 5)

dat.dfPCC <- data.frame(dattype = dat.type, 
                        grpVar.ratio = var.ratio, 
                        dfPCCs = vals)

set_col <- c('red', 'magenta', 'orange', 'blue', 'green')
set_line <- c(1, rep(2, 4))


ggplot() + 
  geom_smooth(data = dat.dfPCC, aes(x = grpVar.ratio, y = dfPCCs,
                                    color = dattype, fill = dattype, 
                                    linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = grpVar.ratio) + 
  theme_bw() + 
  xlab('Inter-group Heterogeneity (k)') + 
  ylab('PCC of Dilution Factors') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


### III. Dilution Effects ####
rep.time <- 20
sv.m <- seq(10, 90, 20)

dilt.PCCs.LSCN <- matrix(0, rep.time, length(sv.m))
dilt.PCCs.CSN <- matrix(0, rep.time, length(sv.m))
dilt.PCCs.L2N <- matrix(0, rep.time, length(sv.m))
dilt.PCCs.PQN <- matrix(0, rep.time, length(sv.m))
dilt.PCCs.QT <- matrix(0, rep.time, length(sv.m))

for (i in c(1: length(sv.m))) {
  
  PCCs.LSCN <- rep(0, rep.time)
  PCCs.CSN <- rep(0, rep.time)
  PCCs.L2N <- rep(0, rep.time)
  PCCs.PQN <- rep(0, rep.time)
  PCCs.QT <- rep(0, rep.time)
  
  for (k in c(1: rep.time)) {
    dat_simu.res <- GenSimuDat(dat.healthy, 
                               diltF.m = sv.m[i])
    dat_simu <- dat_simu.res$SimuDat
    diltF.pres <- 1 / dat_simu.res$Para$diltF
    
    # LSCN
    dat.nor <- LSCN(dat_simu, ita0 = 0.4, 
                    max_ite = 30, print.txt = F)
    coef.LSCN <- dat.nor$NORM.COEF
    
    # CSN
    dat.nor <- CSN_normalize(dat_simu)
    coef.CSN <- dat.nor$norm.coef
    
    # L2N
    dat.nor <- L2_Normalize(dat_simu)
    coef.L2N <- dat.nor$norm.coef
    
    # PQN
    dat.nor <- PQN_normalize(dat_simu)
    coef.PQN <- dat.nor$norm.coef
    
    # quantile
    dat.nor <- Quantile_normalize(dat_simu)
    dat.nor.QT <- dat.nor$dat.norm
    coef.QT <- dat.nor.QT / dat_simu
    pccs <- rep(0, ncol(coef.QT))
    for (j in c(1: ncol(coef.QT))) {
      pccs[j] <- cor(diltF.pres, coef.QT[, j])
    }
    
    PCCs.LSCN[k] <- cor(diltF.pres, coef.LSCN)
    PCCs.CSN[k] <- cor(diltF.pres, coef.CSN)
    PCCs.L2N[k] <- cor(diltF.pres, coef.L2N)
    PCCs.PQN[k] <- cor(diltF.pres, coef.PQN)
    PCCs.QT[k] <- mean(pccs)
  }
  
  
  dilt.PCCs.LSCN[, i] <- PCCs.LSCN
  dilt.PCCs.CSN[, i] <- PCCs.CSN
  dilt.PCCs.L2N[, i] <- PCCs.L2N
  dilt.PCCs.PQN[, i] <- PCCs.PQN
  dilt.PCCs.QT[, i] <- PCCs.QT
  
  cat(i, '\n')
}

# plot
vals <- c(as.vector(dilt.PCCs.LSCN), as.vector(dilt.PCCs.CSN),
          as.vector(dilt.PCCs.L2N), as.vector(dilt.PCCs.PQN),
          as.vector(dilt.PCCs.QT))

dat.type <- factor(rep(c('LSCN', 'CSN', 'L2-norm',
                         'PQN', 'Quantile'),
                       each = rep.time * length(sv.m)),
                   levels = c('LSCN', 'CSN', 'L2-norm', 
                              'PQN', 'Quantile'))

svms <- rep(rep(sv.m, each = rep.time), 5)

dat.dfPCC <- data.frame(dattype = dat.type, 
                        SVM = svms, 
                        dfPCCs = vals)

set_col <- c('red', 'magenta', 'orange', 'blue', 'green')
set_line <- c(1, rep(2, 4))


ggplot() + 
  geom_smooth(data = dat.dfPCC, aes(x = SVM, y = dfPCCs,
                                    color = dattype, fill = dattype, 
                                    linetype = dattype),
              method = lm, formula = y ~ poly(x, 2),
              linewidth = 1, alpha = 0.1) +
  scale_color_manual(values = set_col) +
  scale_fill_manual(values = set_col) + 
  scale_linetype_manual(values = set_line) + 
  scale_x_continuous(breaks = sv.m) + 
  theme_bw() + 
  xlab('Dilution Effects (m)') + 
  ylab('PCC of Dilution Factors') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())


## c. (Real) Correlation between normalization factors ####
# normalize factors to median 1
coef.LSCN <- coef.LSCN / median(coef.LSCN)
coef.CSN <- coef.CSN / median(coef.CSN)
coef.L2N <- coef.L2N / median(coef.L2N)
coef.PQN <- coef.PQN / median(coef.PQN)
coef.QT <- coef.QT / median(coef.QT)
all.coef <- cbind(coef.LSCN, coef.CSN, coef.L2N, 
                  coef.PQN, coef.QT)
method_name <- c('LSCN', 'CSN', 'L2N', 'PQN', 'Quantile')
coef.cor <- rep(0, (length(method_name) - 1))
for (i in c(1: (length(method_name) - 1))) {
  coef.cor[i] <- cor(all.coef[, 1], all.coef[, (i + 1)])
  
  coef.df <- data.frame(x = all.coef[, 1], 
                        y = all.coef[, (i + 1)])
  
  # plot
  x.len <- max(coef.df$x) + 0.05 - min(coef.df$x)
  y.len <- max(coef.df$y) + 0.05 - min(coef.df$y)
  
  dx <- x.len / 8
  dy <- y.len / 18
  
  p.i <- ggplot(coef.df, aes(x = x, y = y)) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2, 
                linewidth = 0.5) + 
    geom_point(size = 3, color = 'blue', alpha = 0.5) + 
    annotate(geom = 'text', 
             x = max(coef.df$x) + 0.05 - dx, 
             y = min(coef.df$y) + dy, 
             label = paste0('r = ', round(coef.cor[i], 2)), 
             size = 8, color = 'black', 
             family = 'sans', fontface = 'bold') + 
    xlim(min(coef.df$x), max(coef.df$x) + 0.05) + 
    ylim(min(coef.df$y), max(coef.df$y) + 0.05) + 
    theme_bw() + 
    xlab(paste0('Normalization Factor (', method_name[1], ')')) + 
    ylab(paste0('Normalization Factor (', method_name[i + 1], ')')) + 
    theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
          axis.text = element_text(size = 22, family = 'sans', face = "bold", 
                                   color = 'black'), 
          legend.position = 'none', 
          panel.border = element_rect(linewidth = 2), 
          panel.grid = element_blank())
  print(p.i)
    
}


## d. (Real) Run Time and Convergent Curve ####
# CSN
dat.nor <- CSN_normalize(dat)
run_time.CSN <- dat.nor$RunTime
run_time.CSN

# L2N
dat.nor <- L2_Normalize(dat)
run_time.L2N <- dat.nor$RunTime
run_time.L2N

# PQN
dat.nor <- PQN_normalize(dat)
run_time.PQN <- dat.nor$RunTime
run_time.PQN

# Quantile
dat.nor <- Quantile_normalize(dat)
run_time.QT <- dat.nor$RunTime
run_time.QT

### Speed (ita0 = 0.4) ####
# LSCN
if (dat.title == 'teadata_DMSO') {
  sq <- 0.9
} else {
  sq <- NA
}
dat.nor <- LSCN(dat, doLog = F, ita0 = 0.4, max_ite = 50, 
                simi.quant = sq, print.txt = T)
run_time.LSCN <- dat.nor$RunTime
run_time.LSCN
conv <- dat.nor$conv.vals

conv.df <- data.frame(iteration = c(1: length(conv)), 
                      value = conv)

ggplot(conv.df) + geom_line(aes(x = iteration, y = value)) + 
  geom_point(aes(x = iteration, y = value), shape = 21) + 
  theme_bw() + 
  xlab('Iteration') + 
  ylab('Value') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

### Accuracy (ita0 = 0.2) ####
# LSCN
if (dat.title == 'teadata_DMSO') {
  sq <- 0.9
} else {
  sq <- NA
}
dat.nor <- LSCN(dat, doLog = F, ita0 = 0.2, max_ite = 100, 
                simi.quant = sq, print.txt = T)
run_time.LSCN <- dat.nor$RunTime
run_time.LSCN
conv <- dat.nor$conv.vals

conv.df <- data.frame(iteration = c(1: length(conv)), 
                      value = conv)

ggplot(conv.df) + geom_line(aes(x = iteration, y = value)) + 
  geom_point(aes(x = iteration, y = value), shape = 21) + 
  theme_bw() + 
  xlab('Iteration') + 
  ylab('Value') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

