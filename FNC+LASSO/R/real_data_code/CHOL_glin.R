library(MASS)
library(hdi)
library(glmnet)
library(bigsnpr)
library(bigstatsr)
library(data.table)
library(magrittr)
library(tidyverse)
library(glue)
library(data.table)
library(foreach)
library(doParallel)
library(glinternet)

fnp.est <- function(z, t, s, side) {
  p <- length(z)
  if (side==1) {rej=sum(z>t)}
  if (side==2) {rej=sum(abs(z)>t)}
  fnp <- 1 - rej/s + side*(p-s)*pnorm(-t, 0, 1)/s
  if (fnp<0) {fnp=0} else { if (fnp>1) {fnp=1} }
  return(fnp)
}

fnp.opt.series <- function(z, epsilon, s, side) {   # input FNP
  if (side==1) {z.order = sort(z, decreasing = T)}
  if (side==2) {z.order = sort(abs(z), decreasing = T)}
  j <- 1
  epsilon <- sort(epsilon, decreasing = T)
  j_hat <- c()
  t_hat <- c()
  for (i in 1:length(epsilon)) {
    while (fnp.est(z, z.order[j], s, side) > epsilon[i]) {
      j <- j + 1
      print(c(i,j))
    }
    j_hat[i] <- j-1
    t_hat[i] <- z.order[j]
  }
  return(list(epsilon=rev(epsilon), t_hat=rev(t_hat), j_hat=rev(j_hat)))
}

MR05 = function(z, c05){
  n <- length(z)
  z_order = sort(abs(z), decreasing = T, index.return = T)$x
  p_order = (1- pnorm(z_order))*2
  H = rep(0, n)
  ind = 1 : n
  H = (ind / n - p_order - c05 * sqrt(p_order/2)) / (1 - p_order)
  pi = max(H[2 : floor(n / 2)], 0)
  return(pi)
}


args = commandArgs(trailingOnly=TRUE)
k = as.numeric(args[1])
loopNum <- 50

setwd('/home3/vvenka23/FNC-Lasso/RedoPaperCoLaus/CHOL0.8')
tmpfile <- tempfile()
snp_readBed("/home3/vvenka23/FNC-Lasso/RedoPaperCoLaus/dataHg38/hg38/fixFiles3/merged_QC_8.bed",backingfile = tmpfile)
obj.bigSNP <- snp_attach(paste0(tmpfile, ".rds"))

samples.target <- fread('/home3/vvenka23/FNC-Lasso/RedoPaperCoLaus/dataHg38/hg38/fixFiles3/merged_QC_8.fam')[,1:2]
colnames(samples.target) <- c("FID", "IID") 

pc10 <- fread('/home3/DrYao/CoLaus/data/top10PC.txt')
pc10.sub <- merge(samples.target, pc10, by="IID", sort = FALSE)
phe <- fread('/home3/DrYao/CoLaus/data/CoLausPheno.csv')[,c('Global_ID', "SubjectAlias","AGE", "SEX","HT_cm_","CHOL", "TRIG")]
sex <- rep(0, dim(phe)[1])
sex[phe$SEX=='M'] <- 1
phe$SEX <- sex

phenotype <- merge(pc10.sub, phe, by.x="IID", by.y="SubjectAlias", sort = FALSE)
phenotype <- phenotype[!is.na(phenotype$CHOL),]
samples.idx <- samples.target$IID %in% phenotype$IID # matched samples indice


G <- obj.bigSNP$genotypes[samples.idx, ]
code <- rep(NA_real_, 256)
code[1:3] <- c(0, 1, 2)
G <- FBM.code256(dim(G)[1], dim(G)[2], code, init = G)

sumstats <- bigreadr::fread2('/home3/vvenka23/FNC-Lasso/RedoPaperCoLaus/CHOL0.8/20686565-GCST000760-EFO_0004574.h.tsv.gz')[,c("hm_rsid", "chromosome", "hm_pos","hm_effect_allele","hm_other_allele", "p_value", "beta", "standard_error")]
names(sumstats) <- c("rsid", "chr", "pos", "a0", "a1", "p", "beta", "se")
sumstats$n_eff <- 100184
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map, join_by_pos = TRUE)

beta <- rep(NA, ncol(G))
beta[info_snp$`_NUM_ID_`] <- info_snp$beta
lpval <- rep(NA, ncol(G))
lpval[info_snp$`_NUM_ID_`] <- -log10((1-pnorm(abs(info_snp$beta/info_snp$se)))*2)
lpval[info_snp$`_NUM_ID_`] <- -log10(info_snp$p)
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

p <-  dim(info_snp)[1]
n <- dim(G)[1] 
t <- rep(NA, ncol(G))
t[info_snp$`_NUM_ID_`] <- info_snp$beta/info_snp$se
tt <- abs(t)
t.ind <- sort(tt, decreasing=T, index.return=T, na.last = TRUE)$ix


library(MASS)
library(hdi) 
library(glmnet)
library(bigsnpr)
library(data.table)
library(magrittr)
library(tidyverse)
  
set.seed(0)
nseed <- ceiling(runif(loopNum)*10000)
set.seed(nseed[k])

train.ratio <- 0.5
ind.train <- sample(nrow(G), round(train.ratio*nrow(G)))
ind.test <- setdiff(rows_along(G), ind.train)

y <- phenotype$CHOL
pc <- phenotype[,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "AGE", "SEX")]
y1 <- y[ind.train]
y2 <- y[ind.test]
sst <- sum((y2-mean(y2))^2)

pc1 <- pc[ind.train,]
pc2 <- pc[ind.test,]
n.test <- length(ind.test)
pi = MR05(scale(t[!is.na(t)]),  0.002976302) #CHOL
  epsilon <- 0.02*c(1:40) 
  rss.lasso <- rep(NA, length(epsilon))
  lasso.tune <-  matrix(0, length(epsilon), p)
  lasso.final <- matrix(0, 1, p)
  R.lasso <- rep(0, length(epsilon))
  n.lasso <- rep(0, length(epsilon))
  fnc.t1 <- Sys.time()
  fnp.tune <- fnp.opt.series(scale(t[!is.na(t)]), epsilon, s = ceiling(pi*p), side = 2)$j_hat
  
  K.fold <- sample(rep_len(1:5, length(ind.train)))
  for (tune in 1:length(epsilon)){
    screen.ind = t.ind[1:fnp.tune[tune]]
    X.fnc <- FBM(n, length(screen.ind), init = G[, screen.ind])
    mod.fnc <- big_spLinReg(X.fnc , y[ind.train], K = 5, ind.train = ind.train, covar.train = covar_from_df(pc1), ind.sets = K.fold)
    
    pred.lasso = predict(mod.fnc, X.fnc, ind.train, covar.row = covar_from_df(pc1))
    rss.lasso[tune] = sum((y1-pred.lasso)^2)/(sum((y1-mean(y1))^2)*(length(ind.train)-summary(mod.fnc)$nb_var-1))
    
    pred.fnc <- predict(mod.fnc, X.fnc, ind.test, covar.row = covar_from_df(pc2))
    
    n.lasso[tune] <- summary(mod.fnc)$nb_var
    R.lasso[tune] <- 1-sum((y2-pred.fnc)^2)/sst
  }
  fnc.t2 <- Sys.time()
  fnc.t <- fnc.t2 - fnc.t1
  R.fnc.lasso <- R.lasso[which.min(rss.lasso)] 
  n.fnc.lasso <- n.lasso[which.min(rss.lasso)] 
  aic.fnc.lasso <- 2*n.fnc.lasso + n.test * log((1-R.fnc.lasso)*sst/n.test)

#############
opt.tune <- which.min(rss.lasso)
  opt.idx <- t.ind[1: fnp.tune[opt.tune]]
  
  X.opt <- FBM(n, length(opt.idx), init = G[, opt.idx])
  mod.opt <- big_spLinReg(X.opt , y[ind.train], K = 5, ind.train = ind.train, covar.train = covar_from_df(pc1), ind.sets = K.fold)
  
  opt.var <- (mod.opt[[1]][[1]][["beta"]]!=0 | mod.opt[[1]][[2]][["beta"]]!=0 | 
                mod.opt[[1]][[3]][["beta"]]!=0 |  mod.opt[[1]][[4]][["beta"]]!=0  | mod.opt[[1]][[5]][["beta"]]!=0)
  opt.var <- opt.var[1:(length(opt.var)-12)]
  
  dim(X.opt[ind.train, opt.var])
  
  cv.glin <- glinternet.cv(cbind(X.opt[ind.train, opt.var], pc1), y[ind.train], nFolds = 5, numCores = 3,
                           numLevels=c(rep(1, sum(opt.var)), rep(1, 11), 2), verbose = TRUE, maxIter = 300)
  
  pred.lasso = predict(cv.glin, cbind(X.opt[ind.train, opt.var], pc1), type = "response", lambdaType="lambdaHat")
  pred.fnc <- predict(cv.glin, cbind(X.opt[ind.test, opt.var], pc2), type = "response", lambdaType="lambdaHat")
  
  coefs <- coef(cv.glin)

  n.glin <- length(coefs$mainEffects$cat) +  length(coefs$mainEffects$cont)
  n.para.glin <- length(coefs$mainEffects$cat) + length(coefs$mainEffects$cont) + length(coefs$interactions$catcont) + length(coefs$interactions$contcont)
  R.glin <- 1-sum((y2-pred.fnc)^2)/sst
  aic.glin <- 2* n.glin + n.test * log(sum((y2-pred.fnc)^2)/n.test)

 
write.table(cbind(k,R.glin,n.glin,aic.glin),file="glin",quote=FALSE,row.names=FALSE,col.names= !file.exists("glin"),append=T) 
 
