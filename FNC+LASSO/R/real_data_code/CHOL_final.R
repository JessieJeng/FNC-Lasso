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

setwd('/Users/yifeihu/newmatching')

snp_readBed("merged_QC.bed")
obj.bigSNP <- snp_attach("merged_QC.rds")

samples.target <- fread('merged_QC.fam')[,1:2]
colnames(samples.target) <- c("FID", "IID") 

pc10 <- fread('top10PC.txt')
pc10.sub <- merge(samples.target, pc10, by="IID", sort = FALSE)
phe <- fread('Pheno.csv')[,c('Global_ID', "SubjectAlias","AGE", "SEX","HT_cm_","CHOL", "TRIG")]
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



sumstats <- 
  bigreadr::fread2('~/Desktop/Project3/CHOL/20686565-GCST000760-EFO_0004574.h.tsv.gz')[,
  c("hm_rsid", "chromosome", "hm_pos", "hm_effect_allele", "hm_other_allele", "p_value", "beta", "standard_error")]


names(sumstats) <- c("rsid", "chr", "pos","a0", "a1","p", "beta", "se")
sumstats$n_eff <- 100184

map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map, join_by_pos = FALSE)


beta <- rep(NA, ncol(G))
beta[info_snp$`_NUM_ID_`] <- info_snp$beta
lpval <- rep(NA, ncol(G))
lpval[info_snp$`_NUM_ID_`] <- -log10((1-pnorm(abs(info_snp$beta/info_snp$se)))*2)
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos


p <-  dim(info_snp)[1]
n <- dim(G)[1] 
t <- rep(NA, ncol(G))
t[info_snp$`_NUM_ID_`] <- info_snp$beta/info_snp$se
tt <- abs(t)
t.ind <- sort(tt, decreasing=T, index.return=T, na.last = TRUE)$ix

X.just.lasso <- FBM(n, p, init=G[,info_snp$`_NUM_ID_`])


f <- function(k) {
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
  ##################################################################################################
  screen.ind.sis <-  t.ind[1:length(ind.train)-1]
  set.seed(0)
  K.fold <- sample(rep_len(1:5, length(ind.train)))
  X.sis <- FBM(n, length(screen.ind.sis), init = G[, screen.ind.sis])
  mod.sis <- big_spLinReg(X.sis, y[ind.train], K=5, ind.train = ind.train, ind.sets = K.fold, covar.train = covar_from_df(pc1))
  pred.sis <- predict(mod.sis, X.sis, ind.test, covar.row = covar_from_df(pc2))
  
  R.sis.lasso <- 1-sum((y2-pred.sis)^2)/sst 
  n.sis.lasso <- summary(mod.sis)$nb_var
  
  rss.sis.lasso <- sum((y2-pred.sis)^2) 
  aic.sis.lasso <- 2 * n.sis.lasso + n.test * log(rss.sis.lasso/n.test)
  ##################################################################################################
  pi = MR05(scale(t[!is.na(t)]), 0.007304347)
  epsilon <- 0.02*c(1:40) 
  rss.lasso <- rep(NA, length(epsilon))
  lasso.tune <-  matrix(0, length(epsilon), p)
  lasso.final <- matrix(0, 1, p)
  R.lasso <- rep(0, length(epsilon))
  n.lasso <- rep(0, length(epsilon))
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
  
  R.fnc.lasso <- R.lasso[which.min(rss.lasso)] 
  n.fnc.lasso <- n.lasso[which.min(rss.lasso)] 
  aic.fnc.lasso <- 2*n.fnc.lasso + n.test * log((1-R.fnc.lasso)*sst/n.test)
  ##################################################################################################
  mod.just.lasso <- big_spLinReg(X.just.lasso, y[ind.train], K = 5, covar.train = covar_from_df(pc1), ind.train = ind.train, ind.sets = K.fold)
  pred.just.lasso <- predict(mod.just.lasso, X.just.lasso, ind.test, covar.row = covar_from_df(pc2))
  
  R.just.lasso <- 1-sum((y2-pred.just.lasso)^2)/sst 
  n.just.lasso <- summary(mod.just.lasso)$nb_var
  aic.just.lasso <- 2*n.just.lasso + n.test * log((1-R.just.lasso)*sst/n.test)
  ##################################################################################################=
  idx.match <- info_snp$`_NUM_ID_`
  t.lm <- rep(NA, ncol(G))
  test <- big_univLinReg(
    G,
    y1,
    ind.train = ind.train,
    ind.col = info_snp$`_NUM_ID_`,
    covar.train = covar_from_df(pc1),
    ncores = 1
  )
  t.lm[info_snp$`_NUM_ID_`] <- test$score
  tt.lm <- abs(t.lm)
  t.lm.ind <- sort(tt.lm, decreasing=T, index.return=T, na.last = TRUE)$ix
  screen.lm.ind.sis <-  t.lm.ind[1:length(ind.train)-1]
  X.lm.sis <- FBM(n, length(screen.lm.ind.sis), init = G[, screen.lm.ind.sis])
  mod.lm.sis <- big_spLinReg(X.lm.sis, y[ind.train], K=5, ind.train = ind.train, covar.train = covar_from_df(pc1), ind.sets = K.fold)
  n.lm.sis <- summary(mod.lm.sis)$nb_var
  pred.lm.sis <- predict(mod.lm.sis, X.lm.sis, ind.test, covar.row = covar_from_df(pc2))
  R.lm.sis <- 1-sum((y2-pred.lm.sis)^2)/sst 
  aic.lm.sis <- 2*n.lm.sis + n.test * log(sum((y2-pred.lm.sis)^2)/n.test)
  ##################################################################################################
  all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train,
                                lpS = lpval, exclude = which(is.na(lpval)),
                                ncores = 3)
  multi_PRS <- snp_grid_PRS(G, all_keep, beta, lpval, ind.row = ind.train,
                            n_thr_lpS = 25, ncores = 3,
                            grid.lpS.thr = 0.9999 * seq_log(max(0.1, 
                                                                min(lpval, na.rm = TRUE)), 
                                                            max(lpval[lpval!=Inf], na.rm =TRUE), 25))
  grid <- attr(all_keep, "grid") %>%
    mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = row_number()) %>%
    unnest(cols = "thr.lp")
  s <- nrow(grid)
  
  grid$t <- big_apply(multi_PRS, 
                      a.FUN = function(X, ind, s, y.train) {
                        single_PRS <- X[, ind]
                        dt = data.frame(y=y.train, single_PRS, pc1)
                        lm.fit = lm(y~., data = dt)
                        (coef(summary(lm.fit))[, "t value"][2])^2
                      }, ind = 1:s, s = s, y.train = y1, a.combine = 'c', block.size = 1, ncores = 1)
  
  max_prs <- grid[which.max(grid$t), ]
  ind.keep <- unlist(purrr::map(all_keep, max_prs$id))
  ind.ct <- ind.keep[lpval[ind.keep] > max_prs$thr.lp]
  opt.prs <- snp_PRS(G, beta[ind.keep], ind.keep = ind.keep, ind.test = ind.test,
                     lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp)
  dt <- data.frame(y2, opt.prs, pc2)
  lm.fit <- lm(y2~., data = dt)
  R.ct.prs <- summary(lm.fit)$r.squared
  n.ct.prs <- length(ind.ct)
  
  ssr.ct <- sum((lm.fit$residuals)^2)
  n.test <- length(ind.test)
  aic.ct.prs <- 2*n.ct.prs + n.test*log(ssr.ct/n.test)
  ##################################################################################################
  tmp <- tempfile(tmpdir = "tmp-data")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
  # Initialize variables for storing the LD score and LD matrix
  corr <- NULL
  ld <- NULL
  # We want to know the ordering of samples in the bed file 
  fam.order <- NULL
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
  # calculate LD
  for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(
      G,
      ind.col = ind.chr2,
      ncores = NCORES,
      infos.pos = POS2[ind.chr2],
      size = 3 / 1000
    )
    if (chr == 1) {
      ld <- Matrix::colSums(corr0^2)
      corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }
  }
  df_beta <- info_snp[, c("beta", "se", "n_eff")]
  colnames(df_beta) <- c("beta", "beta_se", "n_eff") 
  ldsc <- snp_ldsc(   ld, 
                      length(ld), 
                      chi2 = (df_beta$beta / df_beta$beta_se)^2,
                      sample_size = df_beta$n_eff, 
                      blocks = NULL)
  h2_est <- abs(ldsc[["h2"]])
  p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
  h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
  grid.param <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(TRUE))
  # Get adjusted beta from grid model
  beta_grid <- snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)
  # Obtain model PRS
  pred_grid <- big_prodMat(G, beta_grid, ind.col = info_snp$`_NUM_ID_`)
  grid.param$score <- big_univLinReg(as_FBM(pred_grid[ind.train, ]), y1, covar.train = covar_from_df(pc1))$score
  grid.param %>%
    mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
    arrange(desc(score)) %>%
    mutate_at(c("score", "sparsity"), round, digits = 3) %>%
    slice(1:10)
  best_grid_sp <- grid.param %>%
    mutate(id = row_number()) %>%
    filter(sparse) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]
  pred_sp <- big_prodVec(G, best_grid_sp, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
  dt.sp <- data.frame(y2, pred_sp, pc2)
  lm.fit.sp  <- lm(y2~., data = dt.sp)
  R.ld.prs <- summary(lm.fit.sp)$r.squared  
  n.ld.prs <- sum(best_grid_sp!=0)
  aic.ld.prs <- 2*n.ld.prs + n.test*log((1-R.ld.prs)*sst/n.test)
  
  return(matrix(c(R.fnc.lasso, R.sis.lasso, R.just.lasso, R.lm.sis, R.ct.prs, R.ld.prs,
                  n.fnc.lasso, n.sis.lasso, n.just.lasso, n.lm.sis, n.ct.prs, n.ld.prs,
                  aic.fnc.lasso, aic.sis.lasso, aic.just.lasso, aic.lm.sis, aic.ct.prs, aic.ld.prs
                  ), nrow = 18, ncol = 1))
}



loopNum <- 4
cl <- makeCluster(4)
registerDoParallel(cl)
data.out <- matrix(0, 1, 36)
data <- foreach(k = 1:loopNum, .combine = "cbind") %dopar% f(k)
for(j in 1:18) {
  data.out[1, (2*(j-1)+1)] <- mean(data[j, ])
  data.out[1, (2*(j-1)+2)] <- sd(data[j, ])
}

#write.csv(data.out, "valid_CHOL_sum1.csv")
#write.csv(data, "valid_CHOL_data1.csv")

stopImplicitCluster()
stopCluster(cl)



