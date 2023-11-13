########################################################################################################################
############### Transfer Learning with False Negative Control Improves Polygenic Risk Prediction #######################
############### (for cluster computing) ################################################################################
########################################################################################################################

##################################
#### Part 0. Useful Packages
##################################
# check if all necessary packages are installed
list.of.packages <- c("MASS", "hdi", "glmnet", "data.table", "magrittr", "tidyverse", "foreach", "doParallel", "bigsnpr", "glue")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# attach all packages
library(MASS)
library(hdi) 
library(glmnet)
library(data.table)
library(magrittr)
library(tidyverse)
library(glue)
## package for C+T and LDpred methods
library(bigsnpr) 
## packages for cluster computing
library(foreach)
library(doParallel)
library(dplyr)
##################################
#### Part 1. Useful Functions
##################################
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
    }
    j_hat[i] <- j-1
    t_hat[i] <- z.order[j]
  }
  return(list(epsilon=rev(epsilon), t_hat=rev(t_hat), j_hat=rev(j_hat))) 
}

MR05 <- function(z, c05){
  n <- length(z)
  z_order <- sort(z, decreasing = T, index.return = T)$x
  p_order <- (1- pnorm(z_order))
  H <- rep(0, n)
  ind <- 1 : n
  H <- (ind / n - p_order - c05 * sqrt(p_order/2)) / (1 - p_order)
  pi <- max(H[2 : floor(n / 2)], 0)
  return(pi)
}


##################################
#### Part 2. Data Import
##################################
# We need data file: chr21.QC.bed #
setwd('~/Desktop/temp21/')

# target data: exact the genotype matrix of target data from .bed file
system(glue("rm chr21.QC.bk"))
snp_readBed("chr21.QC.bed")
obj.bigSNP <- snp_attach("chr21.QC.rds")
str(obj.bigSNP)
## set size of target data
n_target <- 2000 # training + testing
## recode genotype matrix for 
G <- obj.bigSNP$genotypes
set.seed(1)
sub.sample <-  sample(1 : dim(G)[1], n_target, replace=F)
code <- rep(NA_real_, 256)
code[1:3] <- c(0, 1, 2)
G <- FBM.code256(n_target, dim(G)[2], code, init = G[sub.sample,])
genotype <- as.data.frame(G[1:dim(G)[1], 1:dim(G)[2]])
colnames(genotype) <- obj.bigSNP$map$marker.ID
## extract chr information for PRS methods
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
RSID <- obj.bigSNP$map$marker.ID
dat <- as.data.frame(cbind(CHR,POS))
colnames(dat) <- c("chr","pos")
POS.hg19 <- snp_modifyBuild(info_snp=dat, liftOver = "liftOver",from = "hg18",to = "hg19")$pos
remove.idx <- RSID[which(is.na(POS.hg19))]
##################################
#### Part 3. Simulation Process
##################################
# some basic parameters
p <-  dim(genotype)[2] # number of SNPs
Sigma <- cor(genotype) # correlation matrix
n <- dim(genotype)[1] # number of target data (=n_target)
n_signal <- 50 # number of signals
delta <- 0.7 #0.3, 0.5, 0.7 overlap between base and target data (nested)
n_ref <- 4000 # sample size of base data

# generate critical sequences
alpha <- 1/sqrt(log(p)) #proportion est tends to be larger
set.seed(1)
N_simulate <- 1000
zz <- mvrnorm(N_simulate, rep(0, p), diag(1, p, p))
Vn <- rep(NA, N_simulate)
for(i in 1 : N_simulate){
  z_order <- sort(zz[i,], decreasing = T, index.return = T)$x 
  p_order <- 1- pnorm(z_order)
  ind <- 1 : p
  U <- rep(0, p)
  U <- abs(ind / p - p_order) / sqrt(p_order/2)
  Vn[i] <- max(U[2 : floor(p / 2)], 0)
}
c05 <- quantile(Vn, 1 - alpha)  #c05=0.04617 for p=5089

# standardize genotype matrix for regression model
x <- scale(as.matrix(genotype[1:n,]))

# split sample to training and validation (half - half)
set.seed(1)
train.ratio <- 0.5
ind.train <- sample(nrow(genotype), round(train.ratio*dim(genotype)[1]))
ind.test <- setdiff(rows_along(genotype), ind.train)
x1 <- x[ind.train,]
x2 <- x[ind.test,]
#############################################
#### Part 4. Simulation Function for Cluster
#############################################
f <- function(k) {
  library(MASS)
  library(hdi) 
  library(glmnet)
  library(bigsnpr)
  library(data.table)
  library(magrittr)
  library(tidyverse)
  set.seed(2)
  nseed <- ceiling(runif(loopNum)*10000)
  set.seed(nseed[k])
  
  ### Generate summary statistics ###
  locate.base <- sample(1 : p, n_signal, replace=F) # generate the location of active SNPs
  
  ## simulate summary statistics for screening based on reference data
  intensity <- rep(0, p)
  beta.base <- runif(n_signal, 0.05, 0.15) # generate beta for base data
  
  #negative.base <- sample(1 : length(beta.base), floor(length(beta.base)/2), replace=F) # select some signal as negative
  #beta.base[negative.base] <- -1 * beta.base[negative.base]
  
  intensity[locate.base] = sqrt(n_ref)*beta.base # put beta into intensity vector in which noise = 0
  t <- mvrnorm(1, intensity, Sigma) ## ~4.4mins generate base data
  
  ## estimate signal proportion based on summary of reference data
  ## will translate MR functions to p-value version, instead of z-value version
  pi <- MR05(t, c05)
  # get index for sorted t (t is one-sided so abs is not needed)
  tt <- t # for two-sided
  t.ind <- sort(tt, decreasing=T, index.return=T)$ix
  
  ## Simulate training and testing data
  # sample some signal from base data as signal in target data with ratio delta
  target.index <- c(sample(1:n_signal, floor(n_signal*delta), replace=F)) # select index for target signal
  locate.target <- locate.base[target.index] # select target signal location from base signal location list
  beta.target.prep <- runif(n_signal, 0.05, 0.15) # create a list which length = # signal in base
  beta.target <- beta.target.prep[target.index] # select target signal in this list
  
  #negative.target <- sample(1 : length(beta.target), floor(length(beta.target)/2), replace=F) # select some signal as negative
  #beta.target[negative.target] <- -1 * beta.target[negative.target]
  
  beta <- rep(0, p)
  beta[locate.target] <- beta.target
  y <- x %*% beta + rnorm(n, 0, 1)
  y1 <- y[ind.train]
  y2 <- y[ind.test]
  # calculate SST for testing data
  sst <- sum((y2-mean(y2))^2)
  n.test <- length(y2)
  
  ## FNP+Lasso with optimal epsilon level
  epsilon <- 0.02*c(1:40) # a series of FNP levels
  rss.lasso <- rep(NA, length(epsilon)) # for each levels
  lasso.tune <-  matrix(0, length(epsilon), p)
  lasso.final <- matrix(0, 1, p)
  fnp.tune <- fnp.opt.series(t, epsilon, s = max(ceiling(pi*p), 1), side = 2)$j_hat
  for (tune in 1:length(epsilon)){
    ## fnp screening, remove taking "min"? 
    #fnp.tune[tune] <- min(fnp.tune[tune], (dim(x1)[1]-1))
    ## indices of selected variables
    screen.ind <- t.ind[1:fnp.tune[tune]]
    ## Lasso on fnp reduced set and residual sum of square
    cv.lasso <- cv.glmnet(x1[,screen.ind], y1)
    lasso.tune[tune, screen.ind] = as.numeric(coef(cv.lasso, s=cv.lasso$lambda.min))[-1]
    pred.lasso <- x1 %*% lasso.tune[tune,]
    rss.lasso[tune] <- sum((y1-pred.lasso)^2) 
  } 
  i.star <- which.min(rss.lasso) # optimized index
  lasso.final[1,] <- lasso.tune[i.star,]
  fnp.cut <- fnp.tune[i.star]
  pred.valid.lasso <- x2 %*% lasso.final[1,]
  rss.valid.lasso <- sum((y2-pred.valid.lasso)^2)
  ## prediction R2
  R.lasso <- 1-rss.valid.lasso/sst
  ## variable selection
  lasso.ind <- which(lasso.final[1,]!=0) # selected SNPs
  n.lasso <- length(lasso.ind)
  TP.lasso <- sum(locate.target %in% lasso.ind)
  precision.lasso <- TP.lasso/length(lasso.ind)
  recall.lasso <- TP.lasso/floor(n_signal*delta)
  fm.lasso <- sqrt(precision.lasso*recall.lasso)
  f1.lasso <- 2*precision.lasso*recall.lasso/(precision.lasso+recall.lasso)
  aic.lasso <- 2*n.lasso + n.test * log(rss.valid.lasso/n.test)
  
  ########################################################################################################################
  ########################################################################################################################
  # Compare with SIS+Lasso 
  sis.lasso.final = matrix(0, 1, p)
  screen.ind.sis =  t.ind[1:(dim(x1)[1]-1)]
  ## Lasso on sis reduced set
  cv.lasso1 <- cv.glmnet(x1[, screen.ind.sis], y1)
  sis.lasso.final[1, screen.ind.sis] <- as.numeric(coef(cv.lasso1, s=cv.lasso1$lambda.min))[-1]
  pred.valid.sis.lasso <- x2 %*% sis.lasso.final[1, ]
  rss.valid.sis.lasso <- sum((y2-pred.valid.sis.lasso)^2)  
  ## prediction R2
  R.sis.lasso <- 1-rss.valid.sis.lasso/sst
  ## variable selection
  sis.lasso.ind <- which(sis.lasso.final[1,]!=0) # selected SNPs
  n.sis.lasso <- length(sis.lasso.ind)
  TP.sis.lasso <- sum(locate.target %in% sis.lasso.ind)
  precision.sis.lasso <- TP.sis.lasso/length(sis.lasso.ind)
  recall.sis.lasso <- TP.sis.lasso/floor(n_signal*delta)
  fm.sis.lasso <- sqrt(precision.sis.lasso*recall.sis.lasso)
  f1.sis.lasso <- 2*precision.sis.lasso*recall.sis.lasso/(precision.sis.lasso+recall.sis.lasso)
  aic.sis.lasso <- 2*n.sis.lasso + n.test * log(rss.valid.sis.lasso/n.test)
  
  ########################################################################################################################
  ########################################################################################################################
  # Compare with just Lasso 
  cv.lasso2 <- cv.glmnet(x1, y1)
  just.lasso <- as.numeric(coef(cv.lasso2, s=cv.lasso2$lambda.min))[-1]
  pred.valid.just.lasso <- x2 %*% just.lasso
  rss.valid.just.lasso <- sum((y2-pred.valid.just.lasso)^2)
  ## prediction R2
  R.just.lasso <- 1-rss.valid.just.lasso/sst
  ## variable selection
  just.lasso.ind <- which(just.lasso!=0) # selected SNPs
  n.just.lasso <- length(just.lasso.ind)
  TP.just.lasso <- sum(locate.target %in% just.lasso.ind)
  precision.just.lasso <- TP.just.lasso/length(just.lasso.ind)
  recall.just.lasso <- TP.just.lasso/floor(n_signal*delta)
  fm.just.lasso <- sqrt(precision.just.lasso*recall.just.lasso)
  f1.just.lasso <- 2*precision.just.lasso*recall.just.lasso/(precision.just.lasso+recall.just.lasso)
  aic.just.lasso <- 2*n.just.lasso + n.test * log(rss.valid.just.lasso/n.test)
  
  
  ########################################################################################################################
  ########################################################################################################################
  # Compared with lassosum
  POS2 <- snp_asGeneticPos(CHR, POS.hg19, dir = "tmp-data", ncores = 1)
  df_beta <- data.frame(beta = t, beta_se = 1, n_eff = n_ref)
  corr0 <- snp_cor(G, ncores = 1, infos.pos = POS2, size = 3 / 1000)
  corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))
  
  beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = 1)
  params2 <- attr(beta_lassosum2, "grid_param")
  pred_grid2 <- big_prodMat(G, beta_lassosum2)
  
  params2$score <- apply(pred_grid2[ind.train, ], 2, function(x) {
    if (all(is.na(x))) return(NA)
    summary(lm(y1 ~ x))$coef["x", 3]
  })
  
  best_grid_lassosum2 <- params2 %>%
    mutate(id = row_number()) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_lassosum2[, .]
  
  pred_lassosum <- big_prodVec(G, best_grid_lassosum2, ind.row = ind.test)
  dt.lassosum <- data.frame(y2, pred_lassosum)
  lm.fit.lassosum  <- lm(y2~., data = dt.lassosum)
  R.lassosum.prs = summary(lm.fit.lassosum)$r.squared
  
  lassosum.ind = which(best_grid_lassosum2!=0)
  
  ind.lassosum = lassosum.ind
  n.lassosum.prs = length(ind.lassosum)
  TP.lassosum = sum(locate.target %in% ind.lassosum)
  precision.lassosum.prs = TP.lassosum/length(ind.lassosum)
  recall.lassosum.prs = TP.lassosum/floor(n_signal*delta)
  fm.lassosum.prs <- sqrt(precision.lassosum.prs*recall.lassosum.prs)
  f1.lassosum.prs <- 2*precision.lassosum.prs*recall.lassosum.prs/(precision.lassosum.prs+recall.lassosum.prs)
  rss.lassosum.prs <- sum((summary(lm.fit.lassosum)$residuals)^2)
  aic.lassosum.prs = 2*n.lassosum.prs + n.test * log(rss.lassosum.prs/n.test)
  
  ########################################################################################################################
  ########################################################################################################################
  # Compared with C+T
  ## grid clumping
  p.value = 2*(1-pnorm(abs(t)))
  beta <- t
  lpval <- -log10(p.value)
  all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, exclude = which(is.na(lpval)),
                                lpS = lpval, ncores = 1)
  ## run PRS with different settings
  multi_PRS <- snp_grid_PRS(G, all_keep, beta, lpval, ind.row = ind.train,
                            grid.lpS.thr = 0.9999 * seq_log(max(0.1, min(lpval, na.rm = TRUE)), max(lpval[lpval!=Inf]), 50),
                            n_thr_lpS = 50, ncores = 1)
  grid <- attr(all_keep, "grid") %>%
    mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = row_number()) %>%
    unnest(cols = "thr.lp")
  s <- nrow(grid)
  ## find the highest R2
  grid$r2 <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
    single_PRS <- X[, ind]
    dt = data.frame(y=y.train, x=single_PRS)
    lm.fit = lm(y~x, data = dt)
    (coef(summary(lm.fit))[, "t value"][2])^2
  }, ind = 1:s, s = s, y.train = y1,
  a.combine = 'c', block.size = 1, ncores = 1)
  ## find the best coefficients for SNPs
  max_prs <- grid[which.max(grid$r2), ]
  ind.keep <- unlist(purrr::map(all_keep, max_prs$id))
  ind.ct <- ind.keep[lpval[ind.keep] > max_prs$thr.lp] # selected SNPs by the best model
  ## find the corresponding best PRS
  opt.prs <- snp_PRS(G, beta[ind.keep], ind.keep = ind.keep, ind.test = ind.test,
                     lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp)
  dt <- data.frame(y2, opt.prs)
  ## get the prediction R2 for y ~ PRS
  lm.fit.ct <- lm(y2~., data = dt)
  R.ct.prs <- summary(lm.fit.ct)$r.squared
  ## variable selection
  TP.ct <- sum(locate.target %in% ind.ct)
  n.ct.prs <- length(ind.ct)
  precision.ct.prs <- TP.ct/length(ind.ct)
  recall.ct.prs <- TP.ct/floor(n_signal*delta)
  fm.ct.prs <- sqrt(precision.ct.prs*recall.ct.prs)
  f1.ct.prs <- 2*precision.ct.prs*recall.ct.prs/(precision.ct.prs+recall.ct.prs)
  rss.ct.prs <- sum((summary(lm.fit.ct)$residuals)^2)
  aic.ct.prs <- 2*n.ct.prs + n.test * log(rss.ct.prs/n.test)
  
  ########################################################################################################################
  ########################################################################################################################
  # Compared with LDpred 
  POS2 <- snp_asGeneticPos(CHR, POS.hg19, dir = "tmp-data", ncores = 1)
  df_beta <- data.frame(beta = t, beta_se = 1, n_eff = n_ref)
  corr0 <- snp_cor(G, ncores = 1, infos.pos = POS2, size = 3 / 1000)
  corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))
  ## estimate h2
  ldsc <- snp_ldsc2(corr0, df_beta)
  h2_est <- abs(ldsc[["h2"]])
  ## generate grids of parameters
  h2_seq <- h2_est * seq(0.10, 1.10, 0.2)
  p_seq <- signif(seq_log(1e-3, 1, length.out = 20), 2)
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(TRUE))
  ## get beta for each SNPs
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 1)
  ## generate PRS
  pred_grid <- big_prodMat(G, beta_grid)
  pred_grid[ind.train, ][is.na(pred_grid[ind.train, ])] <- 0
  params$score <- big_univLinReg(as_FBM(pred_grid[ind.train, ]), y1)$score
  params %>%
    mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
    arrange(desc(score)) %>%
    mutate_at(c("score", "sparsity"), round, digits = 3) %>%
    slice(1:10)
  best_grid_sp <- params %>%
    mutate(id = row_number()) %>%
    filter(sparse) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]
  ## find the corresponding best PRS
  pred_sp <- big_prodVec(G, best_grid_sp, ind.row = ind.test)
  ## get the prediction R2 for y ~ PRS
  dt.sp <- data.frame(y2, pred_sp)
  lm.fit.sp  <- lm(y2~., data = dt.sp)
  R.ld.prs <- summary(lm.fit.sp)$r.squared
  p.th <- params$p[which.max(params$score)]
  ldpred.ind <- which(best_grid_sp!=0)
  ind.ld <- ldpred.ind # selected SNPs by the best model
  ## variable selection 
  TP.ldpred <- sum(locate.target %in% ind.ld)
  n.ld.prs <- length(ind.ld)
  precision.ld.prs <- TP.ldpred/length(ind.ld)
  recall.ld.prs <- TP.ldpred/floor(n_signal*delta)
  fm.ld.prs <- sqrt(precision.ld.prs*recall.ld.prs)
  f1.ld.prs <- 2*precision.ld.prs*recall.ld.prs/(precision.ld.prs+recall.ld.prs)
  rss.ld.prs <- sum((summary(lm.fit.sp)$residuals)^2)
  aic.ld.prs <- 2*n.ld.prs + n.test * log(rss.ld.prs/n.test)
  
  return(matrix(c(R.just.lasso, R.lasso, R.sis.lasso, R.lassosum.prs, R.ct.prs, R.ld.prs,
                  precision.just.lasso, precision.lasso, precision.sis.lasso, precision.lassosum.prs, precision.ct.prs, precision.ld.prs,
                  recall.just.lasso, recall.lasso,  recall.sis.lasso, recall.lassosum.prs, recall.ct.prs, recall.ld.prs,
                  aic.just.lasso, aic.lasso, aic.sis.lasso, aic.lassosum.prs, aic.ct.prs, aic.ld.prs,
                  n.just.lasso, n.lasso, n.sis.lasso, n.lassosum.prs, n.ct.prs, n.ld.prs), nrow = 30, ncol = 1))
}


#############################################
#### Part 5. Trigger the Cluster
#############################################
loopNum <- 50 # number of simulation loops
n_signal <- 50 # number of signals

cl <- makeCluster(8)
registerDoParallel(cl)
data.out <- matrix(0, 1, 61)
data <- foreach(k = 1:loopNum, .combine = "cbind") %dopar% f(k)
data.out[1, 1] <- 0
for(j in 1:30) {
  data.out[1, (2*(j-1)+2)] <- mean(data[j, ])
  data.out[1, (2*(j-1)+3)] <- sd(data[j, ])
}

write.csv(data.out, "results_simu1_unimputed_n4000_2000_delta07.csv")
write.csv(data, "data_simu1_unimputed_n4000_2000_delta07.csv")
stopImplicitCluster()
stopCluster(cl)

print(data.out)

