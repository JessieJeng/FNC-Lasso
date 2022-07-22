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
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos


sumstats <- 
  bigreadr::fread2('/Users/yifeihu/CHOL/20686565-GCST000760-EFO_0004574.h.tsv.gz')[,
                                                                                   c("hm_rsid", "chromosome", "hm_pos", "hm_effect_allele", "hm_other_allele", "p_value", "beta", "standard_error")]

names(sumstats) <- c("rsid", "chr", "pos","a0", "a1","p", "beta", "beta_se")
sumstats$n_eff <- 100184

map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map, join_by_pos = FALSE)




ld.t3 <- Sys.time()
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")),add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = 3)
# calculate LD
for (chr in 1:22) {
  # Extract SNPs that are included in the chromosome
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  # Calculate the LD
  corr0 <- snp_cor(
    G,
    ind.col = ind.chr2,
    ncores = 3,
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}


df_beta <- info_snp
df_beta[is.na(df_beta)] <- 0
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]
p_seq <- signif(seq_log(1e-4, 0.99, length.out = 17), 4)
h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
grid.param <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(TRUE))
# Get adjusted beta from grid model
beta_grid <- snp_ldpred2_grid(corr, df_beta, grid.param, ncores = 3)
# Obtain model PRS
pred_grid <- big_prodMat(G, beta_grid, ind.col = info_snp$`_NUM_ID_`)
ld.t4 <- Sys.time()
ld.prep.t <- ld.t4 - ld.t3


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
  
  y <- phenotype$HDLCH
  pc <- phenotype[,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "AGE", "SEX")]
  y1 <- y[ind.train]
  y2 <- y[ind.test]
  sst <- sum((y2-mean(y2))^2)
  
  pc1 <- pc[ind.train,]
  pc2 <- pc[ind.test,]
  
  n.test <- length(ind.test) 
  ##################################################################################################
  ##################################################################################################
  ld.t1 <- Sys.time()
  grid.param$score <- apply(pred_grid[ind.train,], 2, function(x) {
    if (all(is.na(x))) return(NA)
    summary(glm(y1 ~ x + as.matrix(pc1), family = "gaussian"))$coef["x", 3]
  })
  best_beta_grid <- grid.param %>%
    mutate(id = row_number()) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]
  
  ld.t2 <- Sys.time()
  ld.t <- ld.t2 - ld.t1 
  
  pred_sp <- big_prodVec(G, best_beta_grid, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
  dt.sp <- data.frame(y2, pred_sp, pc2)
  lm.fit.sp  <- lm(y2~., data = dt.sp)
  R.ld.prs <- summary(lm.fit.sp)$r.squared  
  n.ld.prs <- sum(best_beta_grid!=0)
  aic.ld.prs <- 2*n.ld.prs + n.test*log((1-R.ld.prs)*sst/n.test)
  
  return(matrix(c(R.ld.prs, n.ld.prs, aic.ld.prs, ld.t), nrow = 4, ncol = 1))
}

loopNum <- 16
data <- matrix(0, 4, loopNum)
for (k in 1:loopNum){
  data[,k] <- f(k)
}


data.out <- matrix(0, 1, 8)
for(j in 1:4) {
  data.out[1, (2*(j-1)+1)] <- mean(data[j, ])
  data.out[1, (2*(j-1)+2)] <- sd(data[j, ])
}



write.csv(data, "CHOL_ld_data.csv")
write.csv(data.out, "CHOL_ld_sum.csv")
