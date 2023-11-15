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

setwd('/home3/vvenka23/FNC-Lasso/RedoPaperCoLaus/LDL0.8')
tmpfile <- tempfile()
snp_readBed("/home3/vvenka23/FNC-Lasso/RedoPaperCoLaus/dataHg38/hg38/fixFiles3/merged_QC_8.bed",backingfile = tmpfile)
obj.bigSNP <- snp_attach(paste0(tmpfile, ".rds"))

samples.target <- fread('/home3/vvenka23/FNC-Lasso/RedoPaperCoLaus/dataHg38/hg38/fixFiles2/merged_QC_8.fam')[,1:2]
colnames(samples.target) <- c("FID", "IID") 

pc10 <- fread('/home3/DrYao/CoLaus/data/top10PC.txt')
pc10.sub <- merge(samples.target, pc10, by="IID", sort = FALSE)
phe <- fread('/home3/DrYao/CoLaus/data/CoLausPheno.csv')[,c('Global_ID', "SubjectAlias","AGE", "SEX","HT_cm_","HDLCH", "LDLCH")]
sex <- rep(0, dim(phe)[1])
sex[phe$SEX=='M'] <- 1
phe$SEX <- sex

phenotype <- merge(pc10.sub, phe, by.x="IID", by.y="SubjectAlias", sort = FALSE)
phenotype <- phenotype[!is.na(phenotype$LDLCH),]
samples.idx <- samples.target$IID %in% phenotype$IID # matched samples indice


G <- obj.bigSNP$genotypes[samples.idx, ]
code <- rep(NA_real_, 256)
code[1:3] <- c(0, 1, 2)
G <- FBM.code256(dim(G)[1], dim(G)[2], code, init = G)

CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
RSID <- obj.bigSNP$map$marker.ID
dat <- as.data.frame(cbind(CHR,POS))
colnames(dat) <- c("chr","pos")
POS.hg19 <- snp_modifyBuild(info_snp=dat, liftOver = "/home3/vvenka23/FNC-Lasso/liftover/liftOver",from = "hg38",to = "hg19")$pos
remove.idx <- paste(CHR[which(is.na(POS.hg19))],POS[which(is.na(POS.hg19))],sep=":")
##################################################################################################

sumstats <-
  bigreadr::fread2('/home3/vvenka23/FNC-Lasso/RedoPaperCoLaus/LDL0.8/20686565-GCST000759-EFO_0004611.h.tsv.gz')[,c("hm_rsid", "chromosome", "hm_pos", "hm_effect_allele", "hm_other_allele", "p_value", "beta", "standard_error")]

names(sumstats) <- c("rsid", "chr", "pos","a0", "a1", "p", "beta", "se")
sumstats$n_eff <- 95454
sumstats$chr <- as.numeric(sumstats$chr)
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a0", "a1")
sumstats$chr.pos <- paste(sumstats$chr,sumstats$pos,sep=":")
sumstats <- sumstats[!sumstats$chr.pos %in% remove.idx,]
info_snp <- snp_match(sumstats, map, join_by_pos = TRUE)

beta <- rep(NA, ncol(G))
info_beta <-  info_snp$beta
info_beta[is.na(info_beta)] <- 0
beta[info_snp$`_NUM_ID_`] <- info_beta
lpval <- rep(NA, ncol(G))
lpval[info_snp$`_NUM_ID_`] <- -log10((1-pnorm(abs(info_snp$beta/info_snp$se)))*2)
lpval[info_snp$`_NUM_ID_`] <- -log10(info_snp$p)


lassosum.t3 <- Sys.time()
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")),add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS.hg19, dir = "tmp-data", ncores = 16)
# calculate LD
for (chr in 1:22) {
  # Extract SNPs that are included in the chromosome
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  ord <- order(POS2[ind.chr2])
  # Calculate the LD
  corr0 <- snp_cor(
    G,
    ind.col = ind.chr2[ord],
    ncores = 16,
    infos.pos = POS2[ind.chr2[ord]],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0)
    df_beta <- info_snp[ind.chr[ord], ]
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
    df_beta <- rbind(df_beta, info_snp[ind.chr[ord], ])
  }
}

df_beta[is.na(df_beta)] <- 0
df_beta['beta_se'] <- df_beta['se']


beta_lassosum2 <- snp_lassosum2(corr, df_beta, nlambda = 10, ncores = 16)
params2 <- attr(beta_lassosum2, "grid_param")
pred_grid2 <- big_prodMat(G, beta_lassosum2, ind.col = info_snp$`_NUM_ID_`)


lassosum.t4 <- Sys.time()

lassosum.prep.t <- lassosum.t4 - lassosum.t3
lassosum.prep.t
print(lassosum.prep.t)

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
  y <- phenotype$LDLCH
  pc <- phenotype[,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "AGE", "SEX")]
  y1 <- y[ind.train]
  y2 <- y[ind.test]
  sst <- sum((y2-mean(y2))^2)
  pc1 <- pc[ind.train,]
  pc2 <- pc[ind.test,]
  n.test <- length(ind.test) 
  ##################################################################################################
  lassosum.t1 <- Sys.time()
  params2$score <- apply(pred_grid2[ind.train, ], 2, function(x) {
    if (all(is.na(x))) return(NA)
    summary(glm(y1 ~ x + as.matrix(pc1), family = "gaussian"))$coef["x", 3]
  })
  best_grid_lassosum2 <- params2 %>%
    mutate(id = row_number()) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_lassosum2[, .]
  lassosum.t2 <- Sys.time()
  lassosum.t <- lassosum.t2 - lassosum.t1 
  pred_lassosum <- big_prodVec(G, best_grid_lassosum2, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
  dt.lassosum <- data.frame(y2, pred_lassosum, pc2)
  lm.fit.lassosum  <- lm(y2~., data = dt.lassosum)
  R.lassosum.prs <- summary(lm.fit.lassosum)$r.squared  
  n.lassosum.prs <- sum(best_grid_lassosum2!=0)
  aic.lassosum.prs <- 2*n.lassosum.prs + n.test*log((1-R.lassosum.prs)*sst/n.test)
  return(matrix(c(R.lassosum.prs, n.lassosum.prs, aic.lassosum.prs, lassosum.t), nrow = 4, ncol = 1))
}

loopNum <- 50
data <- matrix(0, 4, loopNum)
for (k in 1:loopNum){
  data[,k] <- f(k)
  print(k)
}


data.out <- matrix(0, 1, 8)
for(j in 1:4) {
  data.out[1, (2*(j-1)+1)] <- mean(data[j, ])
  data.out[1, (2*(j-1)+2)] <- sd(data[j, ])
}

write.csv(data, "LDL_lassosum_data.csv")
write.csv(data.out, "LDL_lassosum_sum.csv")
