set.seed(1)

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


K.fold <- sample(rep_len(1:5, length(ind.train)))
                 
pi = MR05(scale(t[!is.na(t)]), 0.007304347)
                 
 epsilon <- 0.02*c(1:20) 
 rss.lasso <- rep(NA, length(epsilon))
 lasso.tune <-  matrix(0, length(epsilon), p)
 lasso.final <- matrix(0, 1, p)
 R.lasso <- rep(0, length(epsilon))
 n.lasso <- rep(0, length(epsilon))
 fnp.tune <- fnp.opt.series(scale(t[!is.na(t)]), epsilon, s = ceiling(pi*p), side = 2)$j_hat
 
 for (tune in 1:length(epsilon)){
   screen.ind = t.ind[1: fnp.tune[tune]]
   X.fnc <- FBM(n, length(screen.ind), init = G[, screen.ind])
   mod.fnc <- big_spLinReg(X.fnc , y[ind.train], K = 5, ind.train = ind.train, covar.train = covar_from_df(pc1), ind.sets = K.fold)
   
   pred.lasso = predict(mod.fnc, X.fnc, ind.train, covar.row = covar_from_df(pc1))
   rss.lasso[tune] = sum((y1-pred.lasso)^2)/(sum((y1-mean(y1))^2)*(length(ind.train)-summary(mod.fnc)$nb_var-1))
   
   pred.fnc <- predict(mod.fnc, X.fnc, ind.test, covar.row = covar_from_df(pc2))
   
   n.lasso[tune] <- summary(mod.fnc)$nb_var
   R.lasso[tune] <- 1-sum((y2-pred.fnc)^2)/sst
 }
 
 R.fnc.lasso1 <- R.lasso[which.min(rss.lasso)] 
 n.fnc.lasso1 <- n.lasso[which.min(rss.lasso)] 
 
 
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
 
 #pred.lasso = predict(cv.glin, cbind(X.opt[ind.train, opt.var], pc1), type = "response", lambdaType="lambdaHat1Std")
 #pred.fnc <- predict(cv.glin, cbind(X.opt[ind.test, opt.var], pc2), type = "response", lambdaType="lambdaHat1Std")
 
 pred.lasso = predict(cv.glin, cbind(X.opt[ind.train, opt.var], pc1), type = "response", lambdaType="lambdaHat")
 pred.fnc <- predict(cv.glin, cbind(X.opt[ind.test, opt.var], pc2), type = "response", lambdaType="lambdaHat")
 
 
 
 coefs <- coef(cv.glin)
 
 length(coefs$mainEffects$cont)
 dim(coefs$interactions$contcont)
 
 n.glin <- length(coefs$mainEffects$cat) +  length(coefs$mainEffects$cont)
 n.para.glin <- length(coefs$mainEffects$cat) + length(coefs$mainEffects$cont) + length(coefs$interactions$catcont) + length(coefs$interactions$contcont)
 1-sum((y2-pred.fnc)^2)/sst
 aic.glin <- 2* n.glin + n.test * log(sum((y2-pred.fnc)^2)/n.test)
 
 
 n.fnc.lasso1
 R.fnc.lasso1 
 