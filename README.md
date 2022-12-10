# FNC-Lasso
X. Jessie Jeng and Yifei Hu

This R package implements a transfer learning procedure for polygenic risk score analysis.  

The package allows users to transfer information from an external GWAS dataset to improve the prediction results of polygenic risk score based on a smaller datasets. Our proposed transfer learning algorithm consists of two main compenents: (1)
conducting False Negative Control (FNC) marginal screening to extract useful knowledge from
the base data; and (2) conducting Lasso to integrate the extracted knowledge with
the target training data for accurate trans-data prediction.

# Installation
Install DCOELasso from GitHub using

```r{echo = FALSE, message = FALSE}
library(devtools)
install_github(repo = "JessieJeng/FNC-Lasso")
```

# Example Use


## FNC-Lasso model fitting
The following examples demonstrate how to use the package given a matrix of SNPs (X), a vector of phenotype (y) and a vector of summary statistics from base data (z).  Side is 1 if the test statistics is one-sided and 2 for two-sided test statistics. Epsilon is the series of FNP control levels, which works as a hype-parameter in this algorithm.

```r{echo = FALSE, message = FALSE}
side <- 1
epsilon <- 0.02*c(1:20)
model <- dcoelasso.fit(z, y, X, side, epsilon)
```
## FNC-Lasso model predicting
Given a trained model object, one can easily use predict function to predict the results for a new matrix of SNPs.


```r{echo = FALSE, message = FALSE}
pred <- dcoelasso.predict(model, X.test)
```
