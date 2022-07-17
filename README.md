# DCOE-Lasso
X. Jessie Jeng and Yifei Hu

This R package implements a transfer learning procedure for polygenic risk score analysis.  

The package allows users to transfer information from an external GWAS dataset to improve the prediction results of polygenic risk score based on a smaller datasets. Our proposed transfer learning algorithm consists of two main compenents: (1)
conducting dual control of errors (DCOE) marginal screening to extract useful knowledge from
the base data; and (2) conducting Lasso to integrate the extracted knowledge with
the target training data for accurate trans-data prediction.

# Installation
Install DCOELasso from GitHub using

```r{echo = FALSE, message = FALSE}
library(devtools)
install_github(repo = "JessieJeng/DCOE-Lasso")
```
