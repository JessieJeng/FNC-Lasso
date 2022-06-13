######################################################
# DCOELasso
#' Performs DCOELasso
#'
#' This function implements a transfer learning procedure for polygenic risk score analysis.
######################################################

# Signal Proportion Estimation

signal.est <- function(z, N_simulate=1000, alpha){
  # generate series
  p <- length(z)
  zz <- matrix(0, nrow = N_simulate, ncol = p)
  for (i in 1: N_simulate) {
    for (j in 1:p){
      zz[i, j] <- rnorm(1, 0, 1)
    }
  }
  n <- p
  Vn <- rep(NA, N_simulate)
  for(i in 1 : N_simulate){
    z_order = sort(abs(zz[i,]), decreasing = T, index.return = T)$x
    p_order = (1- pnorm(z_order))*2
    ind = 1 : n
    U = rep(0, n)
    U = abs(ind / n - p_order) / sqrt(p_order/2)
    Vn[i] = max(U[2 : floor(n / 2)], 0)
  }
  c05 <- quantile(Vn, 1 - alpha)
  # get estimation
  z_order <- sort(abs(z), decreasing = T, index.return = T)$x
  p_order <- (1- pnorm(z_order))*2
  H = rep(0, n)
  ind = 1 : n
  H = (ind / n - p_order - c05 * sqrt(p_order/2)) / (1 - p_order)
  pi = max(H[2 : floor(n / 2)], 0)
  return(pi*p)
}


# FNP Estimation
## z: test statistics;
## t: threshold;
## side: 1 for one-sided test/ 2 for two-sided test
fnp.est <- function(z, t, side) {
  s <- signal.est(z)
  p <- length(z)
  if (side==1) {rej=sum(z>t)}
  if (side==2) {rej=sum(abs(z)>t)}
  fnp <- 1 - rej/s + side*(p-s)*pnorm(-t, 0, 1)/s
  if (fnp<0) {fnp=0} else { if (fnp>1) {fnp=1} }
  return(fnp)
}



# DCOE Process
dcoe <- function(z, epsilon, side) {
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


