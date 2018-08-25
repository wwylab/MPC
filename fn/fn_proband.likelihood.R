proband.likelihood <- function(beta, gamma, xi, prelim.proband, range.t, allef) {
  time <- prelim.proband$obj$time
  time[time == 0] <- 1.0e-12
  tilde.t <- prelim.proband$obj$tilde.t
  n <- length(time)
  
  Lambda <- prelim.proband$Ft %*% gamma
  lambda <- (prelim.proband$ft %*% gamma)/range.t
  
  test0 <- rep(0, n)
  test1 <- rep(1, n)
  d <- prelim.proband$obj$D
  
  tmp.vec0 <- cbind(test0, d)
  tmp.vec1 <- cbind(test1, d)
  
  ratio0 <- tmp.vec0 %*% beta
  ratio1 <- tmp.vec1 %*% beta
  
  log.like0 <- log(lambda) + log(xi) + ratio0 - (Lambda * xi * exp(ratio0))
  log.like1 <- log(lambda) + log(xi) + ratio1 - (Lambda * xi * exp(ratio1))
  
  pG <- allef[[1]][2]
  p0 <- (1-pG)^2
  
  tmp <- exp(log.like0 + log(rep(p0, n))) + exp(log.like1 + log(rep(1-p0, n)))  
  
  obj <- log(c(tmp))
  
  obj
}

