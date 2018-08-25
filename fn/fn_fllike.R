fllike <- function(beta, gamma, xii, n1, n2, data1, data2, Ft, ft, G, nG, M, range.t, allef, mRate)
{
  id <- data2$id
  time <- data2$time
  time[time == 0] <- 1.0e-12
  tilde.t <- data2$tilde.t
  test <- data2$test
  d <- data2$d
  dp <- data2$dp
  
  
  ##### Key - likelihood part! ########
  #####################################
  test0 <- rep(0, n2)
  test1 <- rep(1, n2)
  
  xp.test0 <- cbind(test0, dp)
  xp.test1 <- cbind(test1, dp)
  
  xpbeta.test0 <- xp.test0 %*% beta
  xpbeta.test1 <- xp.test1 %*% beta
  
  Lambda <- (Ft %*% gamma)         # cumulative baseline
  lambda <- (ft %*% gamma)/range.t # baseline
  
  unique.id <- unique(id)
  diff.Lambda <- NULL
  for (ii in unique.id) {
    diff.Lambda <- c(diff.Lambda, diff(c(0, Lambda[id == ii])))
  }
  
  
  log.l0 <- log(lambda) + log(xii) + xpbeta.test0
  log.l1 <- log(lambda) + log(xii) + xpbeta.test1
  
  log.L0 <- diff.Lambda * xii * exp(xpbeta.test0)
  log.L1 <- diff.Lambda * xii * exp(xpbeta.test1)
  
  llike0 <- llike1 <- NULL
  for (ii in unique.id) {
    temp.id <- which(id == ii)
    
    temp.llike0 <- sum(log.l0[temp.id]) - log.l0[temp.id[length(temp.id)]] - sum(log.L0[temp.id])
    llike0 <- c(llike0, temp.llike0)
    
    temp.llike1 <- sum(log.l1[temp.id]) - log.l1[temp.id[length(temp.id)]] - sum(log.L1[temp.id])
    llike1 <- c(llike1, temp.llike1)
  }
  
  lik <- cbind(exp(llike0), exp(llike1), exp(llike1))
  
  i.test <- NULL
  for (ii in unique.id) {
    temp.id <- which(id == ii)
    i.test <- c(i.test, as.integer(mean(test[temp.id])))
  }
  id0 <- which(i.test == 0)
  id1 <- which(i.test == 1)
  
  lik[id0,2:3] <- 0
  lik[id1,  1] <- 0
  
  obj <- CalPeelingProbLIK(allef, LIK = lik, ped = data1, counselee.id = 1, 1, mRate)
  obj
}
