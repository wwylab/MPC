get.npp <- function(pi, beta, lambda, phi, lambda.c, xi = NULL, G = NULL, S = NULL) 
{
  # censored time
  cc <- rexp(1, lambda.c)
  
  # frailty
  if (is.null(xi)) {
    xi <- rgamma(1, phi, phi)
  } 
  
  # genotype
  if (is.null(G)) {
    p0 <- 1 - (1 - pi)^2
    G <- rbinom(1, 1, p0)  
  } 
  
  # gender
  if (is.null(S)) {
    S <- 2 - rbinom(1, 1, 1/2)
  } else if (S == 0) {
    S <- 2 
  } else if (S == 1) {
    S <- 1
  } else if (S == 2) {
    S <- 2
  } else stop("geneder should be coded as 1 for male and 2 for female")
  
  # number obs to be sampled
  nn <- 5000
  
  x1 <- cbind(G, 0)
  x2 <- cbind(G, 1)
  
  
  h1 <- lambda * xi * exp(x1 %*% beta)
  h2 <- lambda * xi * exp(x2 %*% beta)
  
  w1 <- rexp(1,      h1)
  w2 <- rexp(nn - 1, h2)
  
  w <- c(w1, w2)
  tt <- cumsum(w)
  
  y <- c(tt[tt < cc], cc)
  
  n <- length(y)
  
  if (n == 1) {x <- x1 
  } else x <- rbind(x1, matrix(rep(x2, n - 1), nrow = n - 1, byrow = T))
  obj <- cbind(y, x[1:n,,drop = F], rep(xi, n))
  colnames(obj) <- c("time", "g", "dp", "xi")
  return(obj)
}