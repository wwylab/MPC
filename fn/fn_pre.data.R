pre.data <- function(data.obj1, data.obj2, range.t, M, n.family = length(data.obj2)) {
  n1 <- NULL
  n2 <- NULL
  obj <- as.list(1:n.family)
  Ft.list <- as.list(1:n.family)
  ft.list <- as.list(1:n.family)
  
  for (f in 1:n.family) {
    temp1 <- data.obj1[[f]]
    temp2 <- data.obj2[[f]]
    
    # family size (c)
    n1[f] <- nrow(temp1)
    n2[f] <- nrow(temp2)
    
    # id
    id <- temp2[,"ID"]
    # test
    test <- temp2[,"test"]
    
    # survival time
    time <- temp2[,"time"] # raw
    tilde.t <- time/range.t # rescaled
    
    # history
    d <- temp2[,"D"]
    dp <- temp2[,"Dp"]
    
    obj[[f]] <- data.frame(id = id, time = time, tilde.t = tilde.t, test = test, d = d, dp = dp)
    
    # Bernstein Basis
    W <- unlist(lapply(1:M, function(k) pbeta(tilde.t, k, M-k+1)))
    w <- unlist(lapply(1:M, function(k) dbeta(tilde.t, k, M-k+1)))    
    
    Ft.list[[f]] <- matrix(W, ncol = M)
    ft.list[[f]] <- matrix(w, ncol = M)
  }
  result <- list(n1 = n1, n2 = n2, obj = obj, Ft = Ft.list, ft = ft.list)
  return(result)
}
