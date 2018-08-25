pre.proband <- function(data.obj2, range.t, M, n.family = length(data.obj2)) {
  proband.data <- NULL
  for (f in 1:n.family) 
  {
    proband.data <- rbind(proband.data, data.obj2[[f]][1,])
  }       
  
  time <- proband.data[,"time"]
  tilde.t <- time/range.t
  
  D <- proband.data[,"D"]
  
  test <- proband.data[,"test"]    
  obj <- data.frame(time = time, tilde.t = tilde.t, D = D, test = test)
  
  # Bernstein Basis
  W <- unlist(lapply(1:M, function(k) pbeta(tilde.t, k, M-k+1)))
  w <- unlist(lapply(1:M, function(k) dbeta(tilde.t, k, M-k+1)))    
  
  Ft <- matrix(W, ncol = M)
  ft <- matrix(w, ncol = M)
  
  result <- list(obj = obj, Ft = Ft, ft = ft)
  return(result)
}