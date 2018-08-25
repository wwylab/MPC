data.gen <- function(seed, N, pi, beta, lambda, lambda.c, phi, miss = T)
{
  set.seed(seed)
  count <- 0
  
  temp2 <- as.list(1:N)
  while(count < N)
  {
    temp.obj <- get.npp(pi, beta, lambda, phi, lambda.c) 
    nt <- nrow(temp.obj)
    if (nt > 1) {
      count <- count + 1
      temp2[[count]] <- temp.obj
      cat("Family", count, "ascertained. \n")
    }
  }
  
  
  
  data1 <- data.frame(ID = 1:30, 
                      Gender = c(1,2,1,2,1,2,
                                 1,2,1,2,
                                 1,2,1,2,
                                 2,1,2,1,
                                 1,2,1,2,1,2,2,1,2,1,2,1),
                      FatherID = c(3,5,0,0,0,0,
                                   1,1,1,1,
                                   3,3,5,5,
                                   0,0,0,0,
                                   11,11,11,
                                   16,16,16,
                                   13,13,13,
                                   18,18,18), 
                      MotherID = c(4,6,0,0,0,0,
                                   2,2,2,2,
                                   4,4,6,6,
                                   0,0,0,0,
                                   15,15,15,
                                   12,12,12,
                                   17,17,17,
                                   14,14,14))
  
  data.obj1 <- data.obj2 <- as.list(1:N)
  
  n1 <- nrow(data1)
  
  # genotype
  for (i in 1:N) {
    
    data.obj1[[i]] <- data1
    
    # ID
    ID <- data1$ID
    
    # gender
    S <- data1$Gender
    
    # genotype
    G <- rep(0, n1)
    
    proband <- temp2[[i]]
    G[1] <- unique(proband[,"g"])
    
    if (G[1] == 1) {
      # spouse family
      spouse.index <- c(2, 5,6, 12,14,15,16,17,18, 25,26,27,28,29,30)
      G[spouse.index] <- 0
      
      # parents
      G[3] <- rbinom(1, 1, 1/2)
      G[4] <- 1 - G[3]
      
      # offspring
      G[7:10] <- rbinom(4, 1, 1/2)
      
      # sibling
      G[11:12] <- rbinom(2, 1, 1/2)
      
      # sib-ofspring1
      if (G[11] == 1) G[19:21] <- rbinom(3,1,1/2)
      
      # sib-ofspring2
      if (G[12] == 1) G[22:24] <- rbinom(3,1,1/2)
    }
    
    # generate data
    d <- dp <- proband[,"dp"]
    if (length(d) > 1) d[2] <- 0
    
    data2 <- cbind(ID = rep(ID[1], nrow(proband)), 
                   time = proband[,"time"],
                   test = proband[,"g"],
                   D = d,
                   Dp = dp)
    
    for (k in 2:n1)
    {
      temp <- get.npp(pi, beta, lambda, phi, lambda.c, proband[1,"xi"], G[k], S[k]) 
      
      d <- dp <- temp[,"dp"]
      if (length(d) > 1) d[2] <- 0
      
      Temp <- cbind(ID = rep(ID[k], nrow(temp)), 
                    time = temp[,"time"],
                    test = temp[,"g"],
                    D = d,
                    Dp = dp)
      
      data2 <- rbind(data2, Temp)
    }
    
    rownames(data2) <- NULL
    
    # genearete missing (50%)
    if (miss)
    {
      m.index <- sample(ID[-1], 15)
      na.index <- !is.na(match(data2[,1], m.index))
      data2[na.index,"test"] <- NA
      
      data.obj2[[i]] <- data2  
    } else {
      data.obj2[[i]] <- data2
    }
    
  }
  
  list(data.obj1 = data.obj1, data.obj2 = data.obj2)
}