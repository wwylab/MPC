posterior.update <- function(beta, gamma, phi, xi, llike, n.family, data.obj1, data.obj2, sigma.beta, a.phi, b.phi,
                              delta.beta, delta.gamma, delta.phi, delta.xi, 
                              prelim.data, prelim.proband, frailty, abc, G, nG, M, range.t, allef, mRate) 
{
  p <- length(beta)
  ac.beta <- rep(0, p)
  ac.gamma <- rep(0, M)
  ac.xi <- rep(0, n.family)
  ac.phi <- 0
  
  
  # updata beta
  for (j in 1:p) {
      new.beta <- beta
      # draw new value
      new.beta[j] <- rnorm(1, mean = beta[j], sd = delta.beta[j])    
      
      # prior likelihood
      lprior     <- dnorm(    beta[j], mean = 0, sd = sigma.beta, log = T)
      new.lprior <- dnorm(new.beta[j], mean = 0, sd = sigma.beta, log = T)
      
      # likelihood
      llike.pro <- rep(0, n.family) 
      if (abc) llike.pro <- proband.likelihood(new.beta, gamma, xi, prelim.proband, range.t, allef)
      new.llike <- FamilyLikelihood(new.beta, gamma, xi, n.family, data.obj1, data.obj2, prelim.data, llike.pro, M, range.t, allef, mRate)
      

      # update sample
      r <- exp((new.llike + new.lprior) - (llike + lprior))
      if (r > runif(1)) 
      {
        beta[j] <- new.beta[j]
        ac.beta[j] <- 1
        llike <- new.llike
      }
    }
    
    # update gamma
    for (m in 1:M) {
      new.gamma <- gamma
      # draw new value
      #tmp <- gamma[m] + rnorm(50, 0, min(c(2*gamma[m], delta.gamma[m])))
      #tmp.id <- min(which(tmp > 0))
      #new.gamma[m] <- tmp[tmp.id]
      
      new.gamma[m] <- exp(rnorm(1, mean = log(gamma[m]), sd = delta.gamma[m]))
      
      
      # likelihood
      llike.pro <- rep(0, n.family)
      if (abc) llike.pro <- proband.likelihood(beta, new.gamma, xi, prelim.proband, range.t, allef)
      new.llike <- FamilyLikelihood(beta, new.gamma, xi, n.family, data.obj1, data.obj2, prelim.data, llike.pro, M, range.t, allef, mRate)
      
      # update sample
      r <- exp(new.llike - llike)
      if (r > runif(1)) {
        gamma[m] <- new.gamma[m]
        ac.gamma[m] <- 1
        llike <- new.llike
      }
    }
    
    
    if (frailty) 
    {
      #############
      # update xi #
      #############
      
      tmp.likelihood <- FamilyLikelihood(beta, gamma, xi, n.family, data.obj1, data.obj2, prelim.data, llike.pro, M, range.t, allef, mRate, sum = F)
      
      for (f in 1:n.family)
      {
        n1 <- prelim.data$n1[[f]]
        n2 <- prelim.data$n2[[f]]
        
        data1 <- data.obj1[[f]]
        data2 <- prelim.data$obj[[f]]
        
        Ft <- prelim.data$Ft[[f]]
        ft <- prelim.data$ft[[f]]
        
        # draw xi
        new.xi <- xi
        tmp <- xi[f] + rnorm(50, 0, delta.xi[f])
        tmp.id <- min(which(tmp > 0))
        new.xi[f] <- tmp[tmp.id]
        
        # compute prior
            lprior <- dgamma(    xi[f], shape = phi, rate = phi, log = T)
        new.lprior <- dgamma(new.xi[f], shape = phi, rate = phi, log = T)
        
       # compute likelihood
        f.llike.pro         <- rep(0, n.family)
        if (abc) f.llike.pro <- proband.likelihood(beta, gamma, new.xi, prelim.proband, range.t, allef)
        new.f.llike.unadjusted <- fllike(beta, gamma, new.xi[f], n1, n2, data1, data2, Ft, ft, G, nG, M, range.t, allef, mRate)
        
        f.like <- tmp.likelihood[f]
        new.f.like <- new.f.llike.unadjusted - f.llike.pro[f]
        
      #  # decision 
        r <- exp((new.f.like + new.lprior) - (f.like + lprior))
        if (r > runif(1)) 
        {
          xi[f] <- new.xi[f]
          ac.xi[f] <- 1
          tmp.likelihood[f] <- new.f.like
        }
      }    
      llike <- sum(tmp.likelihood)       
      
      ##############
      # update phi #
      ##############
      
      new.phi <- phi
      
      # draw sample
      tmp <- phi + rnorm(50, 0, delta.phi)
      tmp.id <- min(which(tmp > 0))
      new.phi <- tmp[tmp.id]
      
      # compuate posterior
      post     <- (n.family*    phi + a.phi - 1) * log(    phi) - n.family * lgamma(    phi) - b.phi *     phi +     phi * (sum(log(xi))-sum(xi))
      new.post <- (n.family*new.phi + a.phi - 1) * log(new.phi) - n.family * lgamma(new.phi) - b.phi * new.phi + new.phi * (sum(log(xi))-sum(xi))
      
      ## decision
      r <- exp(new.post - post)
      if (r > runif(1)) 
      {
        phi <- new.phi
        ac.phi <- 1
      }
    } # end frailty
  
  
  obj <- list(llike = llike, 
              beta = beta, gamma = gamma, phi = phi, xi = xi,
              ac.beta = ac.beta, ac.gamma = ac.gamma, ac.phi = ac.phi, ac.xi = ac.xi)

#  print(llike.pro)

  #obj <- list(llike = llike, 
  #            beta = beta, gamma = gamma,
  #            ac.beta = ac.beta, ac.gamma = ac.gamma)
  
  return(obj)
}
