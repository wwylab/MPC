MCMC_MPC <- function(data.obj1, data.obj2, n.sample, M, pG, mRate, init, delta, frailty, abc, range.t, check = T) {
   ########################
   # mutation information #
   ########################
   maf.rate <- pG
   allef <- list(c(1- maf.rate, maf.rate))
   nloci <- 1
   G <- c(0,1,1)
   nG <- length(G)

   ###################
   # hyper parameter #
   ###################
   sigma.beta <- sigma.alpha <- 100 # for beta prior
   a.phi <- b.phi <- 1              # for phi prior

   n.family <- length(data.obj1) # of families: 189
  
   # prepare_data 
   prelim.data <- pre.data(data.obj1, data.obj2, range.t, M) 
   prelim.proband <- pre.proband(data.obj2, range.t, M)
  
   # initial value
   beta  <- init$beta
   gamma <- init$gamma
   phi   <- init$phi
   xi    <- init$xi
  
   if (!frailty) {
     xi <- rep(1, n.family)
    phi <- 1
   }
  
  llike.pro <- rep(0, n.family)
  if (abc) llike.pro <- proband.likelihood(beta, gamma, xi, prelim.proband, range.t, allef)
  llike <- FamilyLikelihood(beta, gamma, xi, n.family, data.obj1, data.obj2, prelim.data, llike.pro, M, range.t, allef, mRate, sum = T)
  
  # step   
  delta.beta <- delta$beta
  delta.gamma <- delta$gamma
  delta.phi <- delta$phi
  delta.xi <- delta$xi
  
  # storage
  p <- length(beta)
  post.beta  <- acpt.beta  <- matrix(0, p, n.sample)
  post.gamma <- acpt.gamma <- matrix(0, M, n.sample)
  
  if (frailty)
  {
    post.phi <- acpt.phi <- rep(0, n.sample)
  }
  
  likelihood <- NULL
  
  # Start: MCMC
  set.seed(1)
  for (iter in 1:n.sample) {
    ttt <- Sys.time()
    updated <- posterior.update(beta, gamma, phi, xi, llike, n.family, data.obj1, data.obj2, sigma.beta, a.phi, b.phi,
                                delta.beta, delta.gamma, delta.phi, delta.xi, 
                                prelim.data, prelim.proband, frailty, abc, G, nG, M, range.t, allef, mRate)
    post.beta[,iter]  <- beta  <- updated$beta
    post.gamma[,iter] <- gamma <- updated$gamma
    
    acpt.beta[,iter]  <- updated$ac.beta
    acpt.gamma[,iter] <- updated$ac.gamma
    
    if (frailty)   {
      post.phi[iter]   <- phi <- updated$phi
      acpt.phi[iter]   <- updated$ac.phi
      xi <- updated$xi
    }
    
    likelihood[iter]  <- llike <- updated$llike
    if (check) 
    {
      cat(iter, " takes ", round(Sys.time() - ttt, 3), ", like: ", llike, ".\n", sep = "")
    }
  } # End: MCMC
  
  if (frailty) 
  {
    posterior <- list(beta = post.beta, 
                      gamma = post.gamma, 
                      phi = post.phi)
    
    acpt.ratio <- list(beta  = apply(acpt.beta , 1, mean), 
                       gamma = apply(acpt.gamma, 1, mean),
                       phi   = mean(acpt.phi))
  } else {
    posterior <- list(beta = post.beta, 
                      gamma = post.gamma)
    acpt.ratio <- list(beta  = apply(acpt.beta , 1, mean), 
                       gamma = apply(acpt.gamma, 1, mean))
  } 
    
  obj <- list(posterior = posterior, acpt.ratio = acpt.ratio, likelihood = likelihood, range.t = range.t)
}
