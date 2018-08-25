dyn.load("fn/cpp/CalPeelingProbLIK.so") 

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("fn")

N <- 50
phi <- 3
frailty = T


# estimation setting
n.sample <- 30000 # of posterior samples
M <- 2 # degrees of Bernstein polynomials


# Set-up
true.beta <- beta <- c(6,1) # regression coefficients (test, D)
lambda <- 1.0e-4 # baseline intensity
lambda.c <- 0.5  # haszard for censoring time
pG <- pi <- 0.001 # maf prevalence

# Generate Families
temp <- data.gen(seed = 4, N, pi, beta, lambda, lambda.c, phi, miss = T)

data.obj1 <- temp$data.obj1 # disease history
data.obj2 <- temp$data.obj2 # pedigree

# rescaling constant for Bernstein Polynomials
range.t <- max(unlist(lapply(data.obj2, function(x) max(x[,2])))) + 1.0e-3

# Initial Values
init  <- list(beta = true.beta - 5, gamma = (1:M) * 0.1, xi = rep(1.0, N), phi = 10)
delta <- list(beta = c(0.7, 1), gamma = rep(1, M), xi = rep(0.1, N), phi = 1)

# MCMC sampling
obj <- MCMC_MPC(data.obj1, data.obj2, n.sample, M, pG, mRate = 0, init, delta, frailty = T, abc = T, range.t)

# posterior
post.beta <- obj$posterior$beta
post.gamma <- obj$posterior$gamma
post.phi <- obj$posterior$phi

# aceptance ratio
print(obj$acpt.ratio)


write(post.beta,  "posterior/beta.txt")
write(post.gamma, "posterior/gamma.txt")
write(post.phi,   "posterior/phi.txt")

