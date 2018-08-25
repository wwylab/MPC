true.beta <- c(6,1) 
true.lambda <- 1.0e-4 
ture.phi <- 3
range.t <- 15.73679

# load posteriors (from "example.R")
post.phi <- scan("posterior/phi.txt")
N <- length(post.phi)
post.beta  <- matrix(scan("posterior/beta.txt"), ncol = N)
post.gamma <- matrix(scan("posterior/gamma.txt"), ncol = N)
M <- nrow(post.gamma)


par(mfrow = c(1,1))
# beta
for (k in 1:2) 
{
  pdf(paste("figures/beta", k, ".pdf", sep = ""), 5, 5)
  plot(post.beta[k,], col = "gray", type = "l", 
       main = "",
       ylab = paste("beta", k, sep =""))
  abline(h = median(post.beta[k,]), col = 4, lwd = 2)
  abline(h = true.beta[k], lty = 2, col = 2)
  dev.off()
}

# gamma
grid <- seq(0, range.t, length = 100)
W <- unlist(lapply(1:M, function(k) pbeta(grid/range.t, k, M - k + 1)))
Ft <- matrix(W, ncol = M)

post.Lambda <- matrix(0, length(id), 100)
id <- 20000 + (1:1000)*10
for (i in 1:length(id)) {
  gamma.i <- post.gamma[,id[i]]
  post.Lambda[i,] <- c(Ft%*%gamma.i)
}

est.gamma <- apply(post.gamma, 1, median)
Lambda <- c(Ft%*%est.gamma)

pdf("figures/Lambda.pdf", 5, 5)
plot(grid, Lambda, type = "n", xlab = "time")
for (i in 1:length(id)) {
  lines(grid, post.Lambda[i,], col = "gray")
}
lines(grid, Lambda, lwd = 2, col = 4, xlab = "time")
abline(0, lambda, col = 2, lty = 2)
dev.off()

# phi
pdf("figures/phi.pdf", 5, 5)
plot(post.phi, col = "gray", type = "l", ylab = "phi")
abline(h = phi, col = 2, lty = 2)
abline(h = median(post.phi), col = 4, lwd = 2)
dev.off()

