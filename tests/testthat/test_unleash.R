context("UNLEASH")

test_that("unleash works", {
  set.seed(92)
  ## correlated normals
  p <- 101 ## number of genes
  m <- 53 ## number of null genes
  rho <- 0.4
  Sig <- matrix(rho, nrow = p, ncol = p)
  diag(Sig) <- 1
  esig <- eigen(Sig)
  Sig_half <- esig$vectors %*% sqrt(diag(esig$values)) %*% t(esig$vectors)
  stopifnot(all(Sig_half %*% Sig_half - Sig < 10 ^ -12))

  beta <- c(Sig_half %*% stats::rnorm(p))
  beta[1:m] <- 0

  betahat <- beta + stats::rnorm(p)
  sebetahat <- rep(1, length(betahat))
}
)
