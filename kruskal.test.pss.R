kruskal.test.pss <- function(k, deltas, weights = rep(1 / k, k), int.fx2, sig.level = 0.05, power, n) {
  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) {
    abs(x - round(x)) < tol
  }
  if (!is.wholenumber(k) | k < 1) stop("input value for k not a positive integer")
  if (sum(weights) != 1) stop("weights must add up to 1")
  if (k < 1) stop("bad input value for k")
  if (sig.level < 0 | sig.level > 1) stop("input value for sig.level not between 0 and 1")
  d <- 12 * int.fx2 ^ 2 * sum(weights * (deltas - mean(deltas)) ^ 2)
  if (missing(n)) {
    if (power < 0 | power > 1) stop("input value for power not between 0 and 1")
    f <- function(x) qchisq(1 - sig.level, df = k - 1) - qchisq(1 - power, df = k - 1, ncp = x)
    lambda <- uniroot(f, c(1 + 1e-10, 1e+04))$root
    n <- lambda / d
    header <- "Sample Size for Kruskal-Wallis' Rank Sum Test"
  }
  if (missing(power)) {
    if (!is.wholenumber(n) | n < 1) stop("input value for n not a positive integer")
    lambda <- n * d
    beta <- pchisq(qchisq(1 - sig.level, df = k - 1), df = k - 1, ncp = lambda)
    power <- round(1 - beta, 4)
    header <- "Power of Kruskal-Wallis' Rank Sum Test"
  }
  cat("\n", header, "\n", "\n")
  cat("  Sample size:", sum(ceiling(n * weights)), "(", n, ")", ceiling(n * weights), "\n", " Significance level:", sig.level, "\n", " Power:", power, "\n",
      " Location shift vector:", deltas, "\n", " Noncentrality parameter:", round(lambda, 4), "\n", " Number of groups:", k, "\n", "\n")
}
