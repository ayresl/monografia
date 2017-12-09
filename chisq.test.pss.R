chisq.test.pss <- function(effectsize, df, sig.level = 0.05, power, n) {
  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) {
    abs(x - round(x)) < tol
  }
  if (!is.wholenumber(df) | df < 1) stop("input value for df not a positive integer")
  if (effectsize == 0) stop("bad input value for effect size")
  if (df < 1) stop("bad input value for degrees of freedom")
  if (sig.level < 0 | sig.level > 1) stop("input value for sig.level not between 0 and 1")
  w <- abs(effectsize)
  if (missing(n)) {
    if (power < 0 | power > 1) stop("input value for power not between 0 and 1")
    f <- function(x) qchisq(1 - sig.level, df = df) - qchisq(1 - power, df = df, ncp = x)
    lambda <- uniroot(f, c(1 + 1e-10, 1e+04))$root
    n <- lambda / w ^ 2
    header <- "Sample Size for Pearson's Chi-Squared Test"
  }
  if (missing(power)) {
    if (!is.wholenumber(n) | n < 1) stop("input value for n not a positive integer")
    lambda <- n * w ^ 2
    beta <- pchisq(qchisq(1 - sig.level, df = df), df = df, ncp = lambda)
    power <- round(1 - beta, 4)
    header <- "Power of Pearson's Chi-Squared Test"
  }
  cat("\n", header, "\n", "\n")
  cat("  Sample size:", ceiling(n), "(", n, ")", "\n", " Significance level:", sig.level, "\n", " Power:", power, "\n",
      " Effect size:", w, "\n", " Degrees of freedom:", df, "\n", "\n")
}
