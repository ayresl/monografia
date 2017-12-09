wilcox.test.pss <- function(p1, p2, p3, sig.level = 0.01, alternative = c("two.sided", "less", "greater"), power, n) {
  is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) {
    abs(x - round(x)) < tol
  }
  if (p1 < 0 | p1 > 1) stop("input value for p1 not between 0 and 1")
  if (p2 < 0 | p2 > 1) stop("input value for p2 not between 0 and 1")
  if (p3 < 0 | p3 > 1) stop("input value for p3 not between 0 and 1")
  if (sig.level < 0 | sig.level > 1) stop("input value for sig.level not between 0 and 1")
  if (alternative == "two.sided") z.a <- qnorm(1 - sig.level / 2) else z.a <- qnorm(1 - sig.level)
  if (missing(n)) {
    if (power < 0 | power > 1) stop("input value for power not between 0 and 1")
    z.b <- qnorm(1 - power)
    if (alternative == "two.sided") {
      if (p1 < 0.5) {
        f <- function(x) (x * (x + 1)) / 4 - z.a * sqrt((x * (x + 1) * (2 * x + 1)) / 24) - x * (p1 + (x - 1) / 2 * p2) +
                          z.b * sqrt((x * p1 * (1 - p1) + (x * (x - 1)) / 2 * (2 * (p1 - p2) ^ 2 + 3 * p2 * (1 - p2)) + x * (x - 1) * (x - 2) * (p3 - p2 ^ 2)))
        n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
      }
      if (p1 > 0.5) {
        f <- function(x) (x * (x + 1)) / 4 + z.a * sqrt((x * (x + 1) * (2 * x + 1)) / 24) - x * (p1 + (x - 1) / 2 * p2) -
                          z.b * sqrt((x * p1 * (1 - p1) + (x * (x - 1)) / 2 * (2 * (p1 - p2) ^ 2 + 3 * p2 * (1 - p2)) + x * (x - 1) * (x - 2) * (p3 - p2 ^ 2)))
        n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
      }
    }
    if (alternative == "less") {
      f <- function(x) (x * (x + 1)) / 4 - z.a * sqrt((x * (x + 1) * (2 * x + 1)) / 24) - x * (p1 + (x - 1) / 2 * p2) +
                        z.b * sqrt((x * p1 * (1 - p1) + (x * (x - 1)) / 2 * (2 * (p1 - p2) ^ 2 + 3 * p2 * (1 - p2)) + x * (x - 1) * (x - 2) * (p3 - p2 ^ 2)))
      n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
    }
    if (alternative == "greater") {
      f <- function(x) (x * (x + 1)) / 4 + z.a * sqrt((x * (x + 1) * (2 * x + 1)) / 24) - x * (p1 + (x - 1) / 2 * p2) -
                        z.b * sqrt((x * p1 * (1 - p1) + (x * (x - 1)) / 2 * (2 * (p1 - p2) ^ 2 + 3 * p2 * (1 - p2)) + x * (x - 1) * (x - 2) * (p3 - p2 ^ 2)))
      n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
    }
    header <- "Sample Size for Wilcoxon's Signed Rank Test"
  }
  if (missing(power)) {
    if (!is.wholenumber(n) | n < 1) stop("input value for n not a positive integer")
    if (alternative == "two.sided") {
      if (p1 < 0.5) {
        beta <- 1 - pnorm(((n * (n + 1)) / 4 - z.a * sqrt((n * (n + 1) * (2 * n + 1)) / 24) - n * (p1 + (n - 1) / 2 * p2)) /
                            sqrt((n * p1 * (1 - p1) + (n * (n - 1)) / 2 * (2 * (p1 - p2) ^ 2 + 3 * p2 * (1 - p2)) + n * (n - 1) * (n - 2) * (p3 - p2 ^ 2))))
      }
      if (p1 > 0.5) {
        beta <- pnorm(((n * (n + 1)) / 4 + z.a * sqrt((n * (n + 1) * (2 * n + 1)) / 24) - n * (p1 + (n - 1) / 2 * p2)) /
                        sqrt((n * p1 * (1 - p1) + (n * (n - 1)) / 2 * (2 * (p1 - p2) ^ 2 + 3 * p2 * (1 - p2)) + n * (n - 1) * (n - 2) * (p3 - p2 ^ 2))))
      }
    }
    if (alternative == "less") {
      beta <- 1 - pnorm(((n * (n + 1)) / 4 - z.a * sqrt((n * (n + 1) * (2 * n + 1)) / 24) - n * (p1 + (n - 1) / 2 * p2)) /
                          sqrt((n * p1 * (1 - p1) + (n * (n - 1)) / 2 * (2 * (p1 - p2) ^ 2 + 3 * p2 * (1 - p2)) + n * (n - 1) * (n - 2) * (p3 - p2 ^ 2))))
    }
    if (alternative == "greater") {
      beta <- pnorm(((n * (n + 1)) / 4 + z.a * sqrt((n * (n + 1) * (2 * n + 1)) / 24) - n * (p1 + (n - 1) / 2 * p2)) /
                      sqrt((n * p1 * (1 - p1) + (n * (n - 1)) / 2 * (2 * (p1 - p2) ^ 2 + 3 * p2 * (1 - p2)) + n * (n - 1) * (n - 2) * (p3 - p2 ^ 2))))
    }
    power <- trunc((1 - beta) * 10 ^ 4) / 10 ^ 4
    header <- "Power of Wilcoxon's Signed Rank Test"
  }
  cat("\n", header, "\n", "\n")
  cat("  Sample size:", ceiling(n), "(", n, ")", "\n", " Significance level:", sig.level, "\n", " Power:", power, "\n",
      " p1 =", p1, "\n", " p2 =", p2, "\n", " p3 =", p3, "\n")
}
