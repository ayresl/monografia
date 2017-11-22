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
chisq.test.pss(effectsize = 0.1, df = 5, sig.level = 0.01, power = 0.95)
chisq.test.pss(effectsize = 0.1, df = 5, sig.level = 0.01, n = 2000)

sim_samples <- array(NA, c(B, n))
for (i in 1:B) {
  sim_samples[i, ] <- runif(n, -0.3, 0.7)
}
sim_Tplus <- NULL
for (i in 1:B) {
  sim_Tplus[i] <- sum(rank(abs(sim_samples[i, ]))[sim_samples[i, ] > 0])
}
hist(sim_Tplus, breaks = seq(0, n*(n+1)/2, 1), freq = FALSE)
length(sim_Tplus[sim_Tplus >= qsignrank(1-alpha, n)]) / B

sim_samples <- array(NA, c(B, n))
for (i in 1:B) {
  sim_samples[i, ] <- runif(n, -0.7, 0.3)
}
sim_Tplus <- NULL
for (i in 1:B) {
  sim_Tplus[i] <- sum(rank(abs(sim_samples[i, ]))[sim_samples[i, ] > 0])
}
hist(sim_Tplus, breaks = seq(0, n*(n+1)/2, 1), freq = FALSE)
length(sim_Tplus[sim_Tplus >= qsignrank(1-alpha, n)]) / B

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
wilcox.test.pss(p1 = 0.3, p2 = 0.18, p3 = 0.072, sig.level = 0.1, alternative = "less", power = 0.8)
wilcox.test.pss(p1 = 0.3, p2 = 0.18, p3 = 0.072, sig.level = 0.1, alternative = "two.sided", power = 0.8)
wilcox.test.pss(p1 = 0.7, p2 = 0.82, p3 = 0.712, sig.level = 0.1, alternative = "greater", power = 0.8)
wilcox.test.pss(p1 = 0.7, p2 = 0.82, p3 = 0.712, sig.level = 0.1, alternative = "two.sided", power = 0.8)

gpower.wilcoxon.signed.rank.test <- function(d, alpha = 0.05, power = 0.95, tails = 1, parentdist = "Normal") {
  if (alpha < 0 | alpha > 1) stop("input value for alpha not between 0 and 1")
  if (tails == 1) {
    Tails <- "One"
    f <- function(n) ((qt(1 - alpha, n - 1) + qt(power, n - 1)) / d) ^ 2 - n
  }
  if (tails == 2) {
    Tails <- "Two"
    f <- function(n) ((qt(1 - alpha / 2, n - 1) + qt(power, n - 1)) / d) ^ 2 - n
  }
  res0 <- ceiling(uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root)
  header <- "G*Power Imitation"
  if (parentdist == "Normal") res <- ceiling(res0 * pi / 3)
  if (parentdist == "Laplace") res <- ceiling(res0 * 2 / 3)
  if (parentdist == "Logistic") res <- ceiling(res0 * 9 / pi ^ 2)
  if (parentdist == "min ARE") res <- ceiling(res0 * 1 / 0.864)
  cat("\n", header, "\n", "\n")
  cat("t tests - Means: Wilcoxon signed-rank test (one sample case)", "\n", "\n","Options:   A.R.E. method", "\n", "\n", "Analysis:  A priori: Compute required sample size","\n", "Input:  Tail(s) =", Tails, "\n", " Parent distribution =", parentdist, "\n", " Effect size d = ", d, "\n", " alpha err prob = ", alpha, "\n", " Power (1-beta err prob) = ", power, "\n", "Output:  Df = ", res0, "\n", " Total sample size = ", res)
}
gpower.wilcoxon.signed.rank.test(0.3, alpha = 0.04, power = 0.91, tails = 1, parentdist = "Logistic")

pass.one.mean.test <- function(difference, sd, sd.known = FALSE, sig.level = 0.05, power = 0.8, two.sided = FALSE, adjustment = "uniform") {
  if (sd < 0) stop("input value for sd not positive")
  if (sig.level < 0 | sig.level > 1) stop("input value for sig.level not between 0 and 1")
  d <- abs(difference)
  if(sd.known) {
    z.b <- qnorm(power)
    if (!two.sided) z.a <- qnorm(1 - sig.level) else qnorm(1 - sig.level / 2)
    res <- ceiling(((z.a + z.b) * sd / d) ^ 2)
  }
  else {
    if (!two.sided) f <- function(n) ((qt(1 - sig.level, n - 1) + qt(power, n - 1)) * sd / d) ^ 2 - n else f <- function(n) ((qt(1 - sig.level / 2, n - 1) + qt(power, n - 1)) * sd / d) ^ 2 - n
    res <- ceiling(uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root)
  }
  if (!two.sided && !sd.known) header <- "Sample Size for One-Sided Hypothesis Test for One Mean (Unknown Standard Deviation)"
  if (two.sided && !sd.known) header <- "Sample Size for Two-Sided Hypothesis Test for One Mean (Unknown Standard Deviation)"
  if (!two.sided && sd.known) header <- "Sample Size for One-Sided Hypothesis Test for One Mean (Known Standard Deviation)"
  if (two.sided && sd.known) header <- "Sample Size for Two-Sided Hypothesis Test for One Mean (Known Standard Deviation)"
  if (adjustment == "uniform") res <- ceiling(res * 1)
  if (adjustment == "double exponential") res <- ceiling(res * 2 / 3)
  if (adjustment == "logistic") res <- ceiling(res * 9 / pi ^ 2)
  if (adjustment == "normal") res <- ceiling(res * pi / 3)
  cat("\n", header, "\n", "\n")
  cat(" significance level =", sig.level, " power =", power, " sample size =", res, " standard deviation =", sd)
}
pass.one.mean.test(3, 15, sd.known = FALSE, sig.level = 0.04, power = 0.88, two.sided = FALSE, adjustment = "uniform")
pass.one.mean.test(3, 15, sd.known = FALSE, sig.level = 0.04, power = 0.88, two.sided = FALSE, adjustment = "double exponential")
pass.one.mean.test(3, 15, sd.known = FALSE, sig.level = 0.04, power = 0.88, two.sided = FALSE, adjustment = "logistic")
pass.one.mean.test(3, 15, sd.known = FALSE, sig.level = 0.04, power = 0.88, two.sided = FALSE, adjustment = "normal")

wilcox2.test.pss <- function(p1, p2, p3, sig.level = 0.01, alternative = c("two.sided", "less", "greater"), power, n, a = 1) {
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }
  if (p1 < 0 | p1 > 1) stop("input value for p1 not between 0 and 1")
  if (p2 < 0 | p2 > 1) stop("input value for p2 not between 0 and 1")
  if (p3 < 0 | p3 > 1) stop("input value for p3 not between 0 and 1")
  if (sig.level < 0 | sig.level > 1) stop("input value for sig.level not between 0 and 1")
  if (alternative == "two.sided") z.a <- qnorm(1-sig.level/2) else z.a <- qnorm(1-sig.level)
  if (missing(n)) {
    if (power < 0 | power > 1) stop("input value for power not between 0 and 1")
    z.b <- qnorm(1-power)
    if (alternative == "two.sided") {
      if (p1 < 0.5) {
        f <- function(x) (x*(a*x+x+1))/2 - z.a*sqrt((a*x^2*(a*x+x+1))/12) - a*x^2*p1 - (x*(x+1))/2 +
                          z.b*sqrt(a*x^2*(p1*(1-p1)+(x-1)*(p2-p1^2)+(a*x-1)*(p3-p1^2)))
        n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
      }
      if (p1 > 0.5) {
        f <- function(x) (x*(a*x+x+1))/2 + z.a*sqrt((a*x^2*(a*x+x+1))/12) - a*x^2*p1 - (x*(x+1))/2 -
                          z.b*sqrt(a*x^2*(p1*(1-p1)+(x-1)*(p2-p1^2)+(a*x-1)*(p3-p1^2)))
        n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
      }
    }
    if (alternative == "less") {
      f <- function(x) (x*(a*x+x+1))/2 - z.a*sqrt((a*x^2*(a*x+x+1))/12) - a*x^2*p1 - (x*(x+1))/2 +
                        z.b*sqrt(a*x^2*(p1*(1-p1)+(x-1)*(p2-p1^2)+(a*x-1)*(p3-p1^2)))
      n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
    }
    if (alternative == "greater") {
      f <- function(x) (x*(a*x+x+1))/2 + z.a*sqrt((a*x^2*(a*x+x+1))/12) - a*x^2*p1 - (x*(x+1))/2 -
                        z.b*sqrt(a*x^2*(p1*(1-p1)+(x-1)*(p2-p1^2)+(a*x-1)*(p3-p1^2)))
      n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
    }
    header <- "Sample Size for Wilcoxon's Rank-Sum Test"
  }
  if (missing(power)) {
    m <- a * n
    if (!is.wholenumber(n) | n < 1) stop("input value for n not a positive integer")
    if (alternative == "two.sided") {
      if (p1 < 0.5) {
        beta <- 1 - pnorm(((n*(m+n+1))/2 - z.a*sqrt((m*n*(m+n+1))/12) - (m*n*p1+(n*(n+1))/2)) /
                        sqrt(m*n*(p1*(1-p1)+(n-1)*(p2-p1^2)+(m-1)*(p3-p1^2))))
      }
      if (p1 > 0.5) {
        beta <- pnorm(((n*(m+n+1))/2 + z.a*sqrt((m*n*(m+n+1))/12) - (m*n*p1+(n*(n+1))/2)) /
                        sqrt(m*n*(p1*(1-p1)+(n-1)*(p2-p1^2)+(m-1)*(p3-p1^2))))
      }
    }
    if (alternative == "less") {
      beta <- 1 - pnorm(((n*(m+n+1))/2 - z.a*sqrt((m*n*(m+n+1))/12) - (m*n*p1+(n*(n+1))/2)) /
                          sqrt(m*n*(p1*(1-p1)+(n-1)*(p2-p1^2)+(m-1)*(p3-p1^2))))
    }
    if (alternative == "greater") {
      beta <- pnorm(((n*(m+n+1))/2 + z.a*sqrt((m*n*(m+n+1))/12) - (m*n*p1+(n*(n+1))/2)) /
                      sqrt(m*n*(p1*(1-p1)+(n-1)*(p2-p1^2)+(m-1)*(p3-p1^2))))
    }
    power <- trunc((1-beta) * 10^4) / 10^4
    header <- "Power of Wilcoxon's Rank-Sum Test"
  }
  cat("\n", header, "\n", "\n")
  cat("  Sample size (group 1):", ceiling(n), "(", n, ")", "\n", " Sample size (group 2):", ceiling(a*n), "(", n, ")", "\n", " Significance level:", sig.level, "\n", " Power:", power, "\n",
      " p1 =", p1, "\n", " p2 =", p2, "\n", " p3 =", p3, "\n")
}
wilcox2.test.pss(p1 = 0.70, p2 = 0.80, p3 = 0.80, sig.level = 0.05, alternative = "two.sided", n = 54, a = 1)
wilcox2.test.pss(p1 = 0.70, p2 = 0.80, p3 = 0.80, sig.level = 0.05, alternative = "two.sided", power = 0.80, a = 1)
wilcox2.test.pss(p1 = 0.623, p2 = 0.447, p3 = 0.485, sig.level = 0.05, alternative = "greater", power = 0.9, a = 1)
p1 <- 0.4
p2 <- (1 - p1) / (1 - 2 * p1) ^ 2 * (1 - 2 * p1 - p1 * log(1 / p1 - 1))
p3 <- (1 - p1) / (1 - 2 * p1) * (2 * p2 - 1)
wilcox2.test.pss(p1 = p1, p2 = p2, p3 = p3, sig.level = 0.05, two.sided = FALSE, power = 0.9)

expedgeworth <- function(n, alfa) {
  f <- function(x) pnorm(x) + (3 * n ^ 2 + 3 * n - 1) / (10 * n * (n + 1) * (2 * n + 1)) * (x ^ 3 - 3 * x) * dnorm(x) - (1 - alfa)
  t <- uniroot(f, c(0, 5))$root
  c <- n * (n + 1) / 4 + t * sqrt(n * (n + 1) * (2 * n + 1) / 24)
  k <- qsignrank(alfa, n, lower.tail = FALSE)
  c - k
}

for (n in 13:18) print(2 * (1 - psignrank(qsignrank(0.95, n), n)), digits = 22)

P.falciparum <- c(3.6, 2.91, 3, 2.4, 2.04, 2.58, 2.1, 2.9, 2.7, 3.13)
P.vivax <- c(3.71, 2.82, 2.93, 2.34, 2.45, 3.06, 2.27, 3.48, 2.89, 3.28, 1.94, 3.22)
n1 <- length(P.vivax)
n2 <- length(P.falciparum)
wilcox.test(P.falciparum, P.vivax)
estimatep1 <- expand.grid(P.vivax, P.falciparum)
p1hat <- sum(estimatep1$Var2 >= estimatep1$Var1) / (n1 * n2)
objecto <- expand.grid(P.vivax, P.vivax, P.falciparum)
estimatep2 <- objecto[!(objecto$Var1 == objecto$Var2), ]
p2hat <- sum(estimatep2$Var3 >= estimatep2$Var1 & estimatep2$Var3 >= estimatep2$Var2) / (n1 * n2 * (n1 - 1))
objecto <- expand.grid(P.falciparum, P.falciparum, P.vivax)
estimatep3 <- objecto[!(objecto$Var1 == objecto$Var2), ]
p3hat <- sum(estimatep2$Var3 >= estimatep2$Var1 & estimatep2$Var3 >= estimatep2$Var2) / (n1 * n2 * (n2 - 1))
combinedsample <- c(P.falciparum, P.vivax)
hist(combinedsample)
concentracaoplasmatica <- c(21, 30, 63, 82, 124, 125, 139, 142, 149, 152, 177, 182, 195, 216, 227, 229, 231, 256, 257, 259, 267, 263, 270, 289, 291, 298, 346, 351, 369, 371, 372, 404, 405, 406, 412, 419, 420, 421, 478, 496, 497, 498, 520, 530, 541, 592, 564, 580, 623, 648, 666, 688, 797, 802, 853, 910, 1040, 1200)
fit_g  <- fitdist(concentracaoplasmatica, "gamma")
fit_gama <- fitdistr(concentracaoplasmatica, "gamma")
ks.test(concentracaoplasmatica, "pgamma", 2.25, 0.0056)
ks.test(concentracaoplasmatica, "pgamma", shape = 2.25, scale = 180)
x <- seq(0, 1200, length = 10000); fx <- dgamma(x, shape = 2.25, rate = 0.0056)
plot(x, fx, type = "l", lty = 2)
curve(dgamma(x - 100, shape = 2.25, rate = 0.0056), add = TRUE, col = "darkgreen")
sim_samples <- array(NA, c(10000, 3))
for (i in 1:10000) {
  sim_samples[i, ] <- rgamma(3, shape = 9/4, scale = 180)
  sim_samples[i, 3] <- sim_samples[i, 3] + 100
}
p1hat <- sum(sim_samples[, 3] > sim_samples[, 1]) / 10000; p2hat <- sum(sim_samples[, 3] > sim_samples[, 1] & sim_samples[, 3] > sim_samples[, 2]) / 10000; p3hat <- sum(sim_samples[, 3] > sim_samples[, 1] & sim_samples[, 2] + 100 > sim_samples[, 1]) / 10000
cat(p1hat, p2hat, p3hat)
sim_samplesX <- array(NA, c(100000, 92))
for (i in 1:100000) {
  sim_samplesX[i, ] <- rgamma(92, shape = 9/4, scale = 180)
}
sim_samplesY <- array(NA, c(100000, 92))
for (i in 1:100000) {
  sim_samplesY[i, ] <- rgamma(92, shape = 9/4, scale = 180) + 100
}
sim_T <- NULL
for (i in 1:100000) {
  sim_T[i] <- sum(rank(c(sim_samplesX[i, ], sim_samplesY[i, ]))[93:184])
}
hist(sim_T, breaks = seq(92 * (92 + 1) / 2, 92 * 92 + 92 * (92 + 1) / 2, 1), freq = FALSE)
length(sim_T[sim_T >= 92 * (92 + 1) / 2 + qwilcox(0.05, 92, 92, lower.tail = FALSE)]) / 100000

rm(list = ls())
library(readxl)
library(xtable)
vietnamdraftlottery <- read_excel("Desktop/monografia/vietnamdraftlottery.xlsx")
readvietnam <- as.data.frame(vietnamdraftlottery)
vietnam <- readvietnam[, -1]
rownames(vietnam) <- readvietnam[, 1]
vietname <- stack(vietnam)
vietna <- vietname[complete.cases(vietname), ]
row.names(vietna) <- NULL
names(vietna) <- c("chamada", "mes")
ni <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
Ri <- c(6236, 5836, 7000, 6110, 6447, 5872, 5628, 5377, 4719, 5656, 4462, 3768)
Ri / ni - (366 + 1) / 2
(sum(ni * (Ri / ni - (366 + 1) / 2) ^ 2) / mean(ni)) / 12
sum((Ri / ni - (366 + 1) / 2) ^ 2) / 12
(sum(ni * (Ri / ni - (366 + 1) / 2) ^ 2) / 366) / 12
sqrt((sum(ni * (Ri / ni - (366 + 1) / 2) ^ 2) / 366) / 12)
tapply(vietna$chamada, vietna$mes, mean) - (366 + 1) / 2
sqrt(sum(ni * (Ri / ni - (366 + 1) / 2) ^ 2) / 366)
sum((tapply(vietna$chamada, vietna$mes, mean) - (366 + 1) / 2) ^ 2)
sqrt(sum((tapply(vietna$chamada, vietna$mes, mean) - (366 + 1) / 2) ^ 2))
sqrt(sum((tapply(vietna$chamada, vietna$mes, mean) - (366 + 1) / 2) ^ 2) / 12)
kruskal.test(chamada ~ mes, data = vietna)
12 * (366 * (366 + 1)) ^ (-1) * sum(ni * (Ri / ni - (366 + 1) / 2) ^ 2)
366 * 25.69109 / (366 + 1)
pchisq(qchisq(1 - 0.05, df = 11), df = 11, ncp = 25.69109, lower.tail = FALSE)
lambda <- seq(0, 50, length = 10^6)
poder <- pchisq(qchisq(1 - 0.05, df = 11), df = 11, ncp = lambda, lower.tail = FALSE)
plot(lambda, poder, type = "l", xlab = expression(italic(lambda)), ylab = "poder")
lambda <- seq(0, 52.5, length = 5250)
poder <- pchisq(qchisq(1 - 0.01, df = 11), df = 11, ncp = lambda, lower.tail = FALSE)
coordenadas <- data.frame(lambda, poder)
write.table(coordenadas, "/Users/lucas/Desktop/monografia/plots/data/exkw.dat", col.names = FALSE, row.names = FALSE)
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
  cat("  Sample size:", ceiling(n), "(", n, ")", "\n", " Significance level:", sig.level, "\n", " Power:", power, "\n",
      " Location shift vector:", deltas, "\n", " Noncentrality parameter:", round(lambda, 4), "\n", " Number of groups:", k, "\n", "\n")
}