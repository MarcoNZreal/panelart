library(panelart)
Tn <- 300; N <- 10; tt <- 1:Tn
a2N <- 2*N + 2*sqrt(2*N*log(Tn)) + 2*log(Tn)
kind <- function(K, kn) {
  near <- sapply(kn, function(b) sum(abs(K-b) <= 30)); far <- sum(sapply(K, function(k) all(abs(k-kn) > 30)))
  if (length(K) == length(kn) && all(near == 1)) "exact" else if (any(near == 0)) "missed" else if (far > 0 && any(near >= 2)) "far+dup" else if (far > 0) "far" else "dup"
}
## (a) two kinks, slope 0.04, 200 runs, same seeds as the diagnosis
res <- list(alpha_mc=character(0), twoN_thin=character(0), fixed_2N=character(0))
for (i in 1:200) {
  set.seed(23000 + i)
  g1 <- 0.04 * sample(c(-1,1), N, TRUE)
  Y <- outer(pmax(tt-100,0), g1) + outer(pmax(tt-200,0), g1) + matrix(rnorm(Tn*N), Tn, N)
  res$alpha_mc  <- c(res$alpha_mc,  kind(art_kink(Y, "mc", nsim=500)$dates, c(100,200)))
  res$twoN_thin <- c(res$twoN_thin, kind(art_kink(Y, "mc", nsim=500, alpha_thin="2N")$dates, c(100,200)))
  res$fixed_2N  <- c(res$fixed_2N,  kind(art_kink(Y, "fixed", alpha=a2N)$dates, c(100,200)))
}
cat("(a) outcome of 200 runs, two kinks, slope 0.04 (exact / far knot / duplicate / far+dup / missed)\n")
tab <- sapply(res, function(x) table(factor(x, levels=c("exact","far","dup","far+dup","missed")))); print(tab)
## (b) weak knot beside a strong one
h1 <- pmax(tt-100,0); h2 <- pmax(tt-200,0); rows <- list()
for (s2 in c(0.008, 0.012)) {
  det <- matrix(NA, 100, 2)
  for (i in 1:100) {
    set.seed(63000 + i); sg <- sample(c(-1,1), N, TRUE)
    Y <- outer(h1, 0.04*sg) + outer(h2, s2*sg) + matrix(rnorm(Tn*N), Tn, N)
    det[i,1] <- any(abs(art_kink(Y, "mc", nsim=500, alpha_thin="2N")$dates - 200) <= 30)
    det[i,2] <- any(abs(art_kink(Y, "fixed", alpha=a2N)$dates - 200) <= 30)
  }
  rows[[length(rows)+1]] <- data.frame(weak_slope=s2, twoN_thinning=mean(det[,1]), fixed_alpha_2N=mean(det[,2]))
}
cat("\n(b) weak knot found (100 runs)\n"); print(do.call(rbind, rows), row.names=FALSE)
write.csv(cbind(outcome=rownames(tab), tab), "study_kink_fixedalpha.csv", row.names=FALSE)
write.csv(do.call(rbind, rows), "study_kink_fixedalpha_weak.csv", row.names=FALSE)
cat(sprintf("\n2N level used: %.1f\n", a2N))
