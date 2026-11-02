library(panelart)
phi <- 0.7; Tn <- 300; N <- 10; tt <- 1:Tn; s <- 0.06; kn <- c(100,200)
gen <- function(seed) { set.seed(seed); g <- s*sample(c(-1,1),N,TRUE)
  mu <- Reduce(`+`, lapply(kn, function(b) outer(pmax(tt-b,0), g)))
  E <- apply(matrix(rnorm(Tn*N),Tn,N), 2, function(z) as.numeric(stats::filter(z, phi, "recursive", init=rep(0,1)))) 
  mu + E*sqrt(1-phi^2) }
miss <- 0; out <- list()
for (i in 1:40) {
  Y <- gen(77000+i); f <- art_kink(Y, "hac_mc", nsim=500)
  near <- sapply(kn, function(b) sum(abs(f$dates-b) <= 30))
  if (any(near == 0)) { miss <- miss + 1
    h <- f$history; cat(sprintf("run %d: final knots [%s]; passes:", i, paste(f$dates, collapse=",")))
    cat(paste0(" (", h$pass, ": L=", h$L, " om2=", round(h$omega2,1), " knots=[", h$dates, "])"), "\n")
    # gain of true-knot model at the final omega2: would the fixed-penalty search with the TRUE omega2 find both?
    ft <- art_kink(Y, "fixed", alpha=5.67*art_penalty_kink_mc(N,Tn,0.05,500L,6L), sigma2=5.67); 
    cat("   search with the true omega2 = 5.67 finds: [", paste(ft$dates, collapse=","), "]\n") } }
cat("missed-window runs:", miss, "of 40\n")
