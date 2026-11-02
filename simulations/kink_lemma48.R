# Exact check: for AR(1) errors with unit innovation variance, v' Gamma v / omega^2 for unit smooth vectors
rows <- list()
for (phi in c(0.5, 0.8)) for (Tn in c(100, 200, 400, 800, 1600)) {
  tt <- 1:Tn; om2 <- 1/(1-phi)^2
  Gam <- toeplitz(phi^(0:(Tn-1)) / (1-phi^2))
  unit <- function(v) v / sqrt(sum(v^2))
  b1 <- round(0.4*Tn); b2 <- round(0.7*Tn)
  h1 <- pmax(tt-b1,0); h2 <- pmax(tt-b2,0)
  v_hinge <- unit(qr.resid(qr(cbind(1,tt)), h1))                        # hinge on (1,t)
  v_hinge2 <- unit(qr.resid(qr(cbind(1,tt,h2)), h1))                    # hinge with another knot present
  v_step <- unit(qr.resid(qr(cbind(1,tt,h1,h2)), as.numeric(tt > b1)))  # step (derivative direction)
  v_loc <- unit(qr.resid(qr(cbind(1,tt,h2,pmax(tt-(b1+10),0))), h1))   # hinge with a knot only 10 dates away (not smooth at that scale)
  f <- function(v) drop(crossprod(v, Gam %*% v)) / om2 - 1
  rows[[length(rows)+1]] <- data.frame(phi=phi, T=Tn, hinge=f(v_hinge), hinge_2knots=f(v_hinge2), step=f(v_step),
        hinge_close_knot=f(v_loc), T_x_hinge=Tn*f(v_hinge), T_x_step=Tn*f(v_step))
}
r <- do.call(rbind, rows); options(width=170); print(r, row.names=FALSE, digits=3)
write.csv(r, "study_kink_lemma48.csv", row.names=FALSE)
