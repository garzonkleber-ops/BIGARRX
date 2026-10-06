# Ejecutar con el directorio de trabajo en la carpeta BiGaRxApp.
# Simulación de Monte Carlo: ¿qué región de confianza cubre de verdad la media poblacional?
# Población = parámetros de Iris (medias de grupo y covarianza combinada), n_k = 50.
source("R/motor_canonico.R")
set.seed(2026)
X <- as.matrix(iris[, 1:4]); g <- iris$Species
pop <- bigarx_canonical(X, g, region = "ellipse", ndim = 10)
mu <- apply(X, 2, function(v) tapply(v, g, mean))
Wp <- crossprod(X - mu[as.integer(g), ]) / (150 - 3)
L <- chol(Wp)
xi <- pop$rows[, 1:2]                       # posición poblacional de cada media (escala W)
nsim <- as.integer(Sys.getenv("NSIM", "2000"))
Pm <- array(NA, c(nsim, 3, 2)); cubre <- matrix(NA, nsim, 3 * 4)
q <- qchisq(.95, 2)
for (s in seq_len(nsim)) {
  Xs <- mu[as.integer(g), ] + matrix(rnorm(150 * 4), 150) %*% L
  f <- bigarx_canonical(Xs, g, region = "ellipse", ndim = 2)
  # La media poblacional también se centra con la media global de la muestra
  sg <- sign(colSums(f$rows[, 1:2] * xi)); sg[sg == 0] <- 1
  P <- sweep(f$rows[, 1:2], 2, sg, "*")
  Pm[s, , ] <- P
  # posicion "verdadera" en los ejes de la muestra no interesa: se compara con xi
  dif <- xi - P
  rg <- canonical_regions(f); rg$cov <- rg$cov * sg[1] * sg[2]
  elip <- sapply(1:3, function(i) {
    S <- matrix(c(rg$var1[i], rg$cov[i], rg$cov[i], rg$var2[i]), 2)
    drop(t(dif[i, ]) %*% solve(S, dif[i, ])) <= q
  })
  d2 <- rowSums(dif^2)
  chi <- d2 <= q / 50
  tt <- d2 <= (qt(.975, 147) / sqrt(50))^2
  bon <- d2 <= (qt(1 - .05 / (2 * 3), 147) / sqrt(50))^2
  cubre[s, ] <- c(elip, chi, tt, bon)
}
v_sim <- t(sapply(1:3, function(i) c(var1 = var(Pm[, i, 1]), var2 = var(Pm[, i, 2]),
                                     cov = cov(Pm[, i, 1], Pm[, i, 2]))))
rownames(v_sim) <- levels(g)
rg0 <- canonical_regions(pop)
cat("\nVarianzas de las medias proyectadas: simulación vs ecuaciones 3-4 (como en el artículo)\n")
print(round(cbind(v_sim, krz_var1 = rg0$var1, krz_var2 = rg0$var2, krz_cov = rg0$cov), 4))
cob <- matrix(colMeans(cubre), 3, 4, dimnames = list(levels(g), c("elipse_Krz", "circulo_chi2", "circulo_t", "circulo_Bonf")))
cat("\nCobertura real (nominal 95%) de cada región:\n"); print(round(100 * cob, 1))
invisible(NULL)
