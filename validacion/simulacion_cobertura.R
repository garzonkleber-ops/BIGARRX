# =============================================================================
#  Simulacion de cobertura de las regiones de confianza de las medias (revisada)
#  Ejecutar con el directorio de trabajo en la carpeta BiGaRxApp.
#  NSIM (por defecto 2000) y NSIM_BOOT (por defecto 200, con B = 300 replicas)
#  se pueden cambiar con variables de entorno.
# =============================================================================
#  Poblacion: normal multivariante con las medias de grupo y la covarianza
#  combinada W del conjunto de datos (se toman como parametros verdaderos), con
#  los mismos tamanos de grupo.
#
#  Parametro evaluado: xi_i = coordenadas de la media poblacional del grupo i en
#  los ejes canonicos POBLACIONALES, en la escala de W (varianza dentro = 1),
#  centradas en la media global ponderada. Es exactamente x$rows del ajuste a
#  los parametros poblacionales.
#
#  Alineacion: los ejes de cada muestra se orientan con el signo que maximiza la
#  coincidencia con los ejes poblacionales: sg_k = sign(sum_i P_ik * xi_ik).
#  No se rotan; solo se cambia el signo (invertir_ejes()).
#
#  Una region "cubre" si (xi_i - P_i)' Sigma_i^-1 (xi_i - P_i) <= chi2_2 (elipses)
#  o si ||xi_i - P_i|| <= radio (circulos).
# =============================================================================
for (f in c("R/motor_canonico.R", "R/boot_canonico.R", "R/motor_dos_factores.R")) source(f)
nsim <- as.integer(Sys.getenv("NSIM", "2000"))
nsim_boot <- as.integer(Sys.getenv("NSIM_BOOT", "200"))
q <- stats::qchisq(0.95, 2)

cubre_elipse <- function(f, xi, region, theta_est = "uno_mas") {
  f$region <- region
  if (region == "ellipse") f$.krz$theta <- if (theta_est == "ms") f$.krz$d2 / f$.krz$dfH else 1 + f$.krz$d2 / f$.krz$dfH
  S <- .krz_cov(f, 1:2); dif <- xi - f$rows[, 1:2]
  vapply(seq_len(nrow(xi)), function(i) {
    if (min(eigen(S[[i]], symmetric = TRUE, only.values = TRUE)$values) <= 0) return(NA)
    drop(dif[i, ] %*% solve(S[[i]], dif[i, ])) <= q
  }, logical(1))
}

simular <- function(X, g, nombre) {
  g <- droplevels(as.factor(g)); X <- as.matrix(X)
  pop <- bigarx_canonical(X, g, ndim = 10)
  xi <- pop$rows[, 1:2]
  mu <- apply(X, 2, function(v) tapply(v, g, mean))
  K <- nlevels(g); n <- nrow(X); dfE <- n - K; nk <- as.vector(table(g))
  L <- chol(crossprod(X - mu[as.integer(g), ]) / dfE)
  set.seed(2026)
  res <- matrix(NA, nsim, 6 * K); Pm <- array(NA, c(nsim, K, 2))
  for (s in seq_len(nsim)) {
    Xs <- mu[as.integer(g), ] + matrix(stats::rnorm(n * ncol(X)), n) %*% L
    f <- bigarx_canonical(Xs, g, ndim = 10)
    sg <- sign(colSums(f$rows[, 1:2] * xi)); f <- invertir_ejes(f, which(sg < 0))
    Pm[s, , ] <- f$rows[, 1:2]
    d2 <- rowSums((xi - f$rows[, 1:2])^2)
    res[s, ] <- c(cubre_elipse(f, xi, "ellipse", "uno_mas"), cubre_elipse(f, xi, "ellipse", "ms"),
                  cubre_elipse(f, xi, "ellipse_amaro"),
                  d2 <= q / nk, d2 <= (stats::qt(0.975, dfE) / sqrt(nk))^2,
                  d2 <= (stats::qt(1 - 0.05 / (2 * K * (K - 1) / 2), dfE) / sqrt(nk))^2)
  }
  cob <- matrix(100 * colMeans(res, na.rm = TRUE), K,
                dimnames = list(levels(g), c("Krz_W_1+lambda", "Krz_W_lambda_ms", "Amaro_publicada",
                                             "circulo_chi2", "circulo_t", "circulo_Bonf")))
  # Bootstrap (mas lento): elipse normal de las replicas de cada media
  cb_boot <- rep(NA, K); if (nsim_boot > 0) {
    set.seed(7); hit <- matrix(NA, nsim_boot, K)
    for (s in seq_len(nsim_boot)) {
      Xs <- mu[as.integer(g), ] + matrix(stats::rnorm(n * ncol(X)), n) %*% L
      b <- tryCatch(boot_canonical(Xs, g, B = 300, seed = s), error = function(e) NULL)
      if (is.null(b)) next
      sg <- sign(colSums(b$fit$rows[, 1:2] * xi)); b <- invertir_boot(b, which(sg < 0))
      nb <- boot_clouds(b, "media")
      hit[s, ] <- vapply(levels(g), function(k) {
        Pk <- nb[[k]]; S <- stats::cov(Pk); dd <- xi[k, ] - colMeans(Pk)
        drop(dd %*% solve(S, dd)) <= q }, logical(1))
    }
    cb_boot <- 100 * colMeans(hit, na.rm = TRUE)
  }
  cob <- cbind(cob, bootstrap = cb_boot)
  cat(sprintf("\n=== %s: %d simulaciones (%d con bootstrap) - cobertura (%%), nominal 95 %%\n", nombre, nsim, nsim_boot))
  print(round(cob, 1))
  pop$region <- "ellipse"; S <- .krz_cov(pop, 1:2)
  vr <- t(sapply(seq_len(K), function(i) c(var1_real = stats::var(Pm[, i, 1]), var1_formula = S[[i]][1, 1],
                                           var2_real = stats::var(Pm[, i, 2]), var2_formula = S[[i]][2, 2],
                                           cov_real = stats::cov(Pm[, i, 1], Pm[, i, 2]), cov_formula = S[[i]][1, 2])))
  rownames(vr) <- levels(g)
  cat("Varianza real de P_i (simulacion) frente a la formula coherente en escala W (parametros poblacionales):\n")
  print(round(vr, 4))
  invisible(list(cobertura = cob, varianzas = vr))
}

iris_res <- simular(iris[, 1:4], iris$Species, "Iris (n = 150, K = 3)")
cer <- read.csv("datos/ceramica_romana.csv")
cer_res <- simular(cer[, c("al", "fe", "mg", "ca", "na")], cer$sitio, "Ceramica romana (n = 26, K = 4)")
