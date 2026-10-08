# =============================================================================
#  MANOVA biplot para arreglos de tratamientos con dos factores
#  Amaro, Vicente-Villardon y Galindo (2004), Interciencia 29(1), 26-32.
# =============================================================================
#  Modelo completo con interaccion: E = SCPC residual (dentro de las celdas),
#  gl_E = n - (numero de celdas). Biplots disponibles:
#    total       : las celdas como grupos de un MANOVA de una via.
#    filas       : medias del factor A; ejes que maximizan la F del efecto A.
#    columnas    : medias del factor B; ejes que maximizan la F del efecto B.
#    interaccion : ejes que maximizan la F de la interaccion A x B; se proyectan
#                  las medias de las celdas (interpretacion por paralelogramos).
#  Las SCPC de los efectos son secuenciales (tipo I: A, B, A x B). En disenos
#  balanceados coinciden con las de cualquier otro tipo.
#  Si el efecto tiene un solo grado de libertad (un factor con dos niveles), no
#  hay plano propio: como en Amaro et al. (2004, figura 5), las medias de ese
#  factor se proyectan sobre el biplot total.
# =============================================================================

.scpc_efectos <- function(Z, A, B) {
  fit <- stats::manova(Z ~ A * B)
  ss <- summary(fit)$SS
  df <- stats::anova(stats::lm(Z[, 1] ~ A * B))$Df
  names(df) <- c("A", "B", "A:B", "Residuals")
  list(H = list(A = ss[[1]], B = ss[[2]], AB = ss[[3]]), E = ss[["Residuals"]],
       df = c(A = df[[1]], B = df[[2]], AB = df[[3]]), dfE = df[[4]])
}

bigarx_manova2 <- function(X, A, B, efecto = c("total", "filas", "columnas", "interaccion"),
                           transform = c("center", "standardize"), conf = 0.95,
                           region = c("chi2", "t", "bonferroni", "ellipse", "ellipse_amaro", "boot"), ndim = 3,
                           fix_signs = TRUE, na_action = c("fail", "impute"), impute_rank = 2,
                           singular = c("stop", "ginv"), nombre_A = "A", nombre_B = "B",
                           theta_est = c("uno_mas", "ms")) {
  theta_est <- match.arg(theta_est)
  efecto <- match.arg(efecto); transform <- match.arg(transform); region <- match.arg(region)
  na_action <- match.arg(na_action); singular <- match.arg(singular)
  X <- .as_data_matrix(X, allow_na = na_action == "impute")
  if (length(A) != nrow(X) || length(B) != nrow(X)) stop("Los factores deben tener una entrada por fila.", call. = FALSE)
  if (anyNA(A) || anyNA(B)) stop("Los factores no pueden tener NA.", call. = FALSE)
  A <- droplevels(as.factor(A)); B <- droplevels(as.factor(B))
  if (nlevels(A) < 2 || nlevels(B) < 2) stop("Cada factor necesita al menos 2 niveles.", call. = FALSE)
  celda <- interaction(A, B, sep = "-", lex.order = FALSE, drop = TRUE)
  tab <- table(A, B)
  if (any(tab == 0)) stop("Hay combinaciones de niveles sin observaciones: el modelo con interacción necesita todas las celdas.", call. = FALSE)
  n <- nrow(X); p <- ncol(X); nc <- nlevels(celda)
  if (n - nc < 1) stop("Hace falta al menos una réplica por celda para estimar el error.", call. = FALSE)
  if (n - nc < p && singular == "stop") {
    stop(sprintf("Hay pocos grados de libertad del error (%d) para %d variables.", n - nc, p), call. = FALSE)
  }
  imputation <- NULL
  if (anyNA(X)) { imputation <- impute_svd(X, rank = impute_rank, transform = transform); X <- imputation$X }
  pp <- .preprocess(X, transform)
  Z <- pp$Z
  colnames(Z) <- colnames(X)
  sc <- .scpc_efectos(Z, A, B)
  E <- sc$E; dfE <- sc$dfE

  medias <- function(f) { nk <- as.vector(table(f)); list(M = rowsum(Z, f) / nk, nk = nk) }
  proyectado <- FALSE
  if (efecto == "total") {
    mc <- medias(celda)
    out <- .canon_core(Z, celda, mc$M, mc$nk, E, dfE, NULL, nc - 1, conf, region, ndim, fix_signs, singular, colnames(X), theta_est = theta_est)
  } else if (efecto == "interaccion") {
    if (sc$df[["AB"]] < 2) stop("La interacción tiene un solo grado de libertad: no hay plano de interacción.", call. = FALSE)
    mc <- medias(celda)
    out <- .canon_core(Z, celda, mc$M, mc$nk, E, dfE, sc$H$AB, sc$df[["AB"]], conf,
                       if (region == "boot") "chi2" else region, ndim, fix_signs, singular, colnames(X), theta_est = theta_est)
  } else {
    f <- if (efecto == "filas") A else B
    key <- if (efecto == "filas") "A" else "B"
    mf <- medias(f)
    if (sc$df[[key]] >= 2) {
      out <- .canon_core(Z, f, mf$M, mf$nk, E, dfE, sc$H[[key]], sc$df[[key]], conf,
                         if (region == "boot") "chi2" else region, ndim, fix_signs, singular, colnames(X), theta_est = theta_est)
    } else {
      # Un solo grado de libertad: proyectar las medias del factor sobre el biplot total
      mc <- medias(celda)
      out <- .canon_core(Z, celda, mc$M, mc$nk, E, dfE, NULL, nc - 1, conf, "chi2", ndim, fix_signs,
                         singular, colnames(X), theta_est = theta_est)
      P <- mf$M %*% out$.Wmh %*% out$.V
      colnames(P) <- paste0("Eje", seq_len(ncol(P)))
      out$rows <- P[, seq_len(out$ndim), drop = FALSE]
      out$quality_rows <- .quality(P)
      out$groups <- f
      out$n_groups <- stats::setNames(mf$nk, rownames(mf$M))
      reg <- if (region %in% c("ellipse", "ellipse_amaro", "boot")) "chi2" else region
      out$region <- reg
      Kf <- nlevels(f); alpha <- 1 - conf
      out$radius <- stats::setNames(switch(reg,
        chi2 = sqrt(stats::qchisq(conf, 2) / mf$nk),
        t = stats::qt(1 - alpha / 2, dfE) / sqrt(mf$nk),
        bonferroni = stats::qt(1 - alpha / (2 * max(1, Kf * (Kf - 1) / 2)), dfE) / sqrt(mf$nk)), rownames(mf$M))
      out$.M <- mf$M
      proyectado <- TRUE
    }
  }
  out$preprocessing <- pp$info
  out$imputation <- imputation
  out$diseno <- list(tipo = "dos_factores", efecto = efecto, proyectado = proyectado,
                     A = A, B = B, celda = celda, nombre_A = nombre_A, nombre_B = nombre_B,
                     pruebas = manova2_tests(sc, p_eff = out$.krz$p_eff, nombre_A, nombre_B),
                     anova = anova2_univariada(Z, A, B, nombre_A, nombre_B))
  out$.anova_uni <- NULL
  out
}

# Contrastes multivariantes por efecto (las mismas aproximaciones F que summary.manova)
.pruebas_theta <- function(th, p, q, dfr) {
  s <- min(p, q); m <- (abs(p - q) - 1) / 2; nn <- (dfr - p - 1) / 2
  V <- sum(th / (1 + th))
  pil <- c(V, (2 * nn + s + 1) / (2 * m + s + 1) * V / (s - V), s * (2 * m + s + 1), s * (2 * nn + s + 1))
  W <- prod(1 / (1 + th))
  rr <- dfr - (p - q + 1) / 2; u <- (p * q - 2) / 4
  tt <- if (p^2 + q^2 - 5 > 0) sqrt((p^2 * q^2 - 4) / (p^2 + q^2 - 5)) else 1
  wil <- c(W, ((W^(-1 / tt) - 1) * (rr * tt - 2 * u)) / (p * q), p * q, rr * tt - 2 * u)
  Tt <- sum(th)
  hl <- c(Tt, 2 * (s * nn + 1) * Tt / (s^2 * (2 * m + s + 1)), s * (2 * m + s + 1), 2 * (s * nn + 1))
  d1 <- max(p, q); d2 <- dfr - d1 + q
  roy <- c(max(th), d2 / d1 * max(th), d1, d2)
  tab <- rbind(Pillai = pil, Wilks = wil, `Hotelling-Lawley` = hl, Roy = roy)
  bart <- -(dfr - (p - q + 1) / 2) * log(W)
  data.frame(contraste = rownames(tab), estadistico = tab[, 1], F_aprox = tab[, 2], gl1 = tab[, 3],
             gl2 = tab[, 4], p_valor = stats::pf(tab[, 2], tab[, 3], tab[, 4], lower.tail = FALSE),
             chi2_bartlett = c(NA, bart, NA, NA), gl_chi2 = c(NA, p * q, NA, NA),
             p_chi2 = c(NA, stats::pchisq(bart, p * q, lower.tail = FALSE), NA, NA),
             row.names = NULL)
}

manova2_tests <- function(sc, p_eff, nombre_A = "A", nombre_B = "B") {
  Einv <- MASS_ginv(sc$E)
  ef <- c(A = nombre_A, B = nombre_B, AB = paste(nombre_A, "×", nombre_B))
  do.call(rbind, lapply(names(ef), function(k) {
    th <- Re(eigen(Einv %*% sc$H[[k]], only.values = TRUE)$values)
    th <- sort(pmax(th, 0), decreasing = TRUE)[seq_len(min(p_eff, sc$df[[k]]))]
    data.frame(efecto = ef[[k]], .pruebas_theta(th, p_eff, sc$df[[k]], sc$dfE))
  }))
}

MASS_ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A); pos <- s$d > max(tol * s$d[1], 0)
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

anova2_univariada <- function(Z, A, B, nombre_A = "A", nombre_B = "B") {
  ef <- c(nombre_A, nombre_B, paste(nombre_A, "×", nombre_B))
  do.call(rbind, lapply(colnames(Z), function(v) {
    a <- stats::anova(stats::lm(Z[, v] ~ A * B))
    data.frame(variable = v, efecto = ef, SC = a$`Sum Sq`[1:3], gl1 = a$Df[1:3],
               F = a$`F value`[1:3], gl2 = a$Df[4], p_valor = a$`Pr(>F)`[1:3], row.names = NULL)
  }))
}

# -----------------------------------------------------------------------------
# Comparaciones por pares: T^2 de Hotelling con la W combinada del modelo
# -----------------------------------------------------------------------------
hotelling_pares <- function(x, ajuste = c("bonferroni", "holm", "none")) {
  ajuste <- match.arg(ajuste)
  M <- x$.M; nk <- x$n_groups; dfE <- x$.krz$dfE; p <- x$.krz$p_eff
  Mw <- M %*% x$.Wmh                         # d' W^+ d = ||d W^-1/2||^2
  pares <- utils::combn(rownames(M), 2)
  res <- do.call(rbind, lapply(seq_len(ncol(pares)), function(k) {
    i <- pares[1, k]; j <- pares[2, k]
    D2 <- sum((Mw[i, ] - Mw[j, ])^2)
    T2 <- nk[[i]] * nk[[j]] / (nk[[i]] + nk[[j]]) * D2
    gl2 <- dfE - p + 1
    Fv <- gl2 / (dfE * p) * T2
    data.frame(grupo1 = i, grupo2 = j, D2_Mahalanobis = D2, T2 = T2, F = Fv, gl1 = p, gl2 = gl2,
               p_valor = stats::pf(Fv, p, gl2, lower.tail = FALSE))
  }))
  res$p_ajustado <- stats::p.adjust(res$p_valor, method = ajuste)
  res
}

# -----------------------------------------------------------------------------
# Datos composicionales (proporciones que suman una constante)
# -----------------------------------------------------------------------------
es_composicional <- function(X, tol = 1e-3) {
  X <- as.matrix(X)
  if (any(X < 0, na.rm = TRUE) || ncol(X) < 3) return(FALSE)
  s <- rowSums(X)
  all(is.finite(s)) && stats::sd(s) <= tol * max(1, mean(s))
}

# Reemplazo multiplicativo de ceros: delta = 0.65 * minimo positivo de la columna
.reemplazar_ceros <- function(X) {
  X <- as.matrix(X)
  tot <- rowSums(X)
  for (j in seq_len(ncol(X))) {
    z <- X[, j] == 0 & !is.na(X[, j])
    if (any(z)) {
      pos <- X[X[, j] > 0 & !is.na(X[, j]), j]
      X[z, j] <- 0.65 * if (length(pos)) min(pos) else 1e-6
    }
  }
  X / rowSums(X) * tot
}

transformar_composicional <- function(X, metodo = c("ninguna", "clr", "alr"), referencia = NULL) {
  metodo <- match.arg(metodo)
  if (metodo == "ninguna") return(X)
  nms <- colnames(X)
  Xp <- .reemplazar_ceros(X)
  L <- log(Xp)
  if (metodo == "clr") {
    out <- L - rowMeans(L)
    colnames(out) <- nms
  } else {
    if (is.null(referencia) || !referencia %in% nms) referencia <- nms[length(nms)]
    out <- L[, nms != referencia, drop = FALSE] - L[, referencia]
    colnames(out) <- nms[nms != referencia]
  }
  as.data.frame(out, check.names = FALSE)
}

# -----------------------------------------------------------------------------
# Inversion de ejes en el bootstrap (para que coincida con el grafico)
# -----------------------------------------------------------------------------
invertir_boot <- function(b, ejes) {
  if (is.null(b) || !length(ejes)) return(b)
  for (k in ejes) {
    pat <- paste0("^(media|variable|correlacion_variable)\\|.*\\|Eje", k, "$")
    cols <- grep(pat, colnames(b$replicates))
    b$replicates[, cols] <- -b$replicates[, cols]
    rows <- b$table$tipo %in% c("media", "variable", "correlacion_variable") & b$table$eje == paste0("Eje", k)
    t <- b$table
    lo <- -t$superior[rows]; hi <- -t$inferior[rows]
    t$estimacion[rows] <- -t$estimacion[rows]; t$media_boot[rows] <- -t$media_boot[rows]
    t$sesgo[rows] <- -t$sesgo[rows]; t$inferior[rows] <- lo; t$superior[rows] <- hi
    b$table <- t
  }
  b
}

# Regiones bootstrap: resumen numerico de las elipses de las medias
regiones_boot <- function(b, conf = b$conf) {
  nm <- boot_clouds(b, "media")
  q <- stats::qchisq(conf, 2)
  do.call(rbind, lapply(names(nm), function(gr) {
    P <- nm[[gr]]; S <- stats::cov(P); e <- eigen(S, symmetric = TRUE)
    data.frame(grupo = gr, centro1 = mean(P[, 1]), centro2 = mean(P[, 2]), var1 = S[1, 1], var2 = S[2, 2],
               cov = S[1, 2], semieje_mayor = sqrt(e$values[1] * q), semieje_menor = sqrt(max(e$values[2], 0) * q),
               angulo = 0.5 * atan(2 * S[1, 2] / (S[1, 1] - S[2, 2])) * 180 / pi,
               area = pi * sqrt(max(det(S), 0)) * q)
  }))
}
