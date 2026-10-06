# =============================================================================
#  Bootstrap del biplot canonico
# =============================================================================
#  Remuestreo estratificado: en cada replica se extraen, con reemplazo, n_k
#  individuos de cada grupo k (se conservan los tamanos de grupo). Se repite el
#  preprocesamiento y el biplot canonico, y se alinea el signo de cada eje con
#  el analisis original (comparando las coordenadas de las medias de grupo).
#
#  Parametros con intervalo de confianza:
#    inercia (%), valor_propio (theta de HE^-1), correlacion_canonica,
#    media (coordenadas de cada media de grupo en cada eje del plano),
#    variable (coordenadas de cada flecha en cada eje del plano),
#    calidad_variable (calidad de cada variable en el plano),
#    correlacion_variable (correlacion de cada variable con cada eje).
#
#  Referencias: Efron y Tibshirani (1993); Nieto, Galindo, Leiva y
#  Vicente-Galindo (2015); Krzanowski y Radley (1989), bootstrap en el analisis
#  de variables canonicas.
# =============================================================================

boot_canonical <- function(X, groups, B = 1000, transform = c("center", "standardize"),
                           axes = c(1, 2), conf = 0.95, ci = c("percentile", "normal"),
                           seed = NULL, na_action = c("fail", "impute"), impute_rank = 2,
                           progress = NULL, singular = c("stop", "ginv")) {
  singular <- match.arg(singular)
  transform <- match.arg(transform)
  ci <- match.arg(ci)
  na_action <- match.arg(na_action)
  if (!is.numeric(B) || B < 50) stop("`B` debe ser al menos 50.", call. = FALSE)
  fit <- bigarx_canonical(X, groups, transform = transform, ndim = max(3, axes),
                          na_action = na_action, impute_rank = impute_rank, singular = singular)
  if (max(axes) > fit$rank) stop("El plano pedido supera el número de ejes canónicos.", call. = FALSE)
  Xc <- if (!is.null(fit$imputation)) fit$imputation$X else .as_data_matrix(X)
  g <- fit$groups
  P0 <- fit$rows[, axes, drop = FALSE]
  est <- .boot_canon_stats(fit, axes, P0)

  if (!is.null(seed)) set.seed(seed)
  idx_g <- split(seq_len(nrow(Xc)), g)
  reps <- matrix(NA_real_, B, length(est), dimnames = list(NULL, names(est)))
  fallos <- 0L
  for (b in seq_len(B)) {
    idx <- unlist(lapply(idx_g, function(ii) ii[sample.int(length(ii), length(ii), replace = TRUE)]),
                  use.names = FALSE)
    s <- tryCatch({
      fb <- bigarx_canonical(Xc[idx, , drop = FALSE], g[idx], transform = transform,
                             ndim = max(3, axes), fix_signs = FALSE, singular = singular)
      if (fb$rank < max(axes)) NULL else .boot_canon_stats(fb, axes, P0)
    }, error = function(e) NULL)
    if (is.null(s)) fallos <- fallos + 1L else reps[b, ] <- s
    if (is.function(progress) && b %% max(1, B %/% 50) == 0) progress(b / B)
  }
  reps <- reps[stats::complete.cases(reps), , drop = FALSE]
  if (nrow(reps) < 20) stop("Muy pocas réplicas válidas; revisa los tamaños de grupo.", call. = FALSE)

  a <- 1 - conf
  mu <- colMeans(reps)
  se <- apply(reps, 2, stats::sd)
  if (ci == "percentile") {
    lo <- apply(reps, 2, stats::quantile, probs = a / 2, names = FALSE)
    hi <- apply(reps, 2, stats::quantile, probs = 1 - a / 2, names = FALSE)
  } else {
    q <- if (nrow(Xc) >= 25) stats::qnorm(1 - a / 2) else stats::qt(1 - a / 2, nrow(Xc) - 1)
    lo <- est - q * se
    hi <- est + q * se
  }
  info <- do.call(rbind, strsplit(names(est), "|", fixed = TRUE))
  tab <- data.frame(tipo = info[, 1], elemento = info[, 2], eje = info[, 3],
                    estimacion = unname(est), media_boot = unname(mu), ee = unname(se),
                    sesgo = unname(mu - est), inferior = unname(lo), superior = unname(hi),
                    row.names = NULL, stringsAsFactors = FALSE)
  structure(list(fit = fit, table = tab, replicates = reps, B = nrow(reps), fallos = fallos,
                 conf = conf, ci = ci, axes = axes), class = "boot_canonical")
}

# Estadisticos de una replica, con el signo de cada eje alineado con P0.
.boot_canon_stats <- function(f, axes, P0) {
  P <- f$rows[, axes, drop = FALSE]
  sg <- sign(colSums(P * P0)); sg[sg == 0] <- 1
  P <- sweep(P, 2, sg, "*")
  Q <- sweep(f$cols[, axes, drop = FALSE], 2, sg, "*")
  S <- sweep(f$structure[, axes, drop = FALSE], 2, sg, "*")
  ejes <- paste0("Eje", axes)
  grp <- rownames(P); vars <- rownames(Q)
  k <- seq_len(min(3, f$rank))
  ax_all <- paste0("Eje", k)
  out <- c(
    stats::setNames(f$inertia$inercia[k], paste0("inercia|", ax_all, "|", ax_all)),
    stats::setNames(f$canonical$valor_propio[k], paste0("valor_propio|", ax_all, "|", ax_all)),
    stats::setNames(f$canonical$correlacion_canonica[k], paste0("correlacion_canonica|", ax_all, "|", ax_all)),
    stats::setNames(as.vector(P), paste0("media|", rep(grp, 2), "|", rep(ejes, each = length(grp)))),
    stats::setNames(as.vector(Q), paste0("variable|", rep(vars, 2), "|", rep(ejes, each = length(vars)))),
    stats::setNames(rowSums(f$quality_cols[, axes, drop = FALSE]),
                    paste0("calidad_variable|", vars, "|", paste(axes, collapse = "-"))),
    stats::setNames(as.vector(S), paste0("correlacion_variable|", rep(vars, 2), "|",
                                         rep(ejes, each = length(vars))))
  )
  out
}

print.boot_canonical <- function(x, ...) {
  cat(sprintf("Bootstrap del biplot canónico: %d réplicas (estratificado por grupo), IC %s al %s%%\n",
              x$B, x$ci, format(100 * x$conf)))
  t <- x$table[x$table$tipo %in% c("inercia", "correlacion_canonica"), ]
  print(t[, c("tipo", "elemento", "estimacion", "ee", "inferior", "superior")], row.names = FALSE, digits = 4)
  invisible(x)
}

# Nube de replicas (lista de matrices B x 2) para medias o para variables.
boot_clouds <- function(x, tipo = c("media", "variable")) {
  tipo <- match.arg(tipo)
  ejes <- paste0("Eje", x$axes)
  els <- unique(x$table$elemento[x$table$tipo == tipo])
  stats::setNames(lapply(els, function(e) {
    cbind(x$replicates[, paste0(tipo, "|", e, "|", ejes[1])],
          x$replicates[, paste0(tipo, "|", e, "|", ejes[2])])
  }), els)
}

# Elipse que contiene la proporcion `conf` de una nube bivariante (normal).
cloud_ellipse <- function(P, conf = 0.95, npts = 100) {
  mu <- colMeans(P)
  e <- eigen(stats::cov(P), symmetric = TRUE)
  t <- seq(0, 2 * pi, length.out = npts)
  r <- sqrt(stats::qchisq(conf, 2))
  pts <- e$vectors %*% (sqrt(pmax(e$values, 0)) * r * rbind(cos(t), sin(t)))
  data.frame(x = mu[1] + pts[1, ], y = mu[2] + pts[2, ])
}
