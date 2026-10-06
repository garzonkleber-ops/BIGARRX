# =============================================================================
#  Grafico del biplot canonico 100 % personalizable (ggplot2)
# =============================================================================

paletas_disponibles <- c("Dark 3", "Set 2", "Okabe-Ito", "Viridis", "Zissou 1", "Grises")

colores_paleta <- function(n, paleta = "Dark 3") {
  switch(paleta,
    "Okabe-Ito" = rep_len(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                            "#D55E00", "#CC79A7", "#000000"), n),
    "Grises" = grDevices::gray.colors(n, start = 0.15, end = 0.7),
    "Viridis" = grDevices::hcl.colors(n, "Viridis"),
    "Zissou 1" = grDevices::hcl.colors(n, "Zissou 1"),
    "Set 2" = grDevices::hcl.colors(n, "Set 2"),
    grDevices::hcl.colors(n, "Dark 3"))
}

formas <- c("Círculo" = 16, "Triángulo" = 17, "Cuadrado" = 15, "Rombo" = 18,
            "Círculo vacío" = 1, "Triángulo vacío" = 2, "Cuadrado vacío" = 0, "Cruz" = 4)

temas <- c("Minimalista" = "minimal", "Clásico" = "classic", "Blanco y negro" = "bw",
           "Claro" = "light", "Gris" = "gray", "Oscuro" = "dark", "Vacío" = "void")

# Opciones por defecto: cualquier opcion puede sobrescribirse.
opciones_grafico <- function(cb, ...) {
  grupos <- levels(cb$groups)
  vars <- rownames(cb$cols)
  o <- list(
    axes = c(1, 2),
    grupos = grupos, variables = vars, min_calidad = 0,
    ver_individuos = TRUE, ver_medias = TRUE, ver_regiones = TRUE, ver_flechas = TRUE,
    etiq_medias = TRUE, etiq_variables = TRUE, etiq_individuos = FALSE,
    paleta = "Dark 3", colores = NULL,
    ind_tam = 1.6, ind_alfa = 0.45, ind_forma = 16, ind_forma_por_grupo = FALSE,
    media_tam = 4, media_forma = 15,
    reg_grosor = 0.9, reg_linea = "solid", reg_relleno = FALSE, reg_alfa = 0.12,
    flecha_color = "#333333", flecha_colores = NULL, flecha_grosor = 0.7, flecha_punta = 0.2,
    flecha_escala = "auto", flecha_mult = 1,
    etiq_var_tam = 3.4, etiq_var_color = NULL, etiq_media_tam = 4, etiq_ind_tam = 2.4,
    evitar_solape = TRUE,
    titulo = "Biplot canónico", subtitulo = "", xlab = "", ylab = "", pie = "",
    tema = "minimal", tam_base = 13, fuente = "sans", leyenda = "right", titulo_leyenda = "Grupo",
    lineas_cero = TRUE, cuadricula = TRUE, fondo = "#FFFFFF",
    xlim = c(NA, NA), ylim = c(NA, NA),
    boot_medias = TRUE, boot_variables = TRUE, boot_nubes = FALSE, boot_color = "#C2410C"
  )
  mod <- list(...)
  o[names(mod)] <- mod
  o
}

grafico_canonico <- function(cb, o = opciones_grafico(cb), boot = NULL) {
  ax <- o$axes
  ejes <- paste0("Eje", ax)
  grupos_all <- levels(cb$groups)
  cols_grupo <- o$colores
  if (is.null(cols_grupo) || length(cols_grupo) < length(grupos_all)) {
    cols_grupo <- stats::setNames(colores_paleta(length(grupos_all), o$paleta), grupos_all)
  }
  cols_grupo <- cols_grupo[grupos_all]
  gsel <- intersect(o$grupos, grupos_all)

  # Escala de las flechas
  qv <- rowSums(cb$quality_cols[, ax, drop = FALSE])
  vsel <- intersect(o$variables, rownames(cb$cols))
  vsel <- vsel[qv[vsel] >= o$min_calidad]
  k <- if (identical(o$flecha_escala, "auto")) {
    ext <- max(abs(cb$individuals[, ax]))
    0.9 * ext / max(sqrt(rowSums(cb$cols[, ax, drop = FALSE]^2)))
  } else 1
  k <- k * o$flecha_mult

  med <- data.frame(x = cb$rows[, ax[1]], y = cb$rows[, ax[2]], grupo = grupos_all)
  med <- med[med$grupo %in% gsel, ]
  ind <- data.frame(x = cb$individuals[, ax[1]], y = cb$individuals[, ax[2]],
                    grupo = as.character(cb$groups), id = rownames(cb$individuals))
  ind <- ind[ind$grupo %in% gsel, ]
  var <- data.frame(x = cb$cols[vsel, ax[1]] * k, y = cb$cols[vsel, ax[2]] * k, variable = vsel)

  g <- ggplot2::ggplot()
  if (o$lineas_cero) {
    g <- g + ggplot2::geom_hline(yintercept = 0, linetype = 3, colour = "grey55") +
      ggplot2::geom_vline(xintercept = 0, linetype = 3, colour = "grey55")
  }
  # Individuos
  if (o$ver_individuos && nrow(ind)) {
    if (o$ind_forma_por_grupo) {
      g <- g + ggplot2::geom_point(data = ind, ggplot2::aes(x = .data$x, y = .data$y,
                                                            colour = .data$grupo, shape = .data$grupo),
                                   size = o$ind_tam, alpha = o$ind_alfa)
    } else {
      g <- g + ggplot2::geom_point(data = ind, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$grupo),
                                   size = o$ind_tam, alpha = o$ind_alfa, shape = o$ind_forma)
    }
    if (o$etiq_individuos) {
      g <- g + .texto(ind, "id", o, size = o$etiq_ind_tam, colour_by_group = TRUE)
    }
  }
  # Regiones de confianza
  boot_ok <- !is.null(boot) && identical(as.numeric(boot$axes), as.numeric(ax))
  rg <- NULL
  if (o$ver_regiones && nrow(med)) {
    if (identical(cb$region, "boot")) {
      if (boot_ok) {
        nm <- boot_clouds(boot, "media")
        rg <- do.call(rbind, lapply(names(nm), function(gr) data.frame(cloud_ellipse(nm[[gr]], boot$conf), grupo = gr)))
      }
    } else {
      rg <- .regions_df(cb, ax, grupo = grupos_all)
    }
  }
  if (!is.null(rg) && nrow(rg)) {
    rg <- rg[rg$grupo %in% gsel, ]
    if (o$reg_relleno) {
      g <- g + ggplot2::geom_polygon(data = rg, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$grupo,
                                                             group = .data$grupo),
                                     alpha = o$reg_alfa, colour = NA, show.legend = FALSE)
    }
    g <- g + ggplot2::geom_path(data = rg, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$grupo,
                                                        group = .data$grupo),
                                linewidth = o$reg_grosor, linetype = o$reg_linea)
  }
  # Bootstrap: elipses y nubes
  if (boot_ok) {
    if (o$boot_medias && !identical(cb$region, "boot")) {
      nm <- boot_clouds(boot, "media")
      nm <- nm[names(nm) %in% gsel]
      for (gr in names(nm)) {
        if (o$boot_nubes) {
          g <- g + ggplot2::geom_point(data = data.frame(x = nm[[gr]][, 1], y = nm[[gr]][, 2]),
                                       ggplot2::aes(x = .data$x, y = .data$y), colour = cols_grupo[gr],
                                       alpha = 0.12, size = 0.5)
        }
        g <- g + ggplot2::geom_path(data = cloud_ellipse(nm[[gr]], boot$conf),
                                    ggplot2::aes(x = .data$x, y = .data$y), colour = cols_grupo[gr],
                                    linetype = 2, linewidth = 0.6)
      }
    }
    if (o$boot_variables && o$ver_flechas) {
      nv <- boot_clouds(boot, "variable")
      nv <- nv[names(nv) %in% vsel]
      for (v in names(nv)) {
        P <- nv[[v]] * k
        if (o$boot_nubes) {
          g <- g + ggplot2::geom_point(data = data.frame(x = P[, 1], y = P[, 2]),
                                       ggplot2::aes(x = .data$x, y = .data$y), colour = o$boot_color,
                                       alpha = 0.12, size = 0.5)
        }
        g <- g + ggplot2::geom_path(data = cloud_ellipse(P, boot$conf),
                                    ggplot2::aes(x = .data$x, y = .data$y), colour = o$boot_color,
                                    linewidth = 0.6)
      }
    }
  }
  # Medias
  if (o$ver_medias && nrow(med)) {
    g <- g + ggplot2::geom_point(data = med, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$grupo),
                                 shape = o$media_forma, size = o$media_tam)
    if (o$etiq_medias) {
      g <- g + .texto(med, "grupo", o, size = o$etiq_media_tam, colour_by_group = TRUE, bold = TRUE,
                      nudge = TRUE)
    }
  }
  # Flechas de las variables
  if (o$ver_flechas && nrow(var)) {
    fc <- rep(o$flecha_color, nrow(var))
    if (!is.null(o$flecha_colores)) {
      hit <- var$variable %in% names(o$flecha_colores)
      fc[hit] <- o$flecha_colores[var$variable[hit]]
    }
    var$color <- fc
    g <- g + ggplot2::geom_segment(data = var, ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y),
                                   colour = var$color, linewidth = o$flecha_grosor,
                                   arrow = ggplot2::arrow(length = ggplot2::unit(o$flecha_punta, "cm")))
    if (o$etiq_variables) {
      lab <- var
      lab$x <- lab$x * 1.08; lab$y <- lab$y * 1.08
      lc <- if (is.null(o$etiq_var_color)) lab$color else rep(o$etiq_var_color, nrow(lab))
      g <- g + .texto(lab, "variable", o, size = o$etiq_var_tam, colour = lc, bold = TRUE)
    }
  }

  inr <- cb$inertia[ax, ]
  xl <- if (nzchar(o$xlab)) o$xlab else sprintf("Eje %d (%s%%)", ax[1], formatC(inr$inercia[1], format = "f", digits = 1))
  yl <- if (nzchar(o$ylab)) o$ylab else sprintf("Eje %d (%s%%)", ax[2], formatC(inr$inercia[2], format = "f", digits = 1))
  g <- g +
    ggplot2::scale_colour_manual(values = cols_grupo, name = o$titulo_leyenda, drop = TRUE) +
    ggplot2::scale_fill_manual(values = cols_grupo, guide = "none") +
    ggplot2::coord_equal(xlim = if (all(is.na(o$xlim))) NULL else o$xlim,
                         ylim = if (all(is.na(o$ylim))) NULL else o$ylim, clip = "off") +
    ggplot2::labs(title = if (nzchar(o$titulo)) o$titulo else NULL,
                  subtitle = if (nzchar(o$subtitulo)) o$subtitulo else NULL,
                  caption = if (nzchar(o$pie)) o$pie else NULL, x = xl, y = yl)
  if (o$ind_forma_por_grupo) {
    g <- g + ggplot2::scale_shape_manual(values = rep_len(c(16, 17, 15, 18, 1, 2, 0, 4), length(grupos_all)),
                                         name = o$titulo_leyenda)
  }
  tema_fun <- switch(o$tema, classic = ggplot2::theme_classic, bw = ggplot2::theme_bw,
                     light = ggplot2::theme_light, gray = ggplot2::theme_gray,
                     dark = ggplot2::theme_dark, void = ggplot2::theme_void, ggplot2::theme_minimal)
  g <- g + tema_fun(base_size = o$tam_base, base_family = o$fuente) +
    ggplot2::theme(legend.position = o$leyenda,
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.background = ggplot2::element_rect(fill = o$fondo, colour = NA))
  if (!o$cuadricula) g <- g + ggplot2::theme(panel.grid = ggplot2::element_blank())
  g
}

# Texto con o sin ggrepel.
.texto <- function(df, col, o, size, colour = NULL, colour_by_group = FALSE, bold = FALSE, nudge = FALSE) {
  face <- if (bold) "bold" else "plain"
  usa_repel <- isTRUE(o$evitar_solape) && requireNamespace("ggrepel", quietly = TRUE)
  aes_txt <- if (colour_by_group) {
    ggplot2::aes(x = .data$x, y = .data$y, label = .data[[col]], colour = .data$grupo)
  } else {
    ggplot2::aes(x = .data$x, y = .data$y, label = .data[[col]])
  }
  args <- list(data = df, mapping = aes_txt, size = size, fontface = face, show.legend = FALSE,
               family = o$fuente)
  if (!colour_by_group) args$colour <- colour
  if (usa_repel) {
    args$max.overlaps <- Inf
    args$min.segment.length <- 0.3
    args$box.padding <- 0.3
    do.call(ggrepel::geom_text_repel, args)
  } else {
    if (nudge) args$vjust <- -1.1
    do.call(ggplot2::geom_text, args)
  }
}
