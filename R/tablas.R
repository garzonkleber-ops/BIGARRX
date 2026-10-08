# =============================================================================
#  Tablas de resultados con estilo APA (como en jamovi)
#  Autores: Kleber Garzon y Roxana Barrezueta. Licencia GPL (>= 3).
# =============================================================================

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits, big.mark = ""))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "< .001", sub("^0", "", formatC(p, format = "f", digits = 3))))
}

# df: data.frame; nombres de columnas ya en espanol. `formatos`: vector con
# "num", "int", "p" o "txt" por columna (por defecto se infiere).
tabla_apa <- function(df, titulo, nota = NULL, digits = 3, formatos = NULL, rownames = TRUE,
                      nombre_filas = "") {
  if (is.null(formatos)) {
    formatos <- vapply(df, function(v) if (is.numeric(v)) {
      if (all(is.na(v) | v == round(v))) "int" else "num"
    } else "txt", character(1))
  }
  cuerpo <- lapply(seq_len(nrow(df)), function(i) {
    celdas <- lapply(seq_along(df), function(j) {
      v <- df[[j]][i]
      txt <- switch(formatos[j],
        p = fmt_p(v), int = ifelse(is.na(v), "", format(v, big.mark = "")),
        num = fmt_num(v, digits), as.character(v))
      htmltools::tags$td(class = if (formatos[j] == "txt") "txt" else "num", txt)
    })
    if (rownames) celdas <- c(list(htmltools::tags$td(class = "txt rn", rownames(df)[i])), celdas)
    htmltools::tags$tr(celdas)
  })
  cab <- lapply(names(df), function(n) htmltools::tags$th(n))
  if (rownames) cab <- c(list(htmltools::tags$th(nombre_filas)), cab)
  htmltools::tags$div(class = "res-tabla",
    htmltools::tags$div(class = "res-titulo", titulo),
    htmltools::tags$table(htmltools::tags$thead(htmltools::tags$tr(cab)),
                          htmltools::tags$tbody(cuerpo)),
    if (!is.null(nota)) htmltools::tags$div(class = "res-nota", htmltools::HTML(nota)))
}

# Descriptivos por grupo (media y desviacion tipica de cada variable).
descriptivos_grupo <- function(X, g) {
  X <- as.data.frame(X)
  filas <- lapply(levels(g), function(k) {
    sub <- X[g == k, , drop = FALSE]
    data.frame(Grupo = k, Variable = names(X), N = nrow(sub),
               Media = vapply(sub, mean, numeric(1)), DT = vapply(sub, stats::sd, numeric(1)),
               row.names = NULL)
  })
  do.call(rbind, filas)
}
