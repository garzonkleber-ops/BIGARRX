# =============================================================================
#  Verificacion de las regiones elipticas (Krzanowski, 1989; Amaro et al., 2008)
#  Responde a la revision del 6 de octubre de 2026. Ejecutar con el directorio
#  de trabajo en la carpeta BiGaRxApp:  source("validacion/verificar_elipses.R")
# =============================================================================
for (f in c("R/motor_canonico.R", "R/boot_canonico.R", "R/motor_dos_factores.R")) source(f)
X <- iris[, 1:4]; g <- iris$Species          # 150 individuos, 3 grupos, 4 variables; centrado
r4 <- function(df) { df[] <- lapply(df, function(v) if (is.numeric(v)) round(v, 4) else v); df }

cat("\nCORRESPONDENCIA DE SIMBOLOS (una via: gl_E = n - K = 147, gl_H = K - 1 = 2)\n")
cb <- bigarx_canonical(X, g, region = "ellipse", ndim = 10)
k <- cb$.krz
print(r4(data.frame(
  simbolo = c("d_j", "d_j^2", "theta_HE = d^2/gl_E", "sqrt(theta_HE)", "lambda_j = d^2",
              "theta_j Krzanowski = 1 + lambda/(K-1)", "theta_j (cuadrado medio) = lambda/(K-1)",
              "theta usado en la version publicada"),
  significado = c("valor singular de N^1/2 M W^-1/2 (escala W)",
                  "autovalor de W^-1/2 SCH W^-1/2 = autovalor de SCH W^-1",
                  "autovalor de H E^-1 (MANOVA); 'Valor propio' en la app",
                  "valor singular en escala E; lo que el articulo llama 'valor propio' (5,67)",
                  "estimador del parametro de no centralidad (Anderson, 1984)",
                  "estimador de theta (articulo, p. 234): opcion 'uno_mas' (por defecto)",
                  "autovalor de (SCH/gl_H) W^-1: opcion 'ms'",
                  "1 + sqrt(theta_HE)/(K-1): mezcla valor singular y escala E"),
  eje1 = c(sqrt(k$d2[1]), k$d2[1], k$d2[1] / k$dfE, sqrt(k$d2[1] / k$dfE), k$d2[1],
           1 + k$d2[1] / k$dfH, k$d2[1] / k$dfH, k$theta_amaro[1]),
  eje2 = c(sqrt(k$d2[2]), k$d2[2], k$d2[2] / k$dfE, sqrt(k$d2[2] / k$dfE), k$d2[2],
           1 + k$d2[2] / k$dfH, k$d2[2] / k$dfH, k$theta_amaro[2]))))
cat("xi_ij = P_ij en la escala de W (la de x$rows, la del termino 1/n_i y la de los individuos).\n")
cat("gamma = (K-1)/(n-K) =", round(k$gamma, 5), "\n")

cat("\nTABLA DE VERIFICACION - IRIS, region = 'ellipse' (coherente), nivel 95 %\n")
cat("Escala: centroides, varianzas y covarianzas en la escala de W (varianza dentro de los grupos = 1)\n")
rg <- canonical_regions(cb)
print(r4(rg[, c("grupo", "centro1", "centro2", "var1", "var2", "cov", "semieje_mayor", "semieje_menor",
                "direccion_eje_mayor", "definida_positiva", "error_contorno")]))
cat("direccion_eje_mayor: angulo del eje mayor medido desde el eje horizontal, en (-90, 90].\n")

cat("\nComprobacion independiente de cada Sigma_i y de la ecuacion del contorno\n")
S <- .krz_cov(cb, 1:2); q <- qchisq(0.95, 2)
for (i in 1:3) {
  e <- eigen(S[[i]], symmetric = TRUE)
  ce <- .contorno_elipse(cb$rows[i, 1:2], S[[i]], q)
  z <- rbind(ce$x - cb$rows[i, 1], ce$y - cb$rows[i, 2])
  cat(sprintf("  %-10s simetrica=%s  autovalores=(%.5f, %.5f)  max|z'S^-1 z - chi2|=%.1e\n", levels(g)[i],
              isSymmetric(S[[i]]), e$values[1], e$values[2], max(abs(colSums(z * solve(S[[i]], z)) - q))))
}

cat("\nComparacion con theta = lambda/(K-1) y con la version publicada\n")
rg_ms <- canonical_regions(bigarx_canonical(X, g, region = "ellipse", theta_est = "ms"))
rg_am <- canonical_regions(bigarx_canonical(X, g, region = "ellipse_amaro"))
print(r4(data.frame(grupo = rg$grupo, semieje_mayor_W = rg$semieje_mayor, semieje_mayor_W_ms = rg_ms$semieje_mayor,
                    semieje_mayor_publicada = rg_am$semieje_mayor, semieje_menor_W = rg$semieje_menor,
                    semieje_menor_publicada = rg_am$semieje_menor)))
cat("Discrepancias internas del articulo (no se ajusta el codigo a ellas): tabla 4, var. eje 2 de\n",
    "Versicolor = 0,2830 frente a 0,028 en su ecuacion de la elipse; tabla 7, semieje mayor de\n",
    "Virginica = 0,67 frente a 0,79 que se deduce de su propia varianza y su area; los centros\n",
    "de las tablas 6-7 no estan en la escala de W (en que se expresan las varianzas, con 1/n_i).\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  source("R/grafico.R")
  gg <- grafico_canonico(cb, opciones_grafico(cb, titulo = "Iris: elipses de Krzanowski en escala W",
                                              subtitulo = "Centros y contornos en la misma escala (coord_equal)"))
  ggplot2::ggsave("validacion/elipses_iris_escala_W.png", gg, width = 18, height = 13, units = "cm", dpi = 150)
  cat("\nGrafico guardado en validacion/elipses_iris_escala_W.png\n")
}
