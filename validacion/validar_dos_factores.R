# Validacion del MANOVA biplot de dos factores con Seber (1984) / Amaro et al. (2004)
# Ejecutar con el directorio de trabajo en la carpeta BiGaRxApp.
for (f in list.files("R", full.names = TRUE)) source(f)
ces <- read.csv("datos/cesped_seber.csv")
X <- log(ces[, c("X1", "X2", "X3")])
r <- function(x, d = 4) { if (is.data.frame(x)) { x[] <- lapply(x, function(v) if (is.numeric(v)) round(v, d) else v); x } else round(x, d) }
ajustes <- lapply(c(total = "total", filas = "filas", columnas = "columnas", interaccion = "interaccion"),
                  function(e) bigarx_manova2(X, ces$tratamiento, ces$temperatura, efecto = e,
                                             nombre_A = "Tratamiento", nombre_B = "Temperatura"))
cat("Tabla I (Wilks, chi2):  Filas 0,765 9,7781 | Columnas 0,0126 164,1014 | Interaccion 0,8017 8,2883\n")
print(r(ajustes$total$diseno$pruebas[ajustes$total$diseno$pruebas$contraste == "Wilks",
                                      c("efecto", "estadistico", "F_aprox", "p_valor", "chi2_bartlett")]))
prop <- function(cb) { th <- cb$canonical$valor_propio; c(theta_rel = th / th[1], inercia_theta2 = 100 * th^2 / sum(th^2), inercia = cb$inertia$inercia) }
cat("\nTabla II (interaccion): 5,8399 2,9939 0,0847 -> relativos 1 0,5127 0,0145 ; inercia 79,175 20,808 0,017\n")
print(r(prop(ajustes$interaccion)))
cat("\nTabla III (calidades de las medias en el biplot de interaccion), articulo T1-GC 3,1963 44,2412 ; T2-GC 94,3761 0,0173\n")
print(r(ajustes$interaccion$quality_rows[, 1:2], 4))
cat("\nTabla IV (total): 120,1159 15,3571 0,4073 -> relativos 1 0,1279 0,0034 ; inercia 98,3906 1,6083 0,0011\n")
print(r(prop(ajustes$total)))
cat("\nTabla V (total), calidades de las medias y variables (coordenadas principales)\n")
print(r(ajustes$total$quality_rows[, 1:2], 4)); print(r(ajustes$total$quality_cols_principal[, 1:2], 4))
cat("\nTabla VI (columnas): 279,5514 30,4903 0,1483 -> relativos 1 0,1091 0,00053 ; inercia 98,8244 1,1756 0\n")
print(r(prop(ajustes$columnas), 5))
cat("\nTabla VII (columnas): T1 99,6234 0,3748 | T2 13,2106 86,2463 | X1 81,6214 18,3786\n")
print(r(ajustes$columnas$quality_rows[, 1:2], 4)); print(r(ajustes$columnas$quality_cols_principal[, 1:2], 4))
cat("\nFilas (1 g.l.): medias proyectadas sobre el biplot total (figura 5)\n")
print(r(ajustes$filas$rows[, 1:2], 4)); cat("proyectado:", ajustes$filas$diseno$proyectado, "\n")
cat("\nT2 de Hotelling por pares (columnas)\n"); print(r(hotelling_pares(ajustes$columnas)[, -c(1, 2)], 4))
