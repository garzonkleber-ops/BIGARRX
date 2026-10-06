# =============================================================================
#  Comparación de BiGaRx Studio con los artículos de Amaro et al.
#  Ejecutar desde la carpeta BiGaRxApp:  source("comparar_articulos.R")
# =============================================================================
source("R/motor_canonico.R")
if (file.exists("R/boot_canonico.R")) source("R/boot_canonico.R")
r <- function(x, d = 4) round(x, d)
cat("\n===================== 1. IRIS: Amaro, Vicente-Villardón y Galindo (2008) =====================\n")
X <- iris[, 1:4]; g <- iris$Species; n <- nrow(X); K <- 3
cb <- bigarx_canonical(X, g, region = "ellipse", ndim = 10)

cat("\nTabla 1. 'Valores propios' del artículo = valores singulares en la escala de E = sqrt(theta)\n")
t1 <- data.frame(eje = 1:2, articulo = c(5.67, 0.53), app = sqrt(cb$canonical$valor_propio[1:2]),
                 inercia_art = c(99.12, 0.87), inercia_app = cb$inertia$inercia[1:2])
print(r(t1, 3))

cat("\nTabla 2. Calidad de las medias (%)\n")
t2 <- cbind(articulo_e1 = c(99.92, 86.28, 99.22), app_e1 = cb$quality_rows[, 1],
            articulo_e2 = c(0.08, 13.72, 0.78), app_e2 = cb$quality_rows[, 2])
print(r(t2, 2))

cat("\nTabla 3. Calidad de las variables (%)  [el artículo trae un 'eje 3' que no existe con K = 3]\n")
t3 <- cbind(art_e1 = c(4.96, 1.42, 49.85, 40.09), app_e1 = cb$quality_cols[, 1],
            art_e2 = c(9.66, 13.49, 40.89, 54.35), app_e2 = cb$quality_cols[, 2],
            art_e3 = c(83.29, 74.59, 2.81, 0.18),
            art_plano = c(14.62, 14.91, 90.74, 94.44), app_plano = rowSums(cb$quality_cols[, 1:2]))
print(r(t3, 2))

cat("\nComprobación directa: correlación al cuadrado (dentro de grupos) entre cada variable y cada eje\n")
Z <- scale(as.matrix(X), scale = FALSE)
S <- cb$individuals[, 1:2]                   # puntuaciones canónicas de los individuos
wc <- function(a) a - ave(a, g)               # centrar dentro de cada grupo
dir2 <- sapply(1:2, function(k) apply(Z, 2, function(v) cor(wc(v), wc(S[, k]))^2 * 100))
colnames(dir2) <- c("eje1", "eje2"); print(r(dir2, 2))

cat("\nTabla 9. Correlaciones variable-eje: el artículo usa la correlación TOTAL (sin quitar el grupo)\n")
tot <- sapply(1:2, function(k) apply(Z, 2, function(v) cor(v, S[, k])))
t9 <- cbind(art_e1 = c(0.7919, -0.5308, 0.9850, 0.9728), total_e1 = tot[, 1],
            art_e2 = c(0.2176, 0.7580, 0.0460, 0.2229), total_e2 = abs(tot[, 2]) * sign(c(0.2176, 0.7580, 0.0460, 0.2229)),
            app_dentro_e1 = cb$structure[, 1], app_dentro_e2 = cb$structure[, 2])
print(r(t9, 4))

cat("\nTablas 4, 5 y 7. Elipses de Krzanowski\n")
rg <- canonical_regions(cb)
t47 <- data.frame(grupo = c("G1", "G2", "G3"),
                  var1_art = c(0.0215, 0.0213, 0.0214), var1_app = rg$var1,
                  var2_art = c(0.1643, 0.2830, 0.1034), var2_app = rg$var2,
                  cov_art = c(-0.0042, -0.0034, 0.0075), cov_app = rg$cov,
                  mayor_art = c(0.99, 0.42, 0.67), mayor_app = rg$semieje_mayor,
                  menor_art = c(0.35, 0.34, 0.35), menor_app = rg$semieje_menor,
                  ang_art = c(1.7, 23.1, 5.21), ang_app = rg$angulo,
                  area_art = c(1.114, 0.418, 0.873), area_app = rg$area)
print(cbind(t47[1], r(t47[-1], 4)))
cat("Área con las propias cifras del artículo para G3: pi*0.67*0.35 =", r(pi * 0.67 * 0.35, 3),
    " ; pi*0.79*0.35 =", r(pi * 0.79 * 0.35, 3), "\n")

cat("\nCentros de las elipses (tabla 7) frente a las coordenadas de las medias\n")
cen <- data.frame(art_x = c(-2.3, 0.5, 1.7), art_y = c(-0.06, 0.22, -0.15),
                  app_W_x = cb$rows[, 1], app_W_y = cb$rows[, 2],
                  app_E_x = cb$rows[, 1] / sqrt(n - K), app_E_y = cb$rows[, 2] / sqrt(n - K))
print(r(cen, 3))
cat("Cociente artículo / escala E en el eje 1:", r(c(-2.3, 0.5, 1.7) / cen$app_E_x, 2),
    " ; (n-K)^(1/4) =", r((n - K)^0.25, 2), "\n")

cat("\nRadios de los círculos chi2 del artículo: sqrt(chi2_2 / n_i)\n")
cat("  escala W (app):", r(sqrt(qchisq(.95, 2) / 50), 4), "  escala E:", r(sqrt(qchisq(.95, 2) / 50) / sqrt(n - K), 4), "\n")

cat("\nMANOVA y comparación con R base\n")
mt <- manova_tests(cb); print(r(mt, 4))
sm <- summary(manova(as.matrix(X) ~ g), test = "Wilks")$stats
cat("Wilks de R base:", r(sm[1, 2], 6), " F:", r(sm[1, 3], 3), "\n")

cat("\n===================== 2. Validación empírica de las elipses (bootstrap) =====================\n")
set.seed(2026)
bt <- boot_canonical(X, g, B = 2000, axes = c(1, 2), seed = 2026)
vb <- sapply(levels(g), function(k) {
  P <- cbind(bt$replicates[, paste0("media|", k, "|Eje1")], bt$replicates[, paste0("media|", k, "|Eje2")])
  c(var1 = var(P[, 1]), var2 = var(P[, 2]), cov = cov(P)[1, 2])
})
cmpv <- data.frame(grupo = levels(g), var1_krz = rg$var1, var1_boot = vb["var1", ],
                   var2_krz = rg$var2, var2_boot = vb["var2", ], cov_krz = rg$cov, cov_boot = vb["cov", ])
print(cbind(cmpv[1], r(cmpv[-1], 4)))

cat("\n===================== 3. Arroz, Varela (1998), diapositivas 'BIPLOT CANÓNICO' =====================\n")
arroz <- data.frame(g = factor(rep(c("G1", "G2", "G3"), c(3, 3, 4))),
                    v1 = c(45, 50, 50, 60, 60, 65, 70, 65, 60, 65),
                    v2 = c(1.56, 1.60, 1.65, 1.75, 1.70, 1.70, 1.60, 1.60, 1.55, 1.70),
                    v3 = c(13, 12, 13, 15, 14, 14, 15, 13, 15, 14),
                    v4 = c(14, 16, 15, 9, 10, 7, 8, 13, 17, 11))
ca <- bigarx_canonical(arroz[, -1], arroz$g, ndim = 10)
print(r(data.frame(eje = 1:2, articulo = c(3.8, 1.4), app_sqrt_theta = sqrt(ca$canonical$valor_propio[1:2]),
                   inercia_art = c(87.1, 12.9), inercia_app = ca$inertia$inercia[1:2]), 3))

cat("\n===================== 4. Césped, Seber (1984) - Amaro et al. (2004), Biplot TOTAL =====================\n")
temp <- rep(c("T14", "T18", "T22", "T26"), each = 6)
ctrl <- matrix(c(2.2,23.5,1.7, 3.0,27.0,2.3, 3.3,24.5,3.2, 2.2,20.5,1.5, 2.0,19.0,2.0, 3.5,23.5,2.9,
                 21.8,41.5,23.0, 11.0,32.5,15.4, 16.4,46.5,22.8, 13.1,31.0,21.5, 15.4,41.5,20.8, 14.5,46.0,20.3,
                 13.6,29.5,30.8, 6.2,23.5,14.6, 16.7,58.5,36.0, 12.2,40.5,23.9, 8.7,37.0,20.3, 12.3,41.5,27.7,
                 3.0,24.0,10.2, 5.3,26.5,15.6, 3.1,24.5,14.7, 4.8,34.0,20.5, 3.4,22.5,14.3, 7.4,32.0,23.2),
               ncol = 3, byrow = TRUE)
inf <- matrix(c(2.3,23.5,2.0, 3.0,21.0,2.7, 2.3,22.0,1.8, 2.5,22.5,2.4, 2.4,21.5,1.1, 2.7,25.0,2.6,
                10.1,43.5,14.2, 7.6,27.0,14.7, 19.7,32.5,21.4, 4.3,28.5,9.7, 5.2,33.5,12.2, 3.9,24.5,8.2,
                10.0,21.0,23.6, 12.3,49.0,28.1, 4.9,28.5,13.3, 9.6,27.0,24.6, 6.5,29.0,19.3, 13.6,30.5,31.5,
                4.2,25.5,13.3, 2.2,23.5,8.5, 2.8,19.5,11.8, 1.3,21.5,7.8, 4.2,28.5,15.1, 3.0,25.0,11.8),
              ncol = 3, byrow = TRUE)
cesped <- data.frame(trat = rep(c("GC", "GI"), each = 24), temp = c(temp, temp), rbind(ctrl, inf))
names(cesped)[3:5] <- c("X1", "X2", "X3")
celda <- interaction(cesped$temp, cesped$trat, sep = "-")
for (tr in c("sin log", "log")) {
  Xc <- cesped[, 3:5]; if (tr == "log") Xc <- log(Xc)
  ct <- bigarx_canonical(Xc, celda, ndim = 10)
  m <- manova_tests(ct)
  cat(sprintf("\nDatos %s: sqrt(theta) = %s ; inercia = %s ; Wilks = %.4f\n", tr,
              paste(r(sqrt(ct$canonical$valor_propio), 4), collapse = ", "),
              paste(r(ct$inertia$inercia, 4), collapse = ", "), m["Wilks", "estadistico"]))
  cat("   escala W, d = sqrt((n-K)*theta):", r(sqrt(ct$canonical$valor_propio * (48 - 8)), 4), "\n")
}
cat("Artículo (tabla IV): 120.1159, 15.3571, 0.4073 ; inercia 98.3906, 1.6083, 0.0011 ; Wilks global 0.0105\n")
fit2 <- manova(as.matrix(log(cesped[, 3:5])) ~ trat * temp, data = cesped)
cat("\nMANOVA de dos factores (log), Wilks de R base [tabla I del artículo: 0.765, 0.0126, 0.8017]:\n")
print(r(summary(fit2, test = "Wilks")$stats[1:3, 1:2], 4))
fit2b <- manova(as.matrix(cesped[, 3:5]) ~ trat * temp, data = cesped)
cat("Mismo MANOVA sin log:\n"); print(r(summary(fit2b, test = "Wilks")$stats[1:3, 1:2], 4))
