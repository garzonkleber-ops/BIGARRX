# Regiones elípticas del biplot canónico: corrección de escala

Respuesta a la revisión del 6 de octubre de 2026 sobre las elipses de Amaro, Vicente-Villardón y Galindo-Villardón (2008).

## 1. Problema

La versión anterior colocaba cada elipse alrededor del centroide `P_i` en la escala de W (W = E / (n − K)). Pero calculaba sus varianzas con dos cantidades de otra escala:

- **Coordenadas:** las ecuaciones 3 y 4 usaban `xi = P / sqrt(n − K)`, en la escala de E. El término `1/n_i` sí estaba en la escala de W, así que los términos que dependen de ξ quedaban divididos entre 147 y la elipse se reducía casi a un círculo de varianza `1/n_i`.
- **θ:** se calculaba como `1 + (d / sqrt(n − K)) / (K − 1)`, es decir, con un valor singular en la escala de E (el 5,67 que el artículo llama "valor propio"). θ es un **valor propio**.

Esa versión reproducía las tablas 4 a 7 del artículo, que también mezclan escalas: sus centros (−2,3; 0,5; 1,7) no están en la escala de sus varianzas. Se conserva como opción `region = "ellipse_amaro"`, marcada como no coherente.

## 2. Corrección (`region = "ellipse"`, por defecto)

Se usa una sola escala, la de W, con varianza dentro de los grupos igual a 1. Es la escala de `x$rows`, de los individuos y del término `1/n_i`, y se conserva `coord_equal()` en el gráfico.

| Símbolo (artículo) | Cantidad en el código | Iris, eje 1 | Iris, eje 2 |
|---|---|---|---|
| ξ_ij | `P_ij` = `x$rows`, escala de W | −7,608 (setosa) | 0,215 |
| λ_j | `d_j²` = autovalor de SCH·W⁻¹ (= W^−½ SCH W^−½) | 4732,2 | 41,95 |
| θ_j | `1 + λ_j/(K−1)` (opción `uno_mas`, texto del artículo, p. 234) | 2367,1 | 21,98 |
| θ_j (alternativa) | `λ_j/(K−1)` = autovalor de (SCH/gl_H)·W⁻¹ (opción `ms`) | 2366,1 | 20,98 |
| γ | `(K−1)/(n−K)` = gl_H / gl_E | 0,01361 | |
| Valor propio del MANOVA | `θ_HE = d²/(n−K)` = autovalor de H·E⁻¹ | 32,19 | 0,285 |
| "Valor propio" de la tabla 1 del artículo | `sqrt(θ_HE)`: valor singular en escala E | 5,674 | 0,534 |
| θ de la versión publicada | `1 + sqrt(θ_HE)/(K−1)` | 3,837 | 1,267 |

Ecuaciones, con l el eje de interés, j ≠ l y gl_E = n − K:

```
Var(P_il)       = 1/n_i + [ ξ_il²/2 + Σ_j θ_l(θ_j + γθ_l) / (γ(θ_j − θ_l)²) · ξ_ij² ] / gl_E
Cov(P_ij, P_il) = − θ_jθ_l(1 + γ) / (gl_E · γ(θ_j − θ_l)²) · ξ_ij ξ_il
```

El contorno es `(z − P_i)' Σ_i⁻¹ (z − P_i) = χ²(2; 0,95)` y se dibuja con `Σ_i = U D Uᵀ`.

El código comprueba que cada Σ_i sea simétrica y definida positiva. Si no lo es, la región no se dibuja y se avisa; no se fuerzan valores propios a cero. También comprueba que cada punto del contorno cumpla la ecuación (error < 10⁻⁸·χ²).

Para diseños de dos factores se sustituyen K − 1 por gl_H del efecto y n − K por gl_E del modelo completo.

## 3. Verificación con Iris

Iris tiene 150 individuos, 3 grupos y 4 variables, centrados. Todo está en la escala de W, con nivel del 95 % y θ = 1 + λ/(K−1).

| Grupo | Centroide | Var. eje 1 | Var. eje 2 | Cov. | Semieje mayor | Semieje menor | Orientación del eje mayor |
|---|---|---|---|---|---|---|---|
| Setosa | (−7,6076; 0,2151) | 0,2174 | 0,2939 | 0,0078 | 1,3288 | 1,1392 | 84,2° |
| Versicolor | (1,8250; −0,7279) | 0,0375 | 0,0376 | 0,0064 | 0,5129 | 0,4321 | 45,1° |
| Virginica | (5,7826; 0,5128) | 0,1368 | 0,1791 | −0,0142 | 1,0482 | 0,8909 | −73,0° |

Las tres matrices son definidas positivas, y el error máximo de la ecuación del contorno es de 10⁻¹⁴.

La versión anterior daba semiejes de 0,996 / 0,358, 0,422 / 0,347 y 0,792 / 0,353, que coinciden con la tabla de la revisión.

**Discrepancias internas del artículo.** No se ajustó el código para reproducirlas:

- Tabla 4: la varianza del eje 2 de Versicolor es 0,2830, pero su ecuación de la elipse usa 0,028.
- Tabla 7: el semieje mayor de Virginica es 0,67, pero de su propia varianza y de su área se deduce 0,79.
- Tablas 6 y 7: los centros no están en la escala de W.

Script: `validacion/verificar_elipses.R`. Gráfico: `validacion/elipses_iris_escala_W.png`.

## 4. Simulación de cobertura revisada

Script: `validacion/simulacion_cobertura.R`.

- **Población:** normal multivariante con las medias y la W combinada de los datos, y los mismos tamaños de grupo.
- **Parámetro evaluado:** ξ_i, la media poblacional del grupo i proyectada sobre los ejes canónicos **poblacionales**, en la escala de W.
- **Alineación:** cada eje muestral cambia de signo para coincidir con el poblacional, `sign(Σ_i P_ik ξ_ik)`. No se rota.

Cobertura al 95 % nominal: 2000 muestras, y 200 en el caso del bootstrap (B = 300).

| | Krzanowski W (1+λ/(K−1)) | Krzanowski W (λ/(K−1)) | Amaro publicada | Círculo χ² | Círculo t | Bootstrap |
|---|---|---|---|---|---|---|
| Iris setosa | 98,5 | 98,5 | 57,0 | 55,2 | 46,3 | 90,0 |
| Iris versicolor | 99,0 | 99,0 | 94,3 | 92,6 | 83,0 | 92,5 |
| Iris virginica | 98,6 | 98,6 | 67,3 | 65,0 | 54,3 | 88,0 |
| Cerámica A | 96,4 | 96,4 | 59,1 | 48,4 | 39,8 | 94,5 |
| Cerámica C | 92,8 | 92,7 | 71,7 | 60,0 | 49,4 | 86,5 |
| Cerámica I | 96,0 | 96,0 | 56,1 | 44,5 | 37,0 | 93,0 |
| Cerámica L | 96,6 | 96,6 | 56,3 | 47,9 | 40,1 | 90,5 |

**Conclusión.** La baja cobertura que se había atribuido al método de Amaro y Krzanowski se debía a la implementación con escalas mezcladas. Con una sola escala, las elipses quedan cerca del 95 % nominal o por encima. Por eso son la opción recomendada en la app.

**Limitación observada.** Comparando la varianza de la fórmula con la varianza real de P_i:

- En Iris, el eje 1 coincide (0,217 frente a 0,213 en setosa), pero el eje 2 se sobrestima (0,294 frente a 0,0015).
- En la cerámica (n = 26) ocurre al revés: el eje 1 se subestima y el eje 2 se sobrestima.

La cobertura conjunta sigue siendo aceptable. Una hipótesis, no verificada, es que los ejes se estiman con las mismas medias que se proyectan, y la aproximación asintótica no recoge esa dependencia cuando hay pocos grupos o pocos individuos.

## Referencias

- Amaro, I. R., Vicente-Villardón, J. L. y Galindo-Villardón, M. P. (2008). Contribuciones al MANOVA-Biplot: regiones de confianza alternativas. *Investigación Operacional*, 29(3), 231–241.
- Anderson, T. W. (1984). *An Introduction to Multivariate Statistical Analysis* (2.ª ed.). Wiley.
- Krzanowski, W. J. (1989). On confidence regions in canonical variate analysis. *Biometrika*, 76(1), 107–116.
