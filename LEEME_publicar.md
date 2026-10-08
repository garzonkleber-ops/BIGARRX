# BiGaRx Studio — app Shiny estilo jamovi

Biplot canónico / MANOVA biplot con imputación de faltantes, bootstrap y un gráfico 100 % personalizable.

## Estructura

| Archivo | Qué hace |
|---|---|
| `app.R` | Interfaz (pestañas Datos / Análisis, cinta, panel de opciones, resultados) y servidor |
| `R/motor_canonico.R` | Cálculo del biplot canónico, MANOVA, regiones de confianza, imputación (de BiGaRx) |
| `R/boot_canonico.R` | Bootstrap estratificado por grupo |
| `R/grafico.R` | Gráfico ggplot2 con todas las opciones de personalización |
| `R/motor_dos_factores.R` | MANOVA biplot de dos factores, T² de Hotelling por pares, datos composicionales |
| `R/tablas.R` | Tablas con estilo APA |
| `validacion/` | Scripts que reproducen las tablas de Amaro et al. (2004, 2008) y la simulación de cobertura |
| `www/estilo.css` | Diseño visual |
| `datos/` | Ejemplos: Iris, cerámica romana, PISA 2025, césped (Seber, dos factores), arroz (Varela) y uso del suelo simulado (composicional) |
| `instalar_y_preparar.R` | Instala los paquetes y crea `manifest.json` |

## 1. Probar en tu Mac

1. Abre `instalar_y_preparar.R` en RStudio y ejecútalo (instala paquetes y crea `manifest.json`).
2. En la consola: `shiny::runApp("ruta/a/BiGaRxApp")`.

## 2. Publicar en Posit Connect Cloud

1. Crea un repositorio **público** en GitHub (por ejemplo `bigarx-studio`) y sube el
   contenido de la carpeta `BiGaRxApp` (incluido `manifest.json`). Puedes hacerlo
   desde github.com → *Add file → Upload files* arrastrando los archivos y carpetas.
2. Entra en <https://connect.posit.cloud> con tu cuenta.
3. *Publish* → *Shiny* → elige el repositorio, la rama `main` y el archivo principal `app.R`.
4. Pulsa *Publish*. La primera vez tarda unos minutos porque instala los paquetes.
5. Para actualizar: sube los cambios a GitHub y pulsa *Republish* (o activa la
   republicación automática).

Si cambias de versión de R o añades paquetes, vuelve a ejecutar `rsconnect::writeManifest()`
y sube el nuevo `manifest.json`.

## Qué se puede personalizar en el gráfico

- Plano (cualquier par de ejes), qué grupos y **qué variables (flechas)** mostrar, filtro por calidad.
- Flechas: color único o uno por variable, grosor, punta, escala y multiplicador de longitud.
- Grupos: paleta o color manual por grupo, forma/tamaño/opacidad de individuos y medias.
- Regiones: elipses de Krzanowski o círculos chi², t, Bonferroni; nivel, grosor, tipo de línea, relleno.
- Etiquetas de grupos, variables e individuos con tamaños propios y anti‑solapamiento.
- Título, subtítulo, ejes, nota al pie, leyenda, tema, fuente, tamaño de letra, fondo, cuadrícula, límites.
- Bootstrap: elipses de confianza de medias y flechas, nubes de réplicas.
- Modo interactivo (plotly) y exportación PNG/PDF/SVG en cm y DPI elegidos; tablas a Excel y coordenadas a CSV.

## Novedades de la versión 2

- **Diseño de dos factores** (Amaro et al., 2004): biplot total, de filas, de columnas y de interacción, con MANOVA y ANOVA por efecto.
- **Elipses bootstrap** como tipo de región (recomendado para ubicar las medias).
- **Calidad de las variables** con la definición de 2008 (correlación²) o la de 2004 (coordenadas principales), y tabla de correlaciones totales.
- **Datos composicionales**: transformaciones clr (con inversa generalizada) y alr; aviso automático.
- **T² de Hotelling por pares** con ajuste de Bonferroni u Holm.
- **χ² de Bartlett**, valores singulares (escala E), inercia θ² opcional, invertir ejes, logaritmo de las variables.

## Versión 2.1 (7 de octubre de 2026): corrección de las elipses

- Las elipses de Krzanowski se calculan en una sola escala, la de W: centros, varianzas y covarianzas. Es la opción recomendada.
- La versión publicada por Amaro et al. (2008) queda como opción aparte, marcada como no coherente, solo para reproducir sus tablas.
- Detalles, tabla de verificación con Iris y simulación revisada: `validacion/ELIPSES_escala_W.md`.
