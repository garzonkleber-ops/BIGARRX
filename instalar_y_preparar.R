# =============================================================================
#  BiGaRx Studio - preparar la app en tu Mac
#  1) Abre este archivo en RStudio con la carpeta BiGaRxApp como directorio de
#     trabajo (Session > Set Working Directory > To Source File Location).
#  2) Ejecuta todo (Cmd + Shift + Enter).
# =============================================================================
paquetes <- c("shiny", "ggplot2", "ggrepel", "rhandsontable", "colourpicker", "shinyjs",
              "plotly", "readxl", "writexl", "svglite", "rsconnect")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) install.packages(faltan, repos = "https://cloud.r-project.org")

# Probar la app en tu computador (se abre en el navegador; cierra con Esc en la consola)
# shiny::runApp(".")

# Generar el manifiesto que necesita Posit Connect Cloud (vuelve a ejecutarlo
# cada vez que agregues paquetes nuevos a la app)
rsconnect::writeManifest()
cat("Listo: manifest.json creado en", getwd(), "\n")
