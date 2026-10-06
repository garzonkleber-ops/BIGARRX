# =============================================================================
#  BiGaRx Studio - Biplot canonico / MANOVA biplot
#  Aplicacion Shiny con diseno tipo jamovi. Los calculos estan en la carpeta R/
#  (Shiny la carga automaticamente).
# =============================================================================
library(shiny)
library(ggplot2)
library(rhandsontable)
library(shinyjs)
library(colourpicker)

options(shiny.maxRequestSize = 30 * 1024^2)

ejemplos <- c("Iris (especies)" = "iris.csv",
              "Cerámica romana (POTTERY)" = "ceramica_romana.csv",
              "PISA 2025 (regiones)" = "pisa2025.csv",
              "Césped, Seber 1984 (dos factores)" = "cesped_seber.csv",
              "Arroz, Varela 1998 (un factor)" = "arroz_varela.csv",
              "Uso del suelo SIMULADO (composicional)" = "uso_suelo_simulado.csv")

seccion <- function(titulo, ..., abierta = FALSE, id = NULL) {
  tags$details(class = "op-seccion", open = if (abierta) NA else NULL, id = id,
               tags$summary(titulo), div(class = "op-cuerpo", ...))
}

boton_cinta <- function(id, icono, texto) {
  actionButton(id, label = tagList(icon(icono), tags$span(texto)), class = "cinta-btn")
}

# -----------------------------------------------------------------------------
# Interfaz
# -----------------------------------------------------------------------------
ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$link(rel = "stylesheet", href = "estilo.css"),
            tags$title("BiGaRx Studio")),

  # Barra superior con pestanas
  div(class = "barra",
      div(class = "marca", tags$b("BiGaRx"), tags$span("Studio")),
      div(class = "pestanas",
          actionLink("tab_datos", "Datos", class = "pestana activa"),
          actionLink("tab_analisis", "Análisis", class = "pestana")),
      div(class = "barra-der", "Biplot canónico · MANOVA")),

  # Cinta de herramientas
  div(class = "cinta", id = "cinta_datos",
      div(class = "cinta-grupo",
          fileInput("archivo", NULL, accept = c(".csv", ".txt", ".xlsx", ".xls"),
                    buttonLabel = tagList(icon("folder-open"), "Abrir"), placeholder = "CSV o Excel"),
          uiOutput("hoja_excel_ui")),
      div(class = "cinta-grupo",
          selectInput("ejemplo", NULL, choices = c("Ejemplos..." = "", ejemplos), width = "210px")),
      div(class = "cinta-grupo",
          boton_cinta("btn_nueva_col", "table-columns", "Nueva columna"),
          boton_cinta("btn_nueva_fila", "plus", "Nueva fila"),
          downloadButton("bajar_datos", label = "Guardar datos", class = "cinta-btn"))),
  hidden(div(class = "cinta", id = "cinta_analisis",
      div(class = "cinta-grupo",
          boton_cinta("btn_canonico", "chart-line", "Biplot canónico"),
          boton_cinta("btn_manova", "table", "MANOVA"),
          boton_cinta("btn_boot", "shuffle", "Bootstrap")))),

  div(class = "principal",
      # ---- Panel izquierdo: hoja de datos u opciones del analisis ----
      div(class = "panel-izq", id = "panel_izq",
          div(id = "vista_datos",
              div(class = "hoja-info", textOutput("info_datos")),
              rHandsontableOutput("hoja", height = "calc(100vh - 190px)")),
          hidden(div(id = "vista_opciones", class = "opciones",
              div(class = "op-cabecera",
                  tags$span("Biplot canónico / MANOVA"),
                  actionLink("cerrar_opciones", icon("arrow-right"), title = "Ver la hoja de datos")),

              seccion("Variables y diseño", abierta = TRUE,
                  selectizeInput("vars", "Variables dependientes", choices = NULL, multiple = TRUE,
                                 options = list(plugins = list("remove_button"))),
                  radioButtons("diseno", "Diseño", c("Un factor" = "uno", "Dos factores" = "dos"), inline = TRUE),
                  selectInput("grupo", "Factor de agrupación (filas)", choices = NULL),
                  conditionalPanel("input.diseno == 'dos'",
                      selectInput("grupo2", "Segundo factor (columnas)", choices = NULL),
                      radioButtons("efecto", "Biplot a representar",
                                   c("Total (todas las celdas)" = "total", "Filas (1.er factor)" = "filas",
                                     "Columnas (2.º factor)" = "columnas", "Interacción" = "interaccion")),
                      tags$p(class = "op-ayuda", "Amaro, Vicente-Villardón y Galindo (2004). Las sumas de cuadrados",
                             "son secuenciales (tipo I); en diseños balanceados no depende del orden."))),

              seccion("Modelo", abierta = TRUE,
                  radioButtons("transform", "Transformación",
                               c("Centrar columnas" = "center", "Estandarizar columnas" = "standardize"),
                               inline = TRUE),
                  radioButtons("na_action", "Datos faltantes",
                               c("Eliminar filas con NA" = "omit", "Imputar (SVD iterativa)" = "impute"),
                               inline = TRUE),
                  conditionalPanel("input.na_action == 'impute'",
                      numericInput("impute_rank", "Rango de la imputación", 2, min = 1, max = 10)),
                  checkboxInput("log_vars", "Aplicar logaritmo natural a las variables", FALSE),
                  selectInput("composicional", "Datos composicionales (proporciones)",
                              c("Sin transformar" = "ninguna",
                                "clr: log-cociente centrado (+ inversa generalizada)" = "clr",
                                "alr: log-cociente respecto a una variable" = "alr")),
                  conditionalPanel("input.composicional == 'alr'",
                      selectInput("alr_ref", "Variable de referencia (denominador)", choices = NULL)),
                  checkboxInput("ginv", "Permitir W singular (inversa generalizada)", FALSE)),

              seccion("Regiones de confianza", abierta = TRUE,
                  selectInput("region", "Tipo",
                              c("Elipses bootstrap (recomendado)" = "boot",
                                "Elipses (Krzanowski / Amaro 2008)" = "ellipse",
                                "Círculos chi²" = "chi2", "Círculos t" = "t",
                                "Círculos Bonferroni" = "bonferroni"), selected = "ellipse"),
                  sliderInput("conf", "Nivel de confianza", 0.80, 0.99, 0.95, 0.01),
                  tags$p(class = "op-ayuda", "En una simulación con Iris, las elipses de Krzanowski y los",
                         "círculos cubrieron solo el 57-67 % de las medias alejadas del origen; las elipses",
                         "bootstrap reproducen la variabilidad real. Los círculos siguen siendo útiles para",
                         "comparar pares de grupos."),
                  conditionalPanel("input.region == 'boot'",
                      tags$p(class = "op-ayuda", "Usa las réplicas, la semilla y el tipo de intervalo de la sección",
                             "Bootstrap; se calcula automáticamente."))),

              seccion("Resultados a mostrar",
                  checkboxGroupInput("tablas", NULL,
                      c("Descriptivos por grupo" = "desc", "Contrastes MANOVA" = "manova",
                        "ANOVA por variable" = "anova", "Valores propios e inercia" = "eigen",
                        "Medias de grupo (coordenadas y calidad)" = "medias",
                        "Variables (calidad y correlaciones)" = "variables",
                        "Correlaciones totales variable-eje" = "corr_total",
                        "Regiones de confianza" = "regiones",
                        "Comparaciones por pares (T² de Hotelling)" = "pares", "Imputación" = "imput"),
                      selected = c("manova", "anova", "eigen", "medias", "variables", "regiones", "pares", "imput")),
                  selectInput("ajuste_p", "Ajuste de p en las comparaciones por pares",
                              c("Bonferroni" = "bonferroni", "Holm" = "holm", "Ninguno" = "none")),
                  checkboxInput("inercia_amaro", "Añadir la inercia con θ² (convención de Amaro et al., 2004)", FALSE)),

              seccion("Bootstrap", id = "sec_boot",
                  checkboxInput("boot_on", "Activar bootstrap", FALSE),
                  conditionalPanel("input.boot_on",
                      numericInput("boot_B", "Réplicas", 1000, min = 100, max = 10000, step = 100),
                      numericInput("boot_seed", "Semilla", 2026),
                      radioButtons("boot_ci", "Intervalos", c("Percentil" = "percentile", "Normal" = "normal"),
                                   inline = TRUE),
                      sliderInput("boot_conf", "Nivel", 0.80, 0.99, 0.95, 0.01),
                      actionButton("boot_run", "Ejecutar bootstrap", icon = icon("play"), class = "btn-accion"),
                      tags$hr(),
                      checkboxInput("boot_medias", "Elipses bootstrap de las medias", TRUE),
                      checkboxInput("boot_variables", "Elipses bootstrap de las flechas", TRUE),
                      checkboxInput("boot_nubes", "Mostrar las réplicas (nubes)", FALSE),
                      colourpicker::colourInput("boot_color", "Color de las elipses de las flechas", "#C2410C"))),

              seccion("Gráfico: plano y elementos",
                  fluidRow(column(6, numericInput("eje_x", "Eje X", 1, min = 1, max = 10)),
                           column(6, numericInput("eje_y", "Eje Y", 2, min = 1, max = 10))),
                  fluidRow(column(6, checkboxInput("invertir_x", "Invertir eje X", FALSE)),
                           column(6, checkboxInput("invertir_y", "Invertir eje Y", FALSE))),
                  checkboxInput("ver_individuos", "Individuos", TRUE),
                  checkboxInput("ver_medias", "Medias de grupo", TRUE),
                  checkboxInput("ver_regiones", "Regiones de confianza", TRUE),
                  checkboxInput("ver_flechas", "Flechas de las variables", TRUE),
                  radioButtons("modo", "Modo", c("Estático" = "estatico", "Interactivo" = "interactivo"),
                               inline = TRUE)),

              seccion("Gráfico: variables (flechas)",
                  checkboxGroupInput("vars_mostrar", "Variables a mostrar", choices = NULL),
                  radioButtons("calidad_def", "Calidad de las variables",
                               c("Correlación² dentro de los grupos (Amaro et al., 2008)" = "corr",
                                 "Coordenadas principales (Amaro et al., 2004)" = "principal")),
                  sliderInput("min_calidad", "Calidad mínima en el plano (%)", 0, 100, 0, 5),
                  radioButtons("flecha_escala", "Escala de las flechas",
                               c("Automática" = "auto", "Original" = "orig"), inline = TRUE),
                  sliderInput("flecha_mult", "Multiplicador de longitud", 0.1, 3, 1, 0.05),
                  colourpicker::colourInput("flecha_color", "Color de las flechas", "#333333"),
                  checkboxInput("colores_por_var", "Un color por variable", FALSE),
                  uiOutput("colores_var_ui"),
                  sliderInput("flecha_grosor", "Grosor", 0.2, 3, 0.7, 0.1),
                  sliderInput("flecha_punta", "Tamaño de la punta (cm)", 0, 0.6, 0.2, 0.02)),

              seccion("Gráfico: grupos",
                  checkboxGroupInput("grupos_mostrar", "Grupos a mostrar", choices = NULL),
                  selectInput("paleta", "Paleta", paletas_disponibles),
                  checkboxInput("colores_manual", "Elegir el color de cada grupo", FALSE),
                  uiOutput("colores_grupo_ui"),
                  sliderInput("ind_tam", "Tamaño de los individuos", 0.2, 6, 1.6, 0.1),
                  sliderInput("ind_alfa", "Opacidad de los individuos", 0.05, 1, 0.45, 0.05),
                  selectInput("ind_forma", "Forma de los individuos", formas, selected = 16),
                  checkboxInput("ind_forma_por_grupo", "Una forma por grupo", FALSE),
                  sliderInput("media_tam", "Tamaño de las medias", 1, 10, 4, 0.5),
                  selectInput("media_forma", "Forma de las medias", formas, selected = 15),
                  sliderInput("reg_grosor", "Grosor de las regiones", 0.2, 3, 0.9, 0.1),
                  selectInput("reg_linea", "Tipo de línea",
                              c("Continua" = "solid", "Discontinua" = "dashed", "Punteada" = "dotted",
                                "Raya-punto" = "dotdash")),
                  checkboxInput("reg_relleno", "Rellenar las regiones", FALSE),
                  sliderInput("reg_alfa", "Opacidad del relleno", 0.02, 0.6, 0.12, 0.02)),

              seccion("Gráfico: etiquetas",
                  checkboxInput("etiq_medias", "Nombres de los grupos", TRUE),
                  checkboxInput("etiq_variables", "Nombres de las variables", TRUE),
                  checkboxInput("etiq_individuos", "Nombres de los individuos", FALSE),
                  checkboxInput("evitar_solape", "Evitar que se superpongan", TRUE),
                  sliderInput("etiq_media_tam", "Tamaño (grupos)", 1, 10, 4, 0.2),
                  sliderInput("etiq_var_tam", "Tamaño (variables)", 1, 10, 3.4, 0.2),
                  sliderInput("etiq_ind_tam", "Tamaño (individuos)", 1, 8, 2.4, 0.2),
                  checkboxInput("etiq_var_mismo", "Etiquetas de variables del color de su flecha", TRUE),
                  conditionalPanel("!input.etiq_var_mismo",
                      colourpicker::colourInput("etiq_var_color", "Color de las etiquetas de variables", "#000000"))),

              seccion("Gráfico: textos y estilo",
                  textInput("titulo", "Título", "Biplot canónico"),
                  textInput("subtitulo", "Subtítulo", ""),
                  textInput("xlab", "Eje X (vacío = automático)", ""),
                  textInput("ylab", "Eje Y (vacío = automático)", ""),
                  textInput("pie", "Nota al pie", ""),
                  textInput("titulo_leyenda", "Título de la leyenda", "Grupo"),
                  selectInput("tema", "Tema", temas),
                  sliderInput("tam_base", "Tamaño de letra", 8, 24, 13, 1),
                  selectInput("fuente", "Fuente", c("Sans" = "sans", "Serif" = "serif", "Mono" = "mono")),
                  selectInput("leyenda", "Leyenda", c("Derecha" = "right", "Abajo" = "bottom",
                                                      "Arriba" = "top", "Izquierda" = "left", "Sin leyenda" = "none")),
                  checkboxInput("cuadricula", "Cuadrícula", TRUE),
                  checkboxInput("lineas_cero", "Líneas en el origen", TRUE),
                  colourpicker::colourInput("fondo", "Color de fondo", "#FFFFFF"),
                  tags$label("Límites (vacío = automático)"),
                  fluidRow(column(3, numericInput("xmin", "x mín", NA)), column(3, numericInput("xmax", "x máx", NA)),
                           column(3, numericInput("ymin", "y mín", NA)), column(3, numericInput("ymax", "y máx", NA)))),

              seccion("Exportar",
                  fluidRow(column(6, numericInput("ancho", "Ancho (px en pantalla)", 760, 300, 3000, 10)),
                           column(6, numericInput("alto", "Alto (px en pantalla)", 620, 300, 3000, 10))),
                  fluidRow(column(4, numericInput("exp_ancho", "Ancho (cm)", 18, 5, 60)),
                           column(4, numericInput("exp_alto", "Alto (cm)", 15, 5, 60)),
                           column(4, numericInput("exp_dpi", "DPI", 300, 72, 1200))),
                  selectInput("exp_formato", "Formato", c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg")),
                  downloadButton("bajar_grafico", "Descargar gráfico"),
                  downloadButton("bajar_tablas", "Descargar tablas (Excel)"),
                  downloadButton("bajar_coords", "Descargar coordenadas (CSV)"))
          ))),

      # ---- Panel derecho: resultados ----
      div(class = "panel-der",
          div(class = "res-cabecera", "Resultados"),
          uiOutput("resultados"))
  )
)

# -----------------------------------------------------------------------------
# Servidor
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  datos <- reactiveVal(read.csv("datos/iris.csv", check.names = FALSE))
  nombre_datos <- reactiveVal("iris.csv")

  # ---- Pestanas y paneles ----
  mostrar_datos <- function() {
    show("cinta_datos"); hide("cinta_analisis"); show("vista_datos"); hide("vista_opciones")
    addClass("tab_datos", "activa"); removeClass("tab_analisis", "activa")
  }
  mostrar_analisis <- function() {
    hide("cinta_datos"); show("cinta_analisis"); hide("vista_datos"); show("vista_opciones")
    addClass("tab_analisis", "activa"); removeClass("tab_datos", "activa")
  }
  observeEvent(input$tab_datos, mostrar_datos())
  observeEvent(input$cerrar_opciones, mostrar_datos())
  observeEvent(input$tab_analisis, mostrar_analisis())
  observeEvent(input$btn_canonico, mostrar_analisis())
  observeEvent(input$btn_manova, {
    mostrar_analisis()
    updateCheckboxGroupInput(session, "tablas", selected = union(input$tablas, c("desc", "manova", "anova")))
  })
  observeEvent(input$btn_boot, {
    mostrar_analisis()
    updateCheckboxInput(session, "boot_on", value = TRUE)
    runjs("document.getElementById('sec_boot').open = true;")
  })

  # ---- Lectura de datos ----
  leer_csv <- function(ruta) {
    l1 <- readLines(ruta, n = 1, warn = FALSE, encoding = "UTF-8")
    sep <- if (lengths(regmatches(l1, gregexpr(";", l1))) > lengths(regmatches(l1, gregexpr(",", l1)))) ";" else
      if (grepl("\t", l1)) "\t" else ","
    dec <- if (sep == ";") "," else "."
    d <- utils::read.table(ruta, header = TRUE, sep = sep, dec = dec, check.names = FALSE,
                           stringsAsFactors = FALSE, quote = "\"", comment.char = "",
                           fileEncoding = "UTF-8", na.strings = c("", "NA"))
    d
  }
  limpiar <- function(d) {
    d <- as.data.frame(d, check.names = FALSE)
    names(d) <- make.unique(trimws(names(d)))
    names(d)[!nzchar(names(d))] <- paste0("V", which(!nzchar(names(d))))
    d
  }
  output$hoja_excel_ui <- renderUI({
    req(input$archivo)
    if (grepl("\\.xlsx?$", input$archivo$name, ignore.case = TRUE)) {
      hojas <- readxl::excel_sheets(input$archivo$datapath)
      if (length(hojas) > 1) selectInput("hoja_excel", NULL, hojas, width = "150px")
    }
  })
  observeEvent(list(input$archivo, input$hoja_excel), {
    req(input$archivo)
    f <- input$archivo
    d <- tryCatch({
      if (grepl("\\.xlsx?$", f$name, ignore.case = TRUE)) {
        hoja <- if (is.null(input$hoja_excel)) 1 else input$hoja_excel
        as.data.frame(readxl::read_excel(f$datapath, sheet = hoja))
      } else leer_csv(f$datapath)
    }, error = function(e) {
      showNotification(paste("No se pudo leer el archivo:", conditionMessage(e)), type = "error")
      NULL
    })
    req(d)
    datos(limpiar(d)); nombre_datos(f$name)
  })
  observeEvent(input$ejemplo, {
    req(nzchar(input$ejemplo))
    datos(limpiar(read.csv(file.path("datos", input$ejemplo), check.names = FALSE, encoding = "UTF-8")))
    nombre_datos(input$ejemplo)
    if (input$ejemplo == "cesped_seber.csv") {
      updateRadioButtons(session, "diseno", selected = "dos")
      updateCheckboxInput(session, "log_vars", value = TRUE)
      showNotification("Ejemplo de Amaro et al. (2004): diseño de dos factores con logaritmo de las variables.",
                       duration = 6)
    }
  })

  # ---- Hoja editable ----
  output$hoja <- renderRHandsontable({
    d <- datos()
    rhandsontable(d, rowHeaders = TRUE, stretchH = "none", useTypes = TRUE, height = 600) |>
      hot_cols(columnSorting = TRUE, manualColumnResize = TRUE) |>
      hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
  })
  observeEvent(input$hoja, {
    nuevo <- tryCatch(hot_to_r(input$hoja), error = function(e) NULL)
    req(nuevo)
    nuevo <- as.data.frame(nuevo, check.names = FALSE)
    if (!identical(dim(nuevo), dim(datos())) || !isTRUE(all.equal(nuevo, datos(), check.attributes = FALSE))) {
      datos(nuevo)
    }
  })
  observeEvent(input$btn_nueva_fila, {
    d <- datos(); d[nrow(d) + 1, ] <- NA; datos(d)
  })
  observeEvent(input$btn_nueva_col, {
    showModal(modalDialog(title = "Nueva columna",
      textInput("nueva_col_nombre", "Nombre", paste0("V", ncol(datos()) + 1)),
      radioButtons("nueva_col_tipo", "Tipo", c("Numérica" = "num", "Texto (grupo)" = "txt"), inline = TRUE),
      footer = tagList(modalButton("Cancelar"), actionButton("nueva_col_ok", "Crear"))))
  })
  observeEvent(input$nueva_col_ok, {
    d <- datos(); nm <- make.unique(c(names(d), trimws(input$nueva_col_nombre)))[ncol(d) + 1]
    d[[nm]] <- if (input$nueva_col_tipo == "num") NA_real_ else NA_character_
    datos(d); removeModal()
  })
  output$info_datos <- renderText({
    d <- datos()
    num <- names(d)[vapply(d, is.numeric, logical(1))]
    sprintf("%s  ·  %d filas × %d columnas  ·  %d numéricas  ·  clic derecho para insertar o borrar filas",
            nombre_datos(), nrow(d), ncol(d), length(num))
  })
  output$bajar_datos <- downloadHandler(
    filename = function() paste0(tools::file_path_sans_ext(nombre_datos()), "_editado.csv"),
    content = function(file) write.csv(datos(), file, row.names = FALSE, fileEncoding = "UTF-8"))

  # ---- Selectores de variables ----
  observe({
    d <- datos()
    num <- names(d)[vapply(d, is.numeric, logical(1))]
    cat_cand <- names(d)[vapply(d, function(v) !is.numeric(v) || length(unique(stats::na.omit(v))) <= 12, logical(1))]
    sel_v <- intersect(isolate(input$vars), num)
    if (!length(sel_v)) sel_v <- num
    sel_g <- isolate(input$grupo)
    if (is.null(sel_g) || !sel_g %in% cat_cand) sel_g <- cat_cand[!cat_cand %in% num][1]
    if (is.na(sel_g)) sel_g <- cat_cand[1]
    sel_g2 <- isolate(input$grupo2)
    if (is.null(sel_g2) || !sel_g2 %in% cat_cand || identical(sel_g2, sel_g)) {
      sel_g2 <- setdiff(cat_cand[!cat_cand %in% num], sel_g)[1]
      if (is.na(sel_g2)) sel_g2 <- setdiff(cat_cand, sel_g)[1]
    }
    updateSelectizeInput(session, "vars", choices = num, selected = setdiff(sel_v, c(sel_g, sel_g2)))
    updateSelectInput(session, "grupo", choices = cat_cand, selected = sel_g)
    updateSelectInput(session, "grupo2", choices = cat_cand, selected = sel_g2)
  })
  observe({
    v <- input$vars
    updateSelectInput(session, "alr_ref", choices = v, selected = if (isolate(input$alr_ref) %in% v) isolate(input$alr_ref) else v[length(v)])
  })
  # Variables del analisis (pueden cambiar con la transformacion alr)
  vars_modelo <- reactive({
    cb <- ajuste_base(); if (inherits(cb, "error_ajuste")) character() else rownames(cb$cols)
  })
  prev_vars <- reactiveVal(character())
  observe({
    v <- vars_modelo()
    old <- isolate(input$vars_mostrar)
    sel <- if (is.null(old)) v else union(intersect(old, v), setdiff(v, isolate(prev_vars())))
    updateCheckboxGroupInput(session, "vars_mostrar", choices = v, selected = sel)
    prev_vars(v)
  })
  niveles <- reactive({
    cb <- ajuste_base()
    if (!inherits(cb, "error_ajuste")) return(levels(cb$groups))
    req(input$grupo, input$grupo %in% names(datos()))
    levels(droplevels(as.factor(datos()[[input$grupo]])))
  })
  observe({
    g <- niveles()
    updateCheckboxGroupInput(session, "grupos_mostrar", choices = g, selected = g, inline = TRUE)
  })

  # Titulo del grafico segun el diseno (solo si el usuario no lo ha cambiado)
  titulos_auto <- c("Biplot canónico", "Biplot total", "Biplot de filas", "Biplot de columnas", "Biplot de interacción")
  observeEvent(list(input$diseno, input$efecto), {
    if (isolate(input$titulo) %in% titulos_auto) {
      nuevo <- if (identical(input$diseno, "dos")) c(total = "Biplot total", filas = "Biplot de filas",
                                                     columnas = "Biplot de columnas",
                                                     interaccion = "Biplot de interacción")[[input$efecto]] else "Biplot canónico"
      updateTextInput(session, "titulo", value = nuevo)
    }
  }, ignoreInit = TRUE)

  # ---- Datos del modelo (filtrado, transformacion composicional) ----
  datos_modelo <- reactive({
    d <- datos()
    dos <- identical(input$diseno, "dos")
    req(input$grupo, length(input$vars) >= 2, all(input$vars %in% names(d)), input$grupo %in% names(d))
    if (dos) req(input$grupo2, input$grupo2 %in% names(d))
    vars <- setdiff(input$vars, c(input$grupo, if (dos) input$grupo2))
    validate(need(length(vars) >= 2, "Elige al menos dos variables dependientes numéricas."))
    if (dos) validate(need(input$grupo2 != input$grupo, "Los dos factores deben ser distintos."))
    X <- d[, vars, drop = FALSE]
    g <- d[[input$grupo]]
    g2 <- if (dos) d[[input$grupo2]] else NULL
    ok <- !is.na(g) & (if (dos) !is.na(g2) else TRUE)
    if (input$na_action == "omit") ok <- ok & stats::complete.cases(X)
    X <- X[ok, , drop = FALSE]
    if (is.null(rownames(d)) || all(rownames(d) == seq_len(nrow(d)))) {
      idcol <- names(d)[vapply(d, function(v) is.character(v) && !anyDuplicated(v), logical(1))]
      rownames(X) <- if (length(idcol)) d[[idcol[1]]][ok] else which(ok)
    }
    comp <- es_composicional(X)
    if (isTRUE(input$log_vars) && input$composicional == "ninguna") {
      validate(need(all(as.matrix(X) > 0, na.rm = TRUE), "El logaritmo necesita valores positivos en todas las variables."))
      X[] <- lapply(X, log)
    }
    if (input$composicional != "ninguna") {
      validate(need(!anyNA(X), "La transformación log-cociente necesita datos completos: elimina o imputa los NA."),
               need(all(as.matrix(X) >= 0), "La transformación log-cociente necesita valores no negativos."))
      rn <- rownames(X)
      X <- transformar_composicional(X, input$composicional, input$alr_ref)
      rownames(X) <- rn
    }
    list(X = X, g = droplevels(as.factor(g[ok])), g2 = if (dos) droplevels(as.factor(g2[ok])),
         dos = dos, composicional = comp)
  })

  # ---- Ajuste del biplot canonico (uno o dos factores) ----
  ajuste_base <- reactive({
    dm <- datos_modelo()
    sing <- if (isTRUE(input$ginv) || input$composicional == "clr") "ginv" else "stop"
    na <- if (input$na_action == "impute") "impute" else "fail"
    tryCatch(withCallingHandlers({
      cb <- if (dm$dos) {
        bigarx_manova2(dm$X, dm$g, dm$g2, efecto = input$efecto, transform = input$transform,
                       conf = input$conf, region = input$region, ndim = 10, na_action = na,
                       impute_rank = input$impute_rank, singular = sing,
                       nombre_A = input$grupo, nombre_B = input$grupo2)
      } else {
        bigarx_canonical(dm$X, dm$g, transform = input$transform, conf = input$conf, region = input$region,
                         ndim = 10, na_action = na, impute_rank = input$impute_rank, singular = sing)
      }
      cb$.composicional <- dm$composicional
      cb
    },
      warning = function(w) {
        showNotification(conditionMessage(w), type = "warning", duration = 8)
        invokeRestart("muffleWarning")
      }),
      error = function(e) structure(list(msg = conditionMessage(e)), class = "error_ajuste"))
  })
  ejes_invertidos <- reactive(unique(c(if (isTRUE(input$invertir_x)) input$eje_x, if (isTRUE(input$invertir_y)) input$eje_y)))
  ajuste <- reactive({
    cb <- ajuste_base()
    if (inherits(cb, "error_ajuste")) return(cb)
    if (identical(input$calidad_def, "principal")) cb$quality_cols <- cb$quality_cols_principal
    invertir_ejes(cb, ejes_invertidos())
  })

  # ---- Bootstrap ----
  boot_raw <- reactiveVal(NULL)
  boot_disponible <- reactive({
    cb <- ajuste_base()
    !inherits(cb, "error_ajuste") && (is.null(cb$diseno) || identical(cb$diseno$efecto, "total"))
  })
  run_boot <- function() {
    cb <- ajuste_base(); req(!inherits(cb, "error_ajuste"))
    if (!boot_disponible()) {
      showNotification("El bootstrap está disponible para un factor y para el biplot total de dos factores.",
                       type = "warning", duration = 8)
      boot_raw(NULL); return(invisible())
    }
    dm <- datos_modelo()
    g <- if (dm$dos) cb$diseno$celda else dm$g
    sing <- if (isTRUE(input$ginv) || input$composicional == "clr") "ginv" else "stop"
    res <- withProgress(message = "Bootstrap en curso", value = 0, {
      tryCatch(boot_canonical(dm$X, g, B = input$boot_B,
                              transform = input$transform, axes = c(input$eje_x, input$eje_y),
                              conf = if (input$region == "boot") input$conf else input$boot_conf,
                              ci = input$boot_ci, seed = input$boot_seed,
                              na_action = if (input$na_action == "impute") "impute" else "fail",
                              impute_rank = input$impute_rank, singular = sing,
                              progress = function(p) setProgress(p)),
               error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    })
    boot_raw(res)
  }
  observeEvent(input$boot_run, run_boot())
  observeEvent(list(ajuste_base(), input$eje_x, input$eje_y, input$boot_conf), {
    cb <- ajuste_base()
    if (!inherits(cb, "error_ajuste") && identical(input$region, "boot") && input$eje_x != input$eje_y &&
        max(input$eje_x, input$eje_y) <= cb$ndim) run_boot() else boot_raw(NULL)
  })
  boot <- reactive({
    b <- boot_raw(); if (is.null(b)) return(NULL)
    invertir_boot(b, ejes_invertidos())
  })

  # ---- Colores dinamicos ----
  output$colores_grupo_ui <- renderUI({
    req(input$colores_manual)
    g <- niveles(); cols <- colores_paleta(length(g), input$paleta)
    lapply(seq_along(g), function(i) colourpicker::colourInput(paste0("col_g_", i), g[i], cols[i]))
  })
  output$colores_var_ui <- renderUI({
    req(input$colores_por_var)
    v <- vars_modelo(); cols <- colores_paleta(length(v), "Dark 3")
    lapply(seq_along(v), function(i) colourpicker::colourInput(paste0("col_v_", i), v[i], cols[i]))
  })

  opciones <- reactive({
    cb <- ajuste(); req(!inherits(cb, "error_ajuste"))
    g <- levels(cb$groups)
    cols_g <- if (isTRUE(input$colores_manual)) {
      vals <- vapply(seq_along(niveles()), function(i) {
        v <- input[[paste0("col_g_", i)]]; if (is.null(v)) NA_character_ else v
      }, character(1))
      if (anyNA(vals)) NULL else stats::setNames(vals, niveles())[g]
    } else NULL
    cols_v <- if (isTRUE(input$colores_por_var)) {
      vm <- vars_modelo()
      vals <- vapply(seq_along(vm), function(i) {
        v <- input[[paste0("col_v_", i)]]; if (is.null(v)) NA_character_ else v
      }, character(1))
      stats::setNames(vals, vm)[!is.na(vals)]
    } else NULL
    opciones_grafico(cb,
      axes = c(input$eje_x, input$eje_y),
      grupos = if (is.null(input$grupos_mostrar)) character() else input$grupos_mostrar,
      variables = if (is.null(input$vars_mostrar)) character() else input$vars_mostrar,
      min_calidad = input$min_calidad,
      ver_individuos = input$ver_individuos, ver_medias = input$ver_medias,
      ver_regiones = input$ver_regiones, ver_flechas = input$ver_flechas,
      etiq_medias = input$etiq_medias, etiq_variables = input$etiq_variables,
      etiq_individuos = input$etiq_individuos,
      paleta = input$paleta, colores = cols_g,
      ind_tam = input$ind_tam, ind_alfa = input$ind_alfa, ind_forma = as.numeric(input$ind_forma),
      ind_forma_por_grupo = input$ind_forma_por_grupo,
      media_tam = input$media_tam, media_forma = as.numeric(input$media_forma),
      reg_grosor = input$reg_grosor, reg_linea = input$reg_linea, reg_relleno = input$reg_relleno,
      reg_alfa = input$reg_alfa,
      flecha_color = input$flecha_color, flecha_colores = cols_v, flecha_grosor = input$flecha_grosor,
      flecha_punta = input$flecha_punta, flecha_escala = if (input$flecha_escala == "auto") "auto" else "orig",
      flecha_mult = input$flecha_mult,
      etiq_var_tam = input$etiq_var_tam, etiq_media_tam = input$etiq_media_tam,
      etiq_ind_tam = input$etiq_ind_tam,
      etiq_var_color = if (isTRUE(input$etiq_var_mismo)) NULL else input$etiq_var_color,
      evitar_solape = input$evitar_solape && input$modo == "estatico",
      titulo = input$titulo, subtitulo = input$subtitulo, xlab = input$xlab, ylab = input$ylab,
      pie = input$pie, titulo_leyenda = input$titulo_leyenda,
      tema = input$tema, tam_base = input$tam_base, fuente = input$fuente, leyenda = input$leyenda,
      lineas_cero = input$lineas_cero, cuadricula = input$cuadricula, fondo = input$fondo,
      xlim = c(input$xmin, input$xmax), ylim = c(input$ymin, input$ymax),
      boot_medias = input$boot_medias, boot_variables = input$boot_variables,
      boot_nubes = input$boot_nubes, boot_color = input$boot_color)
  })

  grafico <- reactive({
    cb <- ajuste(); req(!inherits(cb, "error_ajuste"))
    validate(need(input$eje_x != input$eje_y, "Elige dos ejes distintos."),
             need(max(input$eje_x, input$eje_y) <= cb$ndim,
                  sprintf("Este análisis tiene %d ejes canónicos.", cb$ndim)))
    grafico_canonico(cb, opciones(), boot = boot())
  })

  output$grafico <- renderPlot(grafico(), width = function() input$ancho, height = function() input$alto,
                               res = 96)
  output$grafico_int <- plotly::renderPlotly({
    pl <- plotly::ggplotly(grafico(), width = input$ancho, height = input$alto)
    vistos <- character()
    for (i in seq_along(pl$x$data)) {
      nm <- pl$x$data[[i]]$name
      if (is.null(nm) || !nzchar(nm)) next
      limpio <- sub("^\\((.*?),[^,]*\\)$", "\\1", nm)
      pl$x$data[[i]]$name <- limpio
      pl$x$data[[i]]$legendgroup <- limpio
      if (limpio %in% vistos) pl$x$data[[i]]$showlegend <- FALSE else vistos <- c(vistos, limpio)
    }
    pl
  })

  # ---- Tablas ----
  tablas_res <- reactive({
    cb <- ajuste(); req(!inherits(cb, "error_ajuste"))
    dm <- datos_modelo()
    ax <- c(input$eje_x, input$eje_y)
    ax <- ax[ax <= cb$ndim]; if (length(ax) < 2) ax <- c(1, 2)
    out <- list()
    Xd <- if (!is.null(cb$imputation)) cb$imputation$X else as.matrix(dm$X)
    gdesc <- if (is.null(cb$diseno)) cb$groups else cb$diseno$celda
    out$desc <- descriptivos_grupo(Xd, gdesc)
    if (is.null(cb$diseno)) {
      mt <- manova_tests(cb)
      m <- cb$manova
      out$manova <- data.frame(Contraste = rownames(mt), Valor = mt$estadistico, F = mt$F_aprox,
                               `gl num.` = mt$gl1, `gl den.` = mt$gl2, p = mt$p_valor, check.names = FALSE)
      out$bartlett <- data.frame(Contraste = "Wilks (aprox. χ² de Bartlett)", Lambda = m$wilks, `χ²` = m$chi2,
                                 gl = m$df, p = m$p_value, check.names = FALSE)
      ua <- univariate_anovas(cb)
      out$anova <- data.frame(Variable = rownames(ua), `SC entre` = ua$SC_entre, `SC dentro` = ua$SC_dentro,
                              F = ua$F, `gl num.` = ua$gl1, `gl den.` = ua$gl2, p = ua$p_valor, check.names = FALSE)
    } else {
      pr <- cb$diseno$pruebas
      out$manova <- data.frame(Efecto = pr$efecto, Contraste = pr$contraste, Valor = pr$estadistico, F = pr$F_aprox,
                               `gl num.` = pr$gl1, `gl den.` = pr$gl2, p = pr$p_valor, check.names = FALSE)
      w <- pr[pr$contraste == "Wilks", ]
      out$bartlett <- data.frame(Efecto = w$efecto, Lambda = w$estadistico, `χ²` = w$chi2_bartlett,
                                 gl = w$gl_chi2, p = w$p_chi2, check.names = FALSE)
      an <- cb$diseno$anova
      out$anova <- data.frame(Variable = an$variable, Efecto = an$efecto, SC = an$SC, F = an$F,
                              `gl num.` = an$gl1, `gl den.` = an$gl2, p = an$p_valor, check.names = FALSE)
    }
    th <- cb$canonical$valor_propio
    out$eigen <- data.frame(Eje = cb$canonical$eje, `Valor propio (HE⁻¹)` = th,
                            `Valor singular (escala E)` = sqrt(th),
                            `Correlación canónica` = cb$canonical$correlacion_canonica,
                            `Inercia (%)` = cb$inertia$inercia, `Acumulada (%)` = cb$inertia$acumulada,
                            check.names = FALSE)
    if (isTRUE(input$inercia_amaro)) out$eigen$`Inercia θ² (%)` <- 100 * th^2 / sum(th^2)
    out$medias <- data.frame(Grupo = rownames(cb$rows), n = cb$n_groups,
                             cb$rows[, ax, drop = FALSE],
                             `Calidad (%)` = rowSums(cb$quality_rows[, ax, drop = FALSE]), check.names = FALSE)
    names(out$medias)[3:4] <- paste("Eje", ax)
    out$variables <- data.frame(Variable = rownames(cb$cols), cb$cols[, ax, drop = FALSE],
                                `Calidad (%)` = rowSums(cb$quality_cols[, ax, drop = FALSE]),
                                cb$structure[, ax, drop = FALSE], check.names = FALSE)
    names(out$variables)[c(2, 3, 5, 6)] <- c(paste("Eje", ax), paste("r dentro, Eje", ax))
    tot <- cb$structure_total[, ax, drop = FALSE]
    out$corr_total <- data.frame(Variable = rownames(tot), tot, tot^2 * 100, check.names = FALSE)
    names(out$corr_total)[-1] <- c(paste("r total, Eje", ax), paste("r² × 100, Eje", ax))
    b <- boot()
    if (identical(cb$region, "boot")) {
      if (!is.null(b)) {
        rg <- regiones_boot(b)
        out$regiones <- data.frame(Grupo = rg$grupo, `Centro X` = rg$centro1, `Centro Y` = rg$centro2,
                                   `Var X` = rg$var1, `Var Y` = rg$var2, Cov = rg$cov,
                                   `Semieje mayor` = rg$semieje_mayor, `Semieje menor` = rg$semieje_menor,
                                   `Ángulo (°)` = rg$angulo, `Área` = rg$area, check.names = FALSE)
      }
    } else {
      rg <- canonical_regions(cb, ax)
      out$regiones <- if (cb$region == "ellipse") {
        data.frame(Grupo = rg$grupo, `Centro X` = rg$centro1, `Centro Y` = rg$centro2, `Var X` = rg$var1,
                   `Var Y` = rg$var2, Cov = rg$cov, `Semieje mayor` = rg$semieje_mayor,
                   `Semieje menor` = rg$semieje_menor, `Ángulo (°)` = rg$angulo, `Área` = rg$area,
                   check.names = FALSE)
      } else {
        data.frame(Grupo = rg$grupo, `Centro X` = rg$centro1, `Centro Y` = rg$centro2, Radio = rg$semieje_mayor,
                   check.names = FALSE)
      }
    }
    if (nrow(cb$.M) >= 2) {
      hp <- hotelling_pares(cb, input$ajuste_p)
      out$pares <- data.frame(`Grupo 1` = hp$grupo1, `Grupo 2` = hp$grupo2, `D² Mahalanobis` = hp$D2_Mahalanobis,
                              `T²` = hp$T2, F = hp$F, `gl num.` = hp$gl1, `gl den.` = hp$gl2, p = hp$p_valor,
                              `p ajustado` = hp$p_ajustado, check.names = FALSE)
    }
    if (!is.null(cb$imputation)) {
      v <- cb$imputation$values
      out$imput <- data.frame(Fila = v$fila, Variable = v$columna, `Valor imputado` = v$valor, check.names = FALSE)
    }
    if (!is.null(b)) {
      t <- b$table
      nom <- c(inercia = "Inercia (%)", valor_propio = "Valor propio", correlacion_canonica = "Correlación canónica",
               media = "Media de grupo", variable = "Flecha (variable)", calidad_variable = "Calidad de variable (%)",
               correlacion_variable = "r variable-eje")
      out$boot <- data.frame("Parámetro" = nom[t$tipo], Elemento = t$elemento, Eje = t$eje,
                             "Estimación" = t$estimacion, EE = t$ee, Sesgo = t$sesgo,
                             Inferior = t$inferior, Superior = t$superior, check.names = FALSE)
    }
    out
  })

  output$resultados <- renderUI({
    cb <- ajuste()
    if (inherits(cb, "error_ajuste")) {
      return(div(class = "res-error", icon("triangle-exclamation"), cb$msg))
    }
    tb <- tablas_res()
    bloque <- function(clave, titulo, nota = NULL, digits = 3, formatos = NULL, mostrar = clave) {
      if (!mostrar %in% input$tablas || is.null(tb[[clave]])) return(NULL)
      tabla_apa(tb[[clave]], titulo, nota, digits = digits, formatos = formatos, rownames = FALSE)
    }
    dos <- !is.null(cb$diseno)
    efectos <- c(total = "Biplot total (celdas)", filas = paste("Biplot de filas:", input$grupo),
                 columnas = paste("Biplot de columnas:", input$grupo2), interaccion = "Biplot de interacción")
    titulo <- if (dos) paste("MANOVA biplot de dos factores ·", efectos[[cb$diseno$efecto]]) else
      "Biplot canónico / MANOVA biplot"
    notas <- tagList(
      if (dos && isTRUE(cb$diseno$proyectado)) div(class = "res-aviso", icon("circle-info"),
        "Este efecto tiene un solo grado de libertad (dos niveles) y no tiene plano propio: sus medias se",
        "proyectan sobre el biplot total, como en Amaro et al. (2004, figura 5)."),
      if (isTRUE(cb$.composicional) && input$composicional == "ninguna") div(class = "res-aviso",
        icon("circle-info"), "Las variables suman una constante en cada fila (datos composicionales): conviene",
        "usar la transformación clr o alr en «Modelo»."),
      if (isTRUE(cb$singular)) div(class = "res-aviso", icon("circle-info"),
        sprintf("W es singular: se usó la inversa generalizada (rango efectivo %d de %d variables).",
                cb$.krz$p_eff, nrow(cb$cols))),
      if (identical(cb$region, "boot") && is.null(boot())) div(class = "res-aviso", icon("circle-info"),
        if (boot_disponible()) "Calculando las elipses bootstrap..." else
          "Las elipses bootstrap solo están disponibles para un factor y para el biplot total."))
    tfm <- c(ninguna = "", clr = " · clr", alr = paste0(" · alr (ref.: ", input$alr_ref, ")"))[[input$composicional]]
    if (isTRUE(input$log_vars) && input$composicional == "ninguna") tfm <- " · logaritmo natural"
    region_txt <- c(ellipse = "elipses de Krzanowski", chi2 = "círculos chi²", t = "círculos t",
                    bonferroni = "círculos Bonferroni", boot = "elipses bootstrap")[[cb$region]]
    calidad_txt <- if (identical(input$calidad_def, "principal"))
      "Calidad: coseno² sobre las coordenadas principales (Amaro et al., 2004)." else
      "Calidad: correlación al cuadrado (dentro de los grupos) entre la variable y el plano (Amaro et al., 2008)."
    tagList(
      div(class = "res-seccion",
          tags$h2(titulo),
          tags$p(class = "res-sub", sprintf(
            "%d individuos · %d %s · %d variables · transformación: %s%s%s",
            nrow(cb$individuals), nrow(cb$rows), if (dos && cb$diseno$efecto %in% c("total", "interaccion")) "celdas" else "grupos",
            nrow(cb$cols), if (cb$preprocessing$transform == "center") "centrar" else "estandarizar", tfm,
            if (!is.null(cb$imputation)) sprintf(" · %d datos imputados", nrow(cb$imputation$missing)) else "")),
          notas),
      bloque("desc", if (dos) "Descriptivos por celda" else "Descriptivos por grupo",
             formatos = c("txt", "txt", "int", "num", "num")),
      bloque("manova", if (dos) "Contrastes multivariantes por efecto (MANOVA de dos factores)" else
               "Contrastes multivariantes (MANOVA)",
             nota = "Aproximaciones F de <i>summary.manova</i> de R.",
             formatos = if (dos) c("txt", "txt", "num", "num", "int", "num", "p") else
               c("txt", "num", "num", "int", "num", "p")),
      bloque("bartlett", "Lambda de Wilks con la aproximación χ² de Bartlett",
             nota = "χ² = −[gl<sub>E</sub> − (p − gl<sub>H</sub> + 1)/2]·ln Λ.", mostrar = "manova",
             formatos = c("txt", "num", "num", "int", "p"), digits = 4),
      bloque("anova", if (dos) "ANOVA de dos factores por variable" else "ANOVA de una vía por variable",
             formatos = if (dos) c("txt", "txt", "num", "num", "int", "int", "p") else
               c("txt", "num", "num", "num", "int", "int", "p")),
      bloque("eigen", "Valores propios, correlaciones canónicas e inercia",
             nota = paste("El «valor singular (escala E)» es lo que Amaro et al. (2008) llaman valor propio (5,67 en Iris).",
                          if (isTRUE(input$inercia_amaro)) "La inercia θ² reproduce las tablas II, IV y VI de Amaro et al. (2004)." else "")),
      div(class = "res-grafico",
          div(class = "res-titulo", if (dos) efectos[[cb$diseno$efecto]] else "Biplot canónico"),
          if (input$modo == "estatico") plotOutput("grafico", height = "auto") else
            plotly::plotlyOutput("grafico_int", height = paste0(input$alto, "px")),
          if (!is.null(boot()) && !identical(cb$region, "boot")) div(class = "res-nota", sprintf(
            "Líneas discontinuas: elipses bootstrap al %s%% de las medias (%d réplicas estratificadas).",
            format(100 * boot()$conf), boot()$B)),
          if (!is.null(boot()) && identical(cb$region, "boot")) div(class = "res-nota", sprintf(
            "Regiones: elipses bootstrap al %s%% (%d réplicas estratificadas por grupo).",
            format(100 * boot()$conf), boot()$B))),
      bloque("medias", "Medias de grupo: coordenadas y calidad de representación"),
      bloque("variables", "Variables: coordenadas, calidad y correlación con los ejes", nota = calidad_txt),
      bloque("corr_total", "Correlaciones totales entre las variables y los ejes",
             nota = "Correlación sobre todos los individuos (sin quitar el efecto de grupo), como en la tabla 9 de Amaro et al. (2008)."),
      bloque("regiones", sprintf("Regiones de confianza al %s%% (%s)", format(100 * cb$conf), region_txt), digits = 4),
      bloque("pares", "Comparaciones por pares: T² de Hotelling",
             nota = sprintf("W combinada del modelo (gl = %d). Ajuste de p: %s.", cb$.krz$dfE,
                            c(bonferroni = "Bonferroni", holm = "Holm", none = "ninguno")[[input$ajuste_p]]),
             formatos = c("txt", "txt", "num", "num", "num", "int", "int", "p", "p")),
      bloque("imput", "Valores imputados (SVD iterativa)"),
      if (!is.null(tb$boot)) tabla_apa(tb$boot, sprintf("Intervalos bootstrap al %s%% (%s, %d réplicas)",
                                                        format(100 * boot()$conf),
                                                        if (boot()$ci == "percentile") "percentil" else "normal",
                                                        boot()$B), rownames = FALSE)
      else if (isTRUE(input$boot_on)) div(class = "res-nota", "Pulsa «Ejecutar bootstrap» en las opciones."),
      div(class = "res-ref", HTML(paste(
        "Referencias: Gabriel (1972); Amaro, Vicente-Villardón y Galindo (2004, 2008);",
        "Krzanowski (1989); Aitchison (1986); Nieto, Galindo, Leiva y Vicente-Galindo (2015). Motor de cálculo: BiGaRx.")))
    )
  })

  # ---- Exportacion ----
  output$bajar_grafico <- downloadHandler(
    filename = function() paste0("biplot_canonico.", input$exp_formato),
    content = function(file) {
      dev <- switch(input$exp_formato, svg = svglite::svglite, pdf = grDevices::cairo_pdf, "png")
      ggsave(file, grafico(), device = dev, width = input$exp_ancho, height = input$exp_alto,
             units = "cm", dpi = input$exp_dpi, bg = input$fondo)
    })
  output$bajar_tablas <- downloadHandler(
    filename = function() "resultados_biplot_canonico.xlsx",
    content = function(file) writexl::write_xlsx(tablas_res(), file))
  output$bajar_coords <- downloadHandler(
    filename = function() "coordenadas_biplot_canonico.csv",
    content = function(file) {
      cb <- ajuste(); req(!inherits(cb, "error_ajuste"))
      write.csv(canonical_tidy(cb), file, row.names = FALSE, fileEncoding = "UTF-8")
    })
}

shinyApp(ui, server)
