# ============================================================
# R_helpers_publicacion.R
# Ayudas genéricas de la cadena para dejar listo un artículo (no dependen de A1):
#   1. Paleta segura para daltonismo (Okabe-Ito) y símbolos distintos para las figuras.
#   2. Etiquetas de veredicto de hipótesis según se formularan antes o después de ver los resultados.
#   3. Lectura de resolución y tamaño de PNG (para comprobar el requisito de la revista).
#   4. Rastreo de cifras: cada número del texto debe existir en alguna fuente de resultados, ficha o protocolo.
#   5. Generador del documento de material suplementario en Word (figuras y tablas suplementarias del plan 08).
#   6. Entrega final de la etapa 12 (entrega_final): manuscrito en la plantilla de la revista, carta en Word y paquete de repositorio, guiada por el bloque submission de la ficha.
# Se carga desde 00_setup.R. Origen: lecciones de la auditoría editorial de A1 (24-25 oct 2026).
# ============================================================

# ------------------------------------------------------------
# 1. Paleta segura para daltonismo
# ------------------------------------------------------------

# Okabe y Ito (2008): distinguible con las formas habituales de daltonismo. Se usa SIEMPRE junto con símbolos o tipos de línea distintos.
OKABE_ITO <- c(blue = "#0072B2", orange = "#D55E00", green = "#009E73", pink = "#CC79A7", yellow = "#E69F00", sky = "#56B4E9", black = "#000000")
cb_values <- function(n) unname(OKABE_ITO[c("blue", "orange", "green", "pink", "yellow", "sky")][seq_len(min(n, 6L))])
cb_shapes <- function(n) c(16, 17, 15, 18, 8, 4)[seq_len(min(n, 6L))]

# Escalas de color y forma emparejadas para una variable categórica (niveles con nombre)
cb_scales <- function(levels, name = NULL) {
  list(ggplot2::scale_colour_manual(name = name, values = stats::setNames(cb_values(length(levels)), levels)),
       ggplot2::scale_shape_manual(name = name, values = stats::setNames(cb_shapes(length(levels)), levels)))
}

# ------------------------------------------------------------
# 2. Lenguaje de los veredictos de hipótesis
# ------------------------------------------------------------

# Si las hipótesis o sus umbrales se fijaron después de ver resultados (analysis.hypotheses_post_hoc: true en la ficha), no se habla de hipótesis
# «apoyadas»: se habla de criterios cumplidos o no cumplidos, y las estimaciones con sus intervalos son la evidencia.
verdict_labels <- function(card) {
  # DESCRIPTIVE: hipótesis cuyo criterio no mide lo que se afirma (analysis.descriptive_hypotheses en la ficha); se informa la estimación sin veredicto binario
  if (isTRUE(card_get(card, "analysis.hypotheses_post_hoc", FALSE))) c(SUPPORTED = "Criterion met", NOT_SUPPORTED = "Criterion not met", DESCRIPTIVE = "Assessed descriptively (no criterion)")
  else c(SUPPORTED = "Supported", NOT_SUPPORTED = "Not supported", DESCRIPTIVE = "Assessed descriptively (no criterion)")
}

# ------------------------------------------------------------
# 3. PNG: tamaño y resolución
# ------------------------------------------------------------

png_info <- function(path) {
  b <- readBin(path, "raw", n = 4096L)
  be32 <- function(i) sum(as.numeric(as.integer(b[i:(i + 3L)])) * 256^(3:0))
  out <- c(width = NA_real_, height = NA_real_, dpi = NA_real_)
  if (length(b) >= 24L) out[c("width", "height")] <- c(be32(17L), be32(21L))
  cand <- which(b == as.raw(0x70)); cand <- cand[cand + 12L <= length(b)]
  idx <- cand[vapply(cand, function(j) all(b[j + 1:3] == as.raw(c(0x48, 0x59, 0x73))), logical(1))]
  if (length(idx) > 0L) {
    i <- idx[1L]
    if (as.integer(b[i + 12L]) == 1L) out["dpi"] <- round(be32(i + 4L) * 0.0254)
  }
  out
}

# ------------------------------------------------------------
# 4. Rastreo de cifras
# ------------------------------------------------------------

.round_half_up <- function(x, d) sign(x) * floor(abs(x) * 10^d + 0.5) / 10^d
.fmt_d <- function(x, d) formatC(.round_half_up(x, d), format = "f", digits = d)

# Reúne los números de las tablas de resultados (CSV) en un conjunto de claves «decimales|tipo|valor» (valor y porcentaje, a 0-4 decimales)
number_pool <- function(csv_files) {
  keys <- character()
  for (f in csv_files[file.exists(csv_files)]) {
    tb <- tryCatch(utils::read.csv(f, check.names = FALSE, colClasses = "character", fileEncoding = "UTF-8-BOM"), error = function(e) NULL)
    if (is.null(tb)) next
    x <- suppressWarnings(as.numeric(unlist(tb, use.names = FALSE)))
    x <- unique(x[is.finite(x) & abs(x) < 1e12])
    for (d in 0:4) keys <- c(keys, paste0(d, "|v|", .fmt_d(x, d)), paste0(d, "|v|", .fmt_d(-x, d)), paste0(d, "|p|", .fmt_d(x * 100, d)))
  }
  unique(keys)
}

# Números que aparecen escritos en la ficha, el protocolo o los hechos del método (parámetros, umbrales, descripciones)
param_numbers <- function(text_files) {
  out <- character()
  for (f in text_files[file.exists(text_files)]) {
    t <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    out <- c(out, gsub(",", "", unlist(regmatches(t, gregexpr("[0-9][0-9,]*(?:\\.[0-9]+)?", t, perl = TRUE)))))
  }
  unique(out)
}

# Extrae las cifras de un texto (quita citas autor-año, enlaces, DOI, versiones, referencias a tablas/figuras/secciones, fechas y años)
trace_numbers <- function(text, pool, params = character(), external = character()) {
  t <- text
  t <- gsub("\\([^()]*?[0-9]{4}[ab]?(?:; [^()]*?[0-9]{4}[ab]?)*\\)", " ", t, perl = TRUE)
  t <- gsub("https?://\\S+|10\\.[0-9]{4,}/\\S+", " ", t, perl = TRUE)
  t <- gsub("\\b[0-9]+\\.[0-9]+\\.[0-9]+(?:\\.[0-9]+)?\\b", " ", t, perl = TRUE)
  t <- gsub("\\b(Table|Figure)s? S?[0-9]+(?:[–-][0-9]+)?\\b|\\b(Section|Sections) [0-9]+(?:\\.[0-9]+)*\\b|\\bH[0-9](?:[–-]H?[0-9])?\\b|\\bREF[0-9]+\\b", " ", t, perl = TRUE)
  t <- gsub("[0-9]{1,2} (January|February|March|April|May|June|July|August|September|October|November|December) [0-9]{4}|(January|February|March|April|May|June|July|August|September|October|November|December) [0-9]{4}", " ", t, perl = TRUE)
  t <- gsub("\\b1:[0-9][0-9,.]*", " ", t, perl = TRUE)
  loc <- stringr::str_locate_all(t, "[-−]?[0-9][0-9,]*(?:\\.[0-9]+)?%?")[[1]]
  if (nrow(loc) == 0L) return(tibble::tibble(token = character(), status = character(), context = character()))
  rows <- lapply(seq_len(nrow(loc)), function(i) {
    tok <- substr(t, loc[i, 1], loc[i, 2])
    ctx <- stringr::str_squish(substr(t, max(1L, loc[i, 1] - 45L), min(nchar(t), loc[i, 2] + 35L)))
    pct <- grepl("%$", tok)
    base <- gsub("[%,−-]", "", tok)
    d <- if (grepl("\\.", base)) nchar(sub("^.*\\.", "", base)) else 0L
    val <- suppressWarnings(as.numeric(base))
    st <- if (is.na(val)) "NO_DATA"
    else if (d == 0L && !pct && val < 100) "NO_DATA"
    else if (d == 0L && val >= 1900 && val <= 2030 && !grepl(",", tok)) "NO_DATA"
    else if (any(c(paste0(d, "|v|", .fmt_d(val, d)), paste0(d, "|p|", .fmt_d(val, d))) %in% pool)) "TRACED"
    else if (base %in% params || sub("0+$", "", sub("\\.$", "", base)) %in% params) "PARAMETER"
    else if (base %in% external) "EXTERNAL"
    else "UNTRACED"
    tibble::tibble(token = tok, status = st, context = ctx)
  })
  dplyr::bind_rows(rows)
}

# ------------------------------------------------------------
# 5. Material suplementario en Word
# ------------------------------------------------------------

.supp_fmt <- function(v) {
  v <- as.character(v)
  vapply(v, function(s) {
    s <- trimws(s)
    if (is.na(s) || s == "" || s == "NA") return("")
    if (s %in% c("TRUE", "FALSE")) return(ifelse(s == "TRUE", "Yes", "No"))
    x <- suppressWarnings(as.numeric(s))
    if (is.na(x)) return(s)
    if (grepl("^-?[0-9]+$", s)) return(format(as.numeric(s), big.mark = ",", scientific = FALSE))
    if (abs(x) >= 1000) return(formatC(x, format = "f", digits = 1, big.mark = ","))
    if (abs(x) < 0.001 && x != 0) return(formatC(x, format = "e", digits = 2))
    formatC(x, format = "f", digits = 3)
  }, character(1), USE.NAMES = FALSE)
}

.supp_header <- function(x) {
  map <- c(ci_low = "Lower limit", ci_high = "Upper limit", n_presences = "Presences (n)", n_background = "Background points (n)", area_km2 = "Area (km2)")
  out <- ifelse(x %in% names(map), map[x], gsub("\\s+", " ", gsub("_+", " ", x)))
  paste0(toupper(substr(out, 1L, 1L)), substring(out, 2L))
}

# Tabla de hipótesis en el suplemento: criterio, valor observado e intervalo, sin los textos largos de la etapa 07 (el enunciado está en el texto principal)
.supp_hypothesis_table <- function(hyp, labels) {
  desc <- hyp$id %in% as.character(unlist(card_get(CARD, "analysis.descriptive_hypotheses", character())))   # sin criterio: se informa la estimación (ficha)
  tibble::tibble(ID = hyp$id, `Metric (subset)` = paste0(hyp$metric, " (", hyp$subset, ")", ifelse(nzchar(hyp$compare_subset %||% ""), paste0(" vs ", hyp$compare_subset), "")),
                 Criterion = ifelse(desc, "None (descriptive)", paste0(ifelse(nzchar(hyp$compare_subset %||% ""), "Difference ", ""), hyp$operator, " ", .supp_fmt(hyp$threshold))),
                 Observed = .supp_fmt(hyp$observed), `Lower limit` = .supp_fmt(hyp$ci_low), `Upper limit` = .supp_fmt(hyp$ci_high),
                 `Compared value` = .supp_fmt(hyp$compared_value), Difference = .supp_fmt(hyp$difference),
                 Result = ifelse(desc, labels[["DESCRIPTIVE"]], ifelse(hyp$verdict %in% names(labels), labels[hyp$verdict], hyp$verdict)))
}

# Genera Supplementary_Materials.docx con las figuras (PNG) y tablas suplementarias (final_id con «_S») de los manifiestos de la etapa 08
build_supplementary_docx <- function(path, title, authors, fig_manifest, tab_manifest, figures_dir, card, hypothesis_csv = NULL, st = list(title = "Title", heading1 = "heading 1", normal = "Normal"), template = NULL) {
  supp_f <- fig_manifest[grepl("_S", fig_manifest$final_id), , drop = FALSE]
  supp_t <- tab_manifest[grepl("_S", tab_manifest$final_id), , drop = FALSE]
  num_s <- function(x) as.integer(sub("^.*_S([0-9]+).*$", "\\1", x))
  supp_f <- supp_f[order(num_s(supp_f$final_id)), , drop = FALSE]
  supp_t <- supp_t[order(num_s(supp_t$final_id)), , drop = FALSE]
  notes <- card_get(card, "supplementary.table_notes", list())
  intro <- clean_chr(card_get(card, "supplementary.intro",
                              "This file contains the supplementary figures and tables cited in the main text. Intervals reported for cross-validation metrics are ranges across folds and not confidence intervals, unless stated."))
  labels <- verdict_labels(card)
  nice <- function(id) gsub("_", " ", id)
  doc <- if (!is.null(template) && file.exists(template)) officer::read_docx(template) else officer::read_docx()
  avail <- officer::styles_info(doc)$style_name
  ok_style <- function(s) if (s %in% avail) s else "Normal"   # la plantilla puede no traer todos los estilos
  st <- lapply(st, ok_style)
  add <- function(txt, style = st$normal) officer::body_add_par(doc, txt, style = style)
  doc <- add("Supplementary Materials", st$title); doc <- add(title); doc <- add(authors); doc <- add(""); doc <- add(intro); doc <- add("")
  doc <- add("Contents", st$heading1)
  for (i in seq_len(nrow(supp_f))) doc <- add(paste0(nice(supp_f$final_id[i]), ": ", sub("(?<=\\.)\\s.*$", "", clean_chr(supp_f$final_caption[i]), perl = TRUE)))
  for (i in seq_len(nrow(supp_t))) doc <- add(paste0(nice(supp_t$final_id[i]), ": ", clean_chr(supp_t$final_title[i])))
  if (nrow(supp_f) > 0L) {
    doc <- officer::body_add_break(doc)
    doc <- add("Supplementary figures", st$heading1)
    for (i in seq_len(nrow(supp_f))) {
      png <- file.path(figures_dir, paste0(supp_f$final_id[i], ".png"))
      inf <- png_info(png)
      w <- 15.5; h <- min(w * inf[["height"]] / inf[["width"]], 19)
      doc <- officer::body_add_img(doc, src = png, width = h * inf[["width"]] / inf[["height"]] / 2.54, height = h / 2.54)
      doc <- add(paste0(nice(supp_f$final_id[i]), ". ", clean_chr(supp_f$final_caption[i])))
      doc <- officer::body_add_break(doc)
    }
  }
  if (nrow(supp_t) > 0L) {
    doc <- officer::body_end_section_portrait(doc)
    doc <- add("Supplementary tables", st$heading1)
    for (i in seq_len(nrow(supp_t))) {
      id <- supp_t$final_id[i]
      csv <- supp_t$csv[i]
      tb <- read_csv_safe(csv)
      if (!is.null(hypothesis_csv) && basename(csv) == basename(hypothesis_csv)) {
        tb <- .supp_hypothesis_table(tb, labels)
      } else {
        keep <- stringr::str_split(supp_t$columns_to_keep[i], "\\s*;\\s*")[[1]]; keep <- keep[nzchar(keep) & keep %in% names(tb)]
        if (length(keep) > 0L) tb <- tb[, keep, drop = FALSE]
        names(tb) <- .supp_header(names(tb))
        tb[] <- lapply(tb, .supp_fmt)
      }
      doc <- add(paste0(nice(id), ". ", sub("\\.$", "", clean_chr(supp_t$final_title[i])), "."))
      ft <- flextable::flextable(as.data.frame(tb)) |> flextable::fontsize(size = ifelse(ncol(tb) > 9L, 6.5, 7.5), part = "all") |>
        flextable::bold(part = "header") |> flextable::bg(part = "header", bg = "#D9D9D9") |> flextable::autofit()
      doc <- flextable::body_add_flextable(doc, ft)
      note <- notes[[id]] %||% notes[[gsub("_", "", id)]]
      if (!is.null(note) && nzchar(clean_chr(note))) doc <- add(paste0("Note. ", clean_chr(note)))
      doc <- add("")
    }
    doc <- officer::body_end_section_landscape(doc)
  }
  print(doc, target = path)
  invisible(path)
}

# ------------------------------------------------------------
# 6. Entrega final (etapa 12)
# ------------------------------------------------------------
# Ejecuta un script de Python con constancia (run_logged.py). `env` = vector "NOMBRE=valor". Devuelve list(ok, salida).
run_python_step <- function(script, purpose, env = character()) {
  py <- Sys.which("python"); if (!nzchar(py)) stop("No se encontró python en el PATH.", call. = FALSE)
  runner <- normalizePath(file.path(ROOT, "..", "codigo_py", "_herramientas", "run_logged.py"), winslash = "/", mustWork = TRUE)
  script <- normalizePath(file.path(ROOT, "..", "codigo_py", script), winslash = "/", mustWork = TRUE)
  # system2(env=) no funciona en Windows: las variables se fijan con Sys.setenv y se restauran al terminar
  vars <- c("PYTHONUTF8=1", paste0("A1_RAIZ=", ROOT), env)
  nombres <- sub("=.*$", "", vars); previos <- Sys.getenv(nombres, unset = NA)
  do.call(Sys.setenv, as.list(stats::setNames(sub("^[^=]*=", "", vars), nombres)))
  on.exit(for (i in seq_along(nombres)) if (is.na(previos[i])) Sys.unsetenv(nombres[i]) else do.call(Sys.setenv, as.list(stats::setNames(previos[i], nombres[i]))), add = TRUE)
  out <- suppressWarnings(system2(py, c(shQuote(runner), "--no-isolate", shQuote(purpose), shQuote(script)), stdout = TRUE, stderr = TRUE))
  list(ok = is.null(attr(out, "status")) || attr(out, "status") == 0L, salida = out)
}

# Un archivo ya existente en la carpeta final se considera editado por el autor (no se pisa) salvo que su hash sea el que dejó esta función en el registro.
# Devuelve la ruta donde escribir: la misma si es seguro, o una versión «_regenerado_<fecha>» si no.
.destino_protegido <- function(destino, registro) {
  if (!file.exists(destino)) return(list(ruta = destino, protegido = FALSE))
  reg <- if (file.exists(registro)) read_csv_safe(registro) else NULL
  h <- md5_file(destino)
  generado_por_cadena <- !is.null(reg) && any(reg$file == basename(destino) & reg$md5 == h)
  if (generado_por_cadena) return(list(ruta = destino, protegido = FALSE))
  list(ruta = file.path(dirname(destino), paste0(tools::file_path_sans_ext(basename(destino)), "_regenerado_", format(Sys.time(), "%Y%m%d_%H%M"), ".docx")), protegido = TRUE)
}
.registrar_entrega <- function(ruta, registro) {
  reg <- if (file.exists(registro)) read_csv_safe(registro) else tibble::tibble(file = character(), md5 = character(), generated = character())
  reg <- reg |> dplyr::filter(file != basename(ruta)) |> dplyr::bind_rows(tibble::tibble(file = basename(ruta), md5 = md5_file(ruta), generated = format(Sys.time(), "%Y-%m-%d %H:%M")))
  write_csv_safe(reg, registro)
}

# Versión vigente del manuscrito en la carpeta final: el DOCX más reciente que no es la carta (puede ser la editada por el autor o una regenerada).
ultima_version_final <- function(final_dir, carta_name) {
  f <- list.files(final_dir, pattern = "[.]docx$", full.names = TRUE)
  f <- f[!grepl("^~[$]", basename(f)) & !grepl("^(Cover_letter|Supplementary|Revisi)", basename(f)) & basename(f) != carta_name]
  if (length(f) == 0L) NULL else f[which.max(file.mtime(f))]
}

# Documentos que los autores ya revisaron en la carpeta final (p. ej. «…revisadoV4.docx», «Cover_letter_…_revisada.docx»): el más reciente es el vigente y la cadena no genera otro encima.
# Devuelve NULL si no hay ninguno. `carta` = TRUE busca cartas; FALSE, manuscritos (sin cartas, suplementos ni copias regeneradas).
documento_revisado <- function(final_dir, carta = FALSE, generado = NULL) {
  f <- list.files(final_dir, pattern = "[.]docx$", full.names = TRUE)
  f <- f[!grepl("^~[$]", basename(f)) & !grepl("^Supplementary", basename(f)) & !grepl("_regenerado_", basename(f)) & !(basename(f) %in% generado)]
  f <- f[if (carta) grepl("^Cover_letter", basename(f)) else !grepl("^Cover_letter", basename(f))]
  if (!carta) f <- f[!grepl("^Revisi", basename(f))]   # los informes de revisión externa no son el manuscrito
  if (length(f) == 0L) NULL else f[which.max(file.mtime(f))]
}

# Entrega final: 1) manuscrito en la plantilla de la revista (ítem C03 del checklist), 2) carta de presentación en Word, 3) paquete de repositorio (ítem C19).
# C03 o C19 en NA omiten el paso con la nota del checklist. Escribe outputs/decisions/12_entrega_final.csv. Devuelve esa tabla.
entrega_final <- function(checklist, card = CARD) {
  g <- function(k) { v <- card_get(card, paste0("submission.", k)); if (is.null(v) || !nzchar(clean_chr(v))) NULL else clean_chr(v) }
  final_dir <- file.path(ROOT, g("final_dir") %||% "01_articulo"); dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
  registro <- file.path(final_dir, ".registro_entrega.csv")
  carta_name <- g("cover_letter_name") %||% "Cover_letter.docx"
  rows <- list()
  add <- function(paso, estado, detalle) rows[[length(rows) + 1L]] <<- tibble::tibble(step = paso, status = estado, detail = as.character(detalle), when = format(Sys.time(), "%Y-%m-%d %H:%M"))
  item <- function(id) checklist[checklist$item_id == id, , drop = FALSE]
  decide <- function(id) { r <- item(id); if (nrow(r) == 0L) "YES" else if (r$status[1] == "NA") "NA" else "YES" }
  nota <- function(id) { r <- item(id); if (nrow(r) == 0L) "" else r$note[1] }

  # 1. Plantilla de la revista (el borrador sale del export de la etapa 11; el documento final vigente se toma de la carpeta final)
  if (decide("C03") == "NA" || is.null(g("template_script"))) {
    add("template", "SKIPPED", paste0("C03 = NA: ", nota("C03")))
  } else {
    destino <- file.path(final_dir, g("final_manuscript_name") %||% "Manuscript_final_template.docx")
    revisado <- if (isTRUE(card_get(card, "submission.regenerate_if_revised", FALSE))) NULL else documento_revisado(final_dir, FALSE, basename(destino))
    d <- if (!is.null(revisado)) list(ruta = revisado, protegido = FALSE) else .destino_protegido(destino, registro)
    r <- if (!is.null(revisado)) list(ok = TRUE, salida = character()) else run_python_step(g("template_script"), "Etapa 12: pasar el último borrador a la plantilla de la revista", env = paste0("A1_SALIDA=", d$ruta))
    if (!is.null(revisado)) {
      add("template", "OK", paste0("Versión revisada por los autores en la carpeta final: ", basename(revisado), " (no se genera otra; para regenerarla desde los borradores, submission.regenerate_if_revised: true)"))
      v <- run_python_step(g("verify_script"), "Etapa 12: verificar el manuscrito en plantilla", env = paste0("A1_SALIDA=", d$ruta))
      ver <- file.path(DIR$evidence, "verificacion_docx_land.txt")
      fallos <- if (file.exists(ver)) grep("^FAIL", readLines(ver, warn = FALSE, encoding = "UTF-8"), value = TRUE) else "sin informe"
      add("template_verification", if (v$ok && length(fallos) == 0L) "OK" else "FAIL", if (length(fallos)) paste(fallos, collapse = " | ") else paste0(basename(ver), " sobre ", basename(revisado)))
    } else if (!r$ok || !file.exists(d$ruta)) {
      add("template", "FAIL", paste(utils::tail(r$salida, 3), collapse = " | "))
    } else {
      .registrar_entrega(d$ruta, registro)
      add("template", if (d$protegido) "WARN" else "OK",
          if (d$protegido) paste0("El documento final ya tenía cambios del autor y no se pisó. Versión nueva: ", basename(d$ruta), ". Compárela con la anterior y traspase lo que corresponda.") else basename(d$ruta))
      if (!is.null(g("verify_script"))) {
        v <- run_python_step(g("verify_script"), "Etapa 12: verificar el manuscrito en plantilla", env = paste0("A1_SALIDA=", d$ruta))
        ver <- file.path(DIR$evidence, "verificacion_docx_land.txt")
        fallos <- if (file.exists(ver)) grep("^FAIL", readLines(ver, warn = FALSE, encoding = "UTF-8"), value = TRUE) else "sin informe"
        add("template_verification", if (v$ok && length(fallos) == 0L) "OK" else "FAIL", if (length(fallos)) paste(fallos, collapse = " | ") else basename(ver))
      }
    }
    vig <- ultima_version_final(final_dir, carta_name)
    add("final_document_current", "INFO", if (is.null(vig)) "sin documento" else paste0(basename(vig), " (", format(file.mtime(vig), "%Y-%m-%d %H:%M"), ")"))
  }

  # 2. Carta de presentación en Word (para revisar y firmar); el borrador 08 ya se comprobó sin marcadores en el bloque 6.5
  if (is.null(g("cover_letter_script"))) {
    add("cover_letter", "SKIPPED", "sin cover_letter_script en la ficha")
  } else {
    carta_rev <- if (isTRUE(card_get(card, "submission.regenerate_if_revised", FALSE))) NULL else documento_revisado(final_dir, TRUE, carta_name)
    d <- .destino_protegido(file.path(final_dir, carta_name), registro)
    r <- if (!is.null(carta_rev)) list(ok = TRUE, salida = character()) else run_python_step(g("cover_letter_script"), "Etapa 12: carta de presentación en Word", env = paste0("A1_CARTA=", d$ruta))
    if (!is.null(carta_rev)) {
      add("cover_letter", "OK", paste0("Carta revisada por los autores en la carpeta final: ", basename(carta_rev), " (no se genera otra)"))
    } else if (!r$ok || !file.exists(d$ruta)) {
      add("cover_letter", "FAIL", paste(utils::tail(r$salida, 3), collapse = " | "))
    } else {
      .registrar_entrega(d$ruta, registro)
      add("cover_letter", if (d$protegido) "WARN" else "OK", if (d$protegido) paste0("La carta ya tenía cambios del autor y no se pisó: ", basename(d$ruta)) else basename(d$ruta))
    }
  }

  # 3. Paquete de repositorio (GitHub y Zenodo): solo si el checklist (C19) dice que hace falta
  if (decide("C19") == "NA" || is.null(g("repository_script"))) {
    add("repository_package", "SKIPPED", paste0("C19 = NA: ", nota("C19")))
  } else {
    r <- run_python_step(g("repository_script"), "Etapa 12: paquete del repositorio (GitHub y Zenodo)")
    add("repository_package", if (r$ok) "OK" else "FAIL", paste(utils::tail(r$salida[nzchar(r$salida)], if (r$ok) 2L else 3L), collapse = " | "))
  }
  res <- dplyr::bind_rows(rows)
  write_csv_safe(res, file.path(DIR$decisions, "12_entrega_final.csv"))
  res
}
