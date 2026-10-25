# ============================================================
# R_helpers.R
# Funciones compartidas del pipeline de investigación geoespacial
# Hereda la lógica de "puertas humanas" del flujo PRISMA (curso UTPL)
# y la extiende a datos raster, vectoriales, nubes de puntos y redacción.
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.atomic(x) && is.na(x))) y else x
}

# ------------------------------------------------------------
# 1. Paquetes
# ------------------------------------------------------------

require_pkgs <- function(pkgs, install = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L && install) {
    message("Instalando paquetes faltantes: ", paste(missing, collapse = ", "))
    utils::install.packages(missing, repos = "https://cloud.r-project.org")
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  }
  if (length(missing) > 0L) {
    stop("Paquetes no disponibles: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# ------------------------------------------------------------
# 2. Texto
# ------------------------------------------------------------

clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- stringi::stri_enc_toutf8(x, is_unknown_8bit = TRUE)
  x[is.na(x)] <- ""
  x <- stringi::stri_replace_all_regex(
    x, "[\\x{0000}-\\x{0008}\\x{000B}\\x{000C}\\x{000E}-\\x{001F}\\x{007F}-\\x{009F}]", ""
  )
  x <- stringr::str_replace_all(x, "[\\r\\n\\t]+", " ")
  stringr::str_squish(x)
}

upper_clean <- function(x) stringr::str_to_upper(clean_chr(x))
lower_clean <- function(x) stringr::str_to_lower(clean_chr(x))

normalise_key <- function(x) {
  x <- stringi::stri_trans_general(clean_chr(x), "Latin-ASCII")
  x <- stringr::str_replace_all(stringr::str_to_lower(x), "[^a-z0-9]+", " ")
  stringr::str_squish(x)
}

normalise_doi <- function(x) {
  x <- lower_clean(x)
  x <- stringr::str_remove(x, "^https?://(?:dx\\.)?doi\\.org/")
  stringr::str_trim(stringr::str_remove(x, "^doi:\\s*"))
}

safe_numeric <- function(x) suppressWarnings(as.numeric(clean_chr(x)))

PLACEHOLDER_REGEX <- "\\[(COMPLETAR|AUTHOR INPUT REQUIRED|CITATION REQUIRED|CITATION NOT FOUND|VERIFICAR|TODO)[^]]*\\]"

# ------------------------------------------------------------
# 3. Archivos, registros y auditoría
# ------------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("No se pudo crear la carpeta: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

assert_file <- function(path, label = "Archivo requerido") {
  if (length(path) != 1L || !nzchar(path) || !file.exists(path)) {
    stop(label, " no encontrado:\n", path, call. = FALSE)
  }
  invisible(path)
}

assert_nonempty_file <- function(path, label = "Archivo de salida") {
  assert_file(path, label)
  if (is.na(file.info(path)$size) || file.info(path)$size <= 0) stop(label, " está vacío:\n", path, call. = FALSE)
  invisible(path)
}

stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

write_log <- function(lines, path = file.path(DIR$logs, "pipeline_log.txt")) {
  ensure_dir(dirname(path))
  cat(paste0("[", stamp(), "] ", as.character(lines), collapse = "\n"), "\n", file = path, append = TRUE)
  invisible(path)
}

read_csv_safe <- function(path, ...) {
  assert_file(path)
  readr::read_csv(
    path, col_types = readr::cols(.default = readr::col_character()),
    locale = readr::locale(encoding = "UTF-8"), na = character(), trim_ws = TRUE,
    show_col_types = FALSE, progress = FALSE, name_repair = "unique", ...
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(dplyr::across(dplyr::everything(), clean_chr))
}

write_csv_safe <- function(data, path) {
  ensure_dir(dirname(path))
  readr::write_csv(data, path, na = "")
  assert_nonempty_file(path)
  invisible(path)
}

read_xlsx_safe <- function(path, sheet) {
  assert_file(path)
  if (!sheet %in% openxlsx::getSheetNames(path)) {
    stop("Hoja '", sheet, "' no encontrada en ", basename(path), call. = FALSE)
  }
  # na.strings: openxlsx convierte por defecto la cadena "NA" en un valor perdido cuando el libro se guardó con Excel; en los libros de decisiones "NA" es un valor
  # válido (p. ej. «no aplica» en el checklist), así que se desactiva esa conversión.
  out <- openxlsx::read.xlsx(path, sheet = sheet, detectDates = FALSE, check.names = FALSE, na.strings = "<<sin_na>>")
  if (is.null(out)) return(tibble::tibble())
  tibble::as_tibble(out) |> dplyr::mutate(dplyr::across(dplyr::everything(), clean_chr))
}

write_xlsx_sheets <- function(sheets, path, freeze = TRUE) {
  ensure_dir(dirname(path))
  wb <- openxlsx::createWorkbook()
  head_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "bottom")
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, as.data.frame(sheets[[nm]]), headerStyle = head_style)
    if (freeze) openxlsx::freezePane(wb, nm, firstRow = TRUE)
    openxlsx::setColWidths(wb, nm, cols = seq_len(max(1L, ncol(sheets[[nm]]))), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  assert_nonempty_file(path)
  invisible(path)
}

add_missing_columns <- function(data, columns, default = "") {
  for (nm in setdiff(columns, names(data))) data[[nm]] <- default
  data
}

# Combina una tabla nueva con una editada por la persona sin perder sus decisiones.
# Las columnas "humanas" se conservan desde la versión anterior cuando la clave coincide.
merge_preserving_human <- function(new, old, key, human_cols) {
  if (is.null(old) || nrow(old) == 0L) return(add_missing_columns(new, human_cols))
  old <- add_missing_columns(old, c(key, human_cols))
  keep <- old[, unique(c(key, human_cols)), drop = FALSE]
  new <- new[, setdiff(names(new), human_cols), drop = FALSE]
  out <- dplyr::left_join(new, keep, by = key)
  add_missing_columns(out, human_cols) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(human_cols), clean_chr))
}

md5_file <- function(path) unname(tools::md5sum(path))

# ------------------------------------------------------------
# 4. Puertas humanas, ficha del proyecto y estado de etapas
# ------------------------------------------------------------

human_gate <- function(message, gate = NULL) {
  prefix <- if (is.null(gate)) "PUERTA DE DECISIÓN HUMANA" else paste0("PUERTA DE DECISIÓN HUMANA [", gate, "]")
  write_log(paste0(prefix, ": ", gsub("\n", " | ", message)))
  stop(paste0("\n", prefix, "\n", message, "\n"), call. = FALSE)
}

read_card <- function(path = CARD_FILE) {
  assert_file(path, "Ficha del proyecto (project_card.yaml)")
  yaml::read_yaml(path, eval.expr = FALSE)
}

card_get <- function(card, key, default = NULL) {
  x <- card
  for (p in strsplit(key, ".", fixed = TRUE)[[1]]) {
    if (is.null(x) || !is.list(x) || is.null(x[[p]])) return(default)
    x <- x[[p]]
  }
  x
}

field_problem <- function(x) {
  if (is.null(x) || length(x) == 0L) return("vacío")
  v <- clean_chr(unlist(x, use.names = FALSE))
  if (length(v) == 0L || all(!nzchar(v))) return("vacío")
  if (any(grepl(PLACEHOLDER_REGEX, v))) return("contiene un marcador [COMPLETAR]")
  ""
}

card_require <- function(card, keys, gate) {
  probs <- vapply(keys, function(k) field_problem(card_get(card, k)), character(1))
  bad <- probs[nzchar(probs)]
  if (length(bad) > 0L) {
    human_gate(paste0("Complete estos campos en project_card.yaml:\n",
                      paste0("  - ", names(bad), ": ", bad, collapse = "\n")), gate)
  }
  invisible(TRUE)
}

card_status <- function(card, key) upper_clean(card_get(card, key, "DRAFT"))

stage_file <- function() file.path(DIR$logs, "stage_status.csv")

mark_stage <- function(stage, status = "PASSED", note = "") {
  f <- stage_file()
  old <- if (file.exists(f)) read_csv_safe(f) else tibble::tibble(stage = character(), status = character(), updated_at = character(), note = character())
  old <- old[old$stage != stage, , drop = FALSE]
  new <- dplyr::bind_rows(old, tibble::tibble(stage = stage, status = status, updated_at = stamp(), note = note)) |>
    dplyr::arrange(stage)
  write_csv_safe(new, f)
  write_log(paste0("Etapa ", stage, " -> ", status, if (nzchar(note)) paste0(" (", note, ")") else ""))
  invisible(new)
}

stage_info <- function(stage) {
  f <- stage_file()
  if (!file.exists(f)) return(NULL)
  s <- read_csv_safe(f)
  s[s$stage == stage, , drop = FALSE]
}

require_stage <- function(stage, current = NULL) {
  s <- stage_info(stage)
  if (is.null(s) || nrow(s) == 0L || s$status[1] != "PASSED") {
    human_gate(paste0("La etapa previa ", stage, " no está completada. Ejecute su script antes de continuar."), stage)
  }
  if (!is.null(current)) {
    cur <- stage_info(current)
    if (!is.null(cur) && nrow(cur) > 0L && cur$updated_at[1] < s$updated_at[1]) {
      warning("La etapa ", stage, " cambió después de la última ejecución de ", current,
              ". Los resultados posteriores pueden estar desactualizados.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

pipeline_status <- function() {
  stages <- c("03_protocolo", "04_datos", "05_preprocesamiento",
              "06_metodo")
  f <- stage_file()
  s <- if (file.exists(f)) read_csv_safe(f) else tibble::tibble(stage = character(), status = character(), updated_at = character(), note = character())
  out <- tibble::tibble(stage = stages) |> dplyr::left_join(s, by = "stage") |>
    dplyr::mutate(status = ifelse(is.na(status) | !nzchar(status), "PENDIENTE", status))
  print(as.data.frame(out[, c("stage", "status", "updated_at")]), row.names = FALSE)
  invisible(out)
}

# ------------------------------------------------------------
# 5. Protocolo congelado
# ------------------------------------------------------------

protocol_payload <- function(card) {
  module <- clean_chr(card_get(card, "stage_03_hypothesis.method_module", ""))
  list(
    selected_question = card_get(card, "stage_01_topic.selected_question_en"),
    contribution_type = card_get(card, "stage_01_topic.contribution_type"),
    gap_statement = card_get(card, "stage_02_literature.gap_statement"),
    objectives = card_get(card, "stage_03_hypothesis.objective_general"),
    objectives_specific = card_get(card, "stage_03_hypothesis.objectives_specific"),
    hypotheses = card_get(card, "stage_03_hypothesis.hypotheses"),
    study_area = card_get(card, "stage_03_hypothesis.study_area"),
    method_module = module,
    method_params = card_get(card, paste0("method_params.", module)),
    project_crs = card_get(card, "project.project_crs"),
    random_seed = card_get(card, "project.random_seed")
  )
}

hash_object <- function(obj) {
  tmp <- tempfile(fileext = ".yaml")
  yaml::write_yaml(obj, tmp)
  h <- md5_file(tmp)
  unlink(tmp)
  h
}

frozen_protocol_file <- function() file.path(DIR$protocol, "03_protocol_frozen.yaml")

read_frozen_protocol <- function() {
  f <- frozen_protocol_file()
  if (!file.exists(f)) human_gate("No existe un protocolo congelado. Ejecute 03_hipotesis_y_protocolo.R.", "03_protocolo")
  yaml::read_yaml(f, eval.expr = FALSE)
}

assert_protocol_unchanged <- function(card) {
  frozen <- read_frozen_protocol()
  current_hash <- hash_object(protocol_payload(card))
  if (!identical(current_hash, frozen$hash)) {
    human_gate(paste0(
      "La pregunta, hipótesis, área o parámetros del método cambiaron después de congelar el protocolo.\n",
      "Ejecute 03_hipotesis_y_protocolo.R para registrar y justificar la desviación.\n",
      "Esto evita ajustar hipótesis o parámetros después de ver resultados (HARKing / p-hacking)."
    ), "desviación de protocolo")
  }
  invisible(frozen)
}

flatten_named <- function(x, prefix = "") {
  if (!is.list(x)) return(stats::setNames(paste(clean_chr(x), collapse = "; "), prefix))
  if (length(x) == 0L) return(stats::setNames("", prefix))
  nms <- names(x) %||% as.character(seq_along(x))
  nms[!nzchar(nms)] <- as.character(seq_along(x))[!nzchar(nms)]
  unlist(lapply(seq_along(x), function(i) {
    flatten_named(x[[i]], if (nzchar(prefix)) paste0(prefix, ".", nms[i]) else nms[i])
  }))
}

# ------------------------------------------------------------
# 6. Tablas markdown y paquetes de evidencia
# ------------------------------------------------------------

md_escape <- function(x) {
  x <- clean_chr(x)
  x <- stringr::str_replace_all(x, "\\|", "\\\\|")
  stringr::str_replace_all(x, "`", "'")
}

md_table <- function(data, digits = 3L, max_rows = 200L) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (nrow(data) == 0L || ncol(data) == 0L) return("Sin registros.")
  truncated <- nrow(data) > max_rows
  if (truncated) data <- data[seq_len(max_rows), , drop = FALSE]
  data[] <- lapply(data, function(x) {
    if (is.numeric(x)) x <- ifelse(is.na(x), "", format(round(x, digits), trim = TRUE, scientific = FALSE))
    md_escape(x)
  })
  header <- paste0("| ", paste(md_escape(names(data)), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(data)), collapse = " | "), " |")
  rows <- apply(data, 1L, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  out <- paste(c(header, sep, rows), collapse = "\n")
  if (truncated) out <- paste0(out, "\n\n(Tabla truncada a ", max_rows, " filas.)")
  out
}

write_md <- function(path, title, sections = list()) {
  ensure_dir(dirname(path))
  lines <- c(paste0("# ", title), "", paste0("Generado: ", stamp()), "")
  for (nm in names(sections)) {
    value <- sections[[nm]]
    if (is.null(value) || length(value) == 0L) value <- "Sin registros."
    if (is.data.frame(value)) value <- md_table(value)
    lines <- c(lines, paste0("## ", nm), "", paste(as.character(value), collapse = "\n"), "")
  }
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
  assert_nonempty_file(path)
  invisible(path)
}

read_text <- function(path) if (file.exists(path)) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n") else ""



# ------------------------------------------------------------
# 8. Métricas de exactitud
# ------------------------------------------------------------

metric_rows <- function(metric, value, subset, ci_low = NA_real_, ci_high = NA_real_, unit = "", note = "") {
  tibble::tibble(metric = metric, value = as.numeric(value), ci_low = as.numeric(ci_low),
                 ci_high = as.numeric(ci_high), subset = subset, unit = unit, note = note)
}

# Filas = clase del mapa/predicción; columnas = clase de referencia
classification_metrics <- function(obs, pred, classes = NULL, subset = "all") {
  obs <- as.character(obs); pred <- as.character(pred)
  classes <- classes %||% sort(unique(c(obs, pred)))
  cm <- unclass(table(factor(pred, levels = classes), factor(obs, levels = classes)))
  dimnames(cm) <- list(map = classes, reference = classes)
  n <- sum(cm)
  p <- cm / n
  oa <- sum(diag(p))
  pe <- sum(rowSums(p) * colSums(p))
  kappa <- if (isTRUE(all.equal(pe, 1))) NA_real_ else (oa - pe) / (1 - pe)
  # Desacuerdo de cantidad y de asignación (Pontius & Millones, 2011)
  quantity <- sum(abs(rowSums(p) - colSums(p))) / 2
  allocation <- sum(2 * pmin(rowSums(p) - diag(p), colSums(p) - diag(p))) / 2
  ua <- diag(cm) / rowSums(cm); pa <- diag(cm) / colSums(cm)
  f1 <- ifelse((ua + pa) > 0, 2 * ua * pa / (ua + pa), NA_real_)
  se <- sqrt(oa * (1 - oa) / n)
  overall <- dplyr::bind_rows(
    metric_rows("overall_accuracy", oa, subset, max(0, oa - 1.96 * se), min(1, oa + 1.96 * se), "proportion"),
    metric_rows("kappa", kappa, subset, note = "Report with disagreement components; kappa alone is discouraged"),
    metric_rows("quantity_disagreement", quantity, subset, unit = "proportion"),
    metric_rows("allocation_disagreement", allocation, subset, unit = "proportion"),
    metric_rows("macro_f1", mean(f1, na.rm = TRUE), subset, unit = "proportion"),
    metric_rows("n_samples", n, subset, unit = "count")
  )
  per_class <- tibble::tibble(
    subset = subset, class = classes,
    n_reference = as.integer(colSums(cm)), n_map = as.integer(rowSums(cm)),
    users_accuracy = as.numeric(ua), producers_accuracy = as.numeric(pa), f1 = as.numeric(f1)
  )
  list(confusion = cm, overall = overall, per_class = per_class)
}

# Estimación de área ajustada y exactitud (Olofsson et al., 2014)
olofsson_estimates <- function(cm, mapped_area, subset = "all") {
  cm <- as.matrix(cm)
  classes <- rownames(cm)
  A <- as.numeric(mapped_area[classes]); A[is.na(A)] <- 0
  W <- A / sum(A)
  n_i <- rowSums(cm)
  f <- sweep(cm, 1, ifelse(n_i > 0, n_i, NA_real_), "/"); f[is.na(f)] <- 0
  p <- sweep(f, 1, W, "*")
  oa <- sum(diag(p))
  denom <- ifelse(n_i > 1, n_i - 1, NA_real_)
  se_oa <- sqrt(sum(W^2 * diag(f) * (1 - diag(f)) / denom, na.rm = TRUE))
  se_prop <- vapply(seq_along(classes), function(k) {
    sqrt(sum(W^2 * f[, k] * (1 - f[, k]) / denom, na.rm = TRUE))
  }, numeric(1))
  total <- sum(A)
  list(
    overall = metric_rows("area_weighted_overall_accuracy", oa, subset, max(0, oa - 1.96 * se_oa), min(1, oa + 1.96 * se_oa), "proportion"),
    per_class = tibble::tibble(
      subset = subset, class = classes, mapped_area = A,
      adjusted_area = colSums(p) * total, adjusted_area_ci95 = 1.96 * se_prop * total,
      users_accuracy = diag(p) / rowSums(p), producers_accuracy = diag(p) / colSums(p)
    )
  )
}

continuous_metrics <- function(obs, pred, subset = "all", unit = "", R = 1000L, seed = 2026L) {
  ok <- is.finite(obs) & is.finite(pred)
  o <- obs[ok]; p <- pred[ok]
  n <- length(o)
  calc <- function(o, p) {
    e <- p - o
    c(rmse = sqrt(mean(e^2)), mae = mean(abs(e)), mean_error_bias = mean(e),
      nmad = 1.4826 * stats::median(abs(e - stats::median(e))),
      r2 = if (stats::sd(o) > 0 && stats::sd(p) > 0) stats::cor(o, p)^2 else NA_real_,
      le95_abs_error = stats::quantile(abs(e), 0.95, names = FALSE))
  }
  if (n < 3L) return(metric_rows(c("n_samples", "rmse", "mae", "mean_error_bias", "nmad", "r2"), c(n, rep(NA, 5)), subset, note = "n < 3"))
  point <- calc(o, p)
  set.seed(seed)
  boot <- t(replicate(R, { i <- sample.int(n, n, replace = TRUE); calc(o[i], p[i]) }))
  ci <- apply(boot, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  units <- ifelse(names(point) == "r2", "", unit)
  dplyr::bind_rows(
    metric_rows("n_samples", n, subset, unit = "count"),
    metric_rows(names(point), point, subset, ci[1, ], ci[2, ], units, "95% bootstrap percentile CI")
  )
}

# Validación cruzada por bloques espaciales (reduce el optimismo por autocorrelación)
spatial_block_folds <- function(x, y, block_size, k, seed) {
  block_id <- paste(floor(x / block_size), floor(y / block_size), sep = "_")
  blocks <- unique(block_id)
  if (length(blocks) < k) {
    stop("Solo hay ", length(blocks), " bloques espaciales para ", k,
         " folds. Reduzca block_size_m o cv_folds.", call. = FALSE)
  }
  set.seed(seed)
  fold_of_block <- stats::setNames(sample(rep_len(seq_len(k), length(blocks))), blocks)
  unname(fold_of_block[block_id])
}

# ------------------------------------------------------------
# 9. Figuras y registros de productos
# ------------------------------------------------------------

theme_paper <- function(base_size = 9) {
  family <- if (exists("FIG_FONT", inherits = TRUE)) FIG_FONT else "serif"
  ggplot2::theme_minimal(base_family = family, base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

save_fig <- function(plot, id, width_cm = 17, height_cm = 11, dpi = 600, dir = DIR$figures) {
  ensure_dir(dir)
  svg <- file.path(dir, paste0(id, ".svg"))
  png <- file.path(dir, paste0(id, ".png"))
  ggplot2::ggsave(svg, plot, width = width_cm, height = height_cm, units = "cm", device = svglite::svglite, bg = "white")
  ggplot2::ggsave(png, plot, width = width_cm, height = height_cm, units = "cm", dpi = dpi, bg = "white")
  assert_nonempty_file(png, "Figura PNG")
  c(svg = svg, png = png)
}

register_output <- function(registry_file, row) {
  row <- tibble::as_tibble(lapply(row, function(v) clean_chr(v)))
  old <- if (file.exists(registry_file)) read_csv_safe(registry_file) else NULL
  if (!is.null(old) && nrow(old) > 0L) old <- old[old$id != row$id, , drop = FALSE]
  write_csv_safe(dplyr::bind_rows(old, row), registry_file)
}

figure_registry <- function() file.path(DIR$figures, "figure_registry.csv")
table_registry <- function() file.path(DIR$tables, "table_registry.csv")

register_figure <- function(id, paths, caption_draft, source, width_cm, height_cm) {
  register_output(figure_registry(), list(
    id = id, png = paths[["png"]], svg = paths[["svg"]], caption_draft = caption_draft,
    source_script = source, width_cm = width_cm, height_cm = height_cm, created_at = stamp()
  ))
}

register_table <- function(id, path, description, source) {
  register_output(table_registry(), list(
    id = id, csv = path, description = description, source_script = source, created_at = stamp()
  ))
}

# Muestra un raster para ggplot sin cargarlo completo en memoria
raster_to_df <- function(r, max_cells = 250000L) {
  if (terra::ncell(r) > max_cells) r <- terra::spatSample(r, size = max_cells, method = "regular", as.raster = TRUE)
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "value"
  df
}


message("R_helpers.R cargado correctamente.")
