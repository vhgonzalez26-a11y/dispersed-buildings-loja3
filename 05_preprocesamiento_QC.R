# ============================================================
# 05_preprocesamiento_QC.R
# Armonización espacial (CRS, grilla, recorte, máscara), vectores de
# referencia y control de calidad con decisión humana por banda.
# Salidas : data/processed/stack.tif, data/processed/vectors/*.gpkg
#           outputs/decisions/05_qc_decisions.csv
# Puerta  : cada banda ACCEPT o REJECT (REVIEW detiene el flujo)
# ============================================================

if (!exists("ROOT")) {
  ROOT <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), error = function(e) NULL)
  if (is.null(ROOT) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ROOT <- dirname(rstudioapi::getActiveDocumentContext()$path)
  }
  if (is.null(ROOT) || !nzchar(ROOT)) ROOT <- getwd()
}
source(file.path(ROOT, "00_setup.R"), encoding = "UTF-8")

cat("\n== ETAPA 05: PREPROCESAMIENTO Y CONTROL DE CALIDAD ==\n\n")
require_stage("04_datos", "05_preprocesamiento")
assert_protocol_unchanged(CARD)
require_pkgs(c("terra", "sf"))

PREP_CONFIG <- list(
  force_reprocess = FALSE,
  na_flag_pct = 30,           # % de NA dentro del AOI que activa alerta
  qc_file = file.path(DIR$decisions, "05_qc_decisions.csv"),
  stack_file = file.path(DIR$processed, "stack.tif"),
  quicklook_dir = file.path(DIR$figures, "quicklooks")
)
ensure_dir(PREP_CONFIG$quicklook_dir)

project_crs <- card_get(CARD, "project.project_crs")
prov <- read_csv_safe(file.path(DIR$decisions, "04_data_provenance.csv")) |>
  dplyr::mutate(role = lower_clean(role), data_nature = lower_clean(data_nature)) |>
  dplyr::filter(role != "exclude")
abs_path <- function(f) file.path(ROOT, f)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(basename(x)))
facts <- character()

# ------------------------------------------------------------
# 1. Área de estudio (AOI)
# ------------------------------------------------------------

# ADAPTACIÓN A1 (21-oct-2026): un GeoPackage con varias capas necesita indicar cuál se usa. En la columna `notes` de
# 04_data_provenance.csv escriba layer=<nombre> (p. ej. layer=aoi_amplio). Sin esa marca se lee la primera capa, como antes.
layer_of <- function(notes) {
  m <- stringr::str_match(clean_chr(notes), "layer=([^;|,]+)")[, 2]
  ifelse(is.na(m), NA_character_, stringr::str_squish(m))
}
read_vec <- function(f, notes = "") {
  lyr <- layer_of(notes)
  if (is.na(lyr)) sf::st_read(abs_path(f), quiet = TRUE) else sf::st_read(abs_path(f), layer = lyr, quiet = TRUE)
}

aoi_rows <- prov |> dplyr::filter(role == "aoi", data_type == "vector")
if (nrow(aoi_rows) == 0L) human_gate("Asigne role = aoi a un vector del área de estudio en 04_data_provenance.csv.", "05_preprocesamiento")

aoi <- do.call(rbind, lapply(seq_len(nrow(aoi_rows)), function(i) read_vec(aoi_rows$file[i], aoi_rows$notes[i])[, 0])) |>
  sf::st_transform(project_crs) |> sf::st_make_valid() |> sf::st_union() |> sf::st_as_sf()
aoi_file <- file.path(DIR$processed_vectors, "aoi.gpkg")
sf::st_write(aoi, aoi_file, delete_dsn = TRUE, quiet = TRUE)
aoi_v <- terra::vect(aoi)
aoi_area_km2 <- as.numeric(sf::st_area(aoi)) / 1e6
facts <- c(facts, sprintf("Study area polygon reprojected to %s; area = %.2f km2.", project_crs, aoi_area_km2))

# ------------------------------------------------------------
# 2. Vectores de referencia (reproyección y conteo dentro del AOI)
# ------------------------------------------------------------

vector_roles <- c(reference_training = "training", reference_validation = "validation",
                  checkpoints = "checkpoints", field_measurements = "field_measurements", auxiliary = "auxiliary")
vector_qc <- list()
for (role_name in names(vector_roles)) {
  rows <- prov |> dplyr::filter(role == role_name, data_type == "vector")
  if (nrow(rows) == 0L) next
  layers <- lapply(seq_len(nrow(rows)), function(i) sf::st_transform(read_vec(rows$file[i], rows$notes[i]), project_crs))
  count_outside <- function(v) tryCatch(sum(lengths(sf::st_intersects(v, aoi)) == 0L), error = function(e) NA_integer_)
  if (role_name == "auxiliary") {
    # ADAPTACIÓN A1: las capas auxiliares tienen contenido distinto (límites, máscaras, vías); no se fusionan, solo se controla su calidad.
    qc_aux <- purrr::map_dfr(seq_along(layers), function(i) tibble::tibble(
      role = role_name, output = "(no merged: quality control only)", n_features = nrow(layers[[i]]), n_outside_aoi = count_outside(layers[[i]]),
      fields = paste(setdiff(names(layers[[i]]), attr(layers[[i]], "sf_column")), collapse = "; "), file = rows$file[i]))
    vector_qc[[role_name]] <- qc_aux
    facts <- c(facts, sprintf("auxiliary: %d layer(s) checked against the study area (%s).", nrow(rows),
                              paste(sprintf("%s: %d features, %s outside", basename(rows$file), qc_aux$n_features, qc_aux$n_outside_aoi), collapse = "; ")))
    next
  }
  common <- Reduce(intersect, lapply(layers, names))
  v <- do.call(rbind, lapply(layers, function(x) x[, common]))
  out <- file.path(DIR$processed_vectors, paste0(vector_roles[[role_name]], ".gpkg"))
  sf::st_write(v, out, delete_dsn = TRUE, quiet = TRUE)
  n_out <- count_outside(v)
  vector_qc[[role_name]] <- tibble::tibble(role = role_name, output = basename(out), n_features = nrow(v),
                                           n_outside_aoi = n_out, fields = paste(setdiff(common, attr(v, "sf_column")), collapse = "; "),
                                           file = paste(rows$file, collapse = "; "))
  facts <- c(facts, sprintf("%s: %d features from %d file(s), %s outside the study area.", role_name, nrow(v), nrow(rows), n_out))
}
vector_qc <- dplyr::bind_rows(vector_qc)
if (nrow(vector_qc) > 0L) write_csv_safe(vector_qc, file.path(DIR$tables, "05_vector_qc.csv"))

# ------------------------------------------------------------
# 3. Raster predictores: recorte nativo -> reproyección a grilla común -> máscara
# ------------------------------------------------------------

predictors <- prov |> dplyr::filter(role == "predictor", data_type == "raster")
qc_tbl <- tibble::tibble()

if (nrow(predictors) > 0L) {
  ref_row <- predictors |> dplyr::filter(upper_clean(is_reference_grid) == "YES")
  if (nrow(ref_row) == 0L) ref_row <- predictors[1, ]
  resample_method <- function(nature) if (identical(nature, "categorical")) "near" else "bilinear"

  crop_native <- function(r) {
    aoi_native <- terra::project(aoi_v, terra::crs(r))
    terra::crop(r, aoi_native, snap = "out")
  }

  r_ref <- crop_native(terra::rast(abs_path(ref_row$file[1])))
  template <- terra::project(r_ref, project_crs, method = resample_method(ref_row$data_nature[1]))
  template <- terra::crop(template, aoi_v, snap = "out")
  template <- terra::rast(template)  # solo geometría
  facts <- c(facts, sprintf("Reference grid: %s, resolution %.3f x %.3f (map units), %s.",
                            basename(ref_row$file[1]), terra::res(template)[1], terra::res(template)[2], project_crs))

  aoi_cells <- terra::global(terra::rasterize(aoi_v, template, field = 1), "notNA")[1, 1]

  for (i in seq_len(nrow(predictors))) {
    f <- predictors$file[i]
    out_file <- file.path(DIR$processed_layers, paste0(safe_name(f), ".tif"))
    method <- resample_method(predictors$data_nature[i])
    if (!file.exists(out_file) || PREP_CONFIG$force_reprocess) {
      cat("Procesando ", basename(f), " (", method, ")...\n", sep = "")
      r <- crop_native(terra::rast(abs_path(f)))
      band_names <- names(r)
      r <- terra::project(r, template, method = method)
      r <- terra::mask(r, aoi_v)
      names(r) <- paste0(safe_name(f), "__", gsub("[^A-Za-z0-9]+", "_", band_names))
      terra::writeRaster(r, out_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
    }
    r <- terra::rast(out_file)
    stats <- terra::global(r, c("min", "max", "mean", "sd"), na.rm = TRUE)
    not_na <- terra::global(r, "notNA")[, 1]
    ql <- file.path(PREP_CONFIG$quicklook_dir, paste0(safe_name(f), ".png"))
    grDevices::png(ql, width = 1600, height = 1200, res = 200)
    terra::plot(r[[1]], main = names(r)[1]); terra::lines(aoi_v)
    grDevices::dev.off()
    qc_tbl <- dplyr::bind_rows(qc_tbl, tibble::tibble(
      band = names(r), layer_file = basename(out_file), source_file = f,
      resampling = method, na_pct_in_aoi = round(100 * (1 - not_na / aoi_cells), 2),
      min = stats$min, max = stats$max, mean = stats$mean, sd = stats$sd, quicklook = basename(ql)
    ))
  }
  facts <- c(facts, "Predictor rasters were cropped in their native CRS, projected onto the reference grid (bilinear for continuous, nearest neighbour for categorical layers) and masked to the study area.")

  qc_tbl <- qc_tbl |> dplyr::mutate(
    auto_flag = dplyr::case_when(
      na_pct_in_aoi > PREP_CONFIG$na_flag_pct ~ "CHECK: high NA share",
      is.finite(sd) & sd == 0 ~ "CHECK: constant band",
      TRUE ~ "OK"
    )
  )
  old_qc <- if (file.exists(PREP_CONFIG$qc_file)) read_csv_safe(PREP_CONFIG$qc_file) else NULL
  qc <- merge_preserving_human(qc_tbl |> dplyr::mutate(dplyr::across(dplyr::everything(), as.character)),
                               old_qc, "band", c("qc_decision", "qc_note"))
  write_csv_safe(qc, PREP_CONFIG$qc_file)

  decision <- upper_clean(qc$qc_decision)
  if (any(!decision %in% c("ACCEPT", "REJECT"))) {
    human_gate(paste0(
      "Revise los quicklooks en outputs/figures/quicklooks y complete qc_decision (ACCEPT/REJECT)\n",
      "en outputs/decisions/05_qc_decisions.csv. Pendientes o REVIEW: ", sum(!decision %in% c("ACCEPT", "REJECT")), "\n",
      "Bandas con alerta automática: ", sum(qc$auto_flag != "OK"), ". Toda REJECT debe tener qc_note."
    ), "05_preprocesamiento")
  }
  if (any(decision == "REJECT" & !nzchar(qc$qc_note))) human_gate("Cada banda REJECT requiere qc_note.", "05_preprocesamiento")

  accepted <- qc$band[decision == "ACCEPT"]
  if (length(accepted) == 0L) human_gate("No hay bandas aceptadas.", "05_preprocesamiento")
  stack <- terra::rast(lapply(unique(qc$layer_file[decision == "ACCEPT"]), function(lf) terra::rast(file.path(DIR$processed_layers, lf))))
  stack <- stack[[accepted]]
  terra::writeRaster(stack, PREP_CONFIG$stack_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
  writeLines(names(stack), file.path(DIR$processed, "stack_band_names.txt"))
  facts <- c(facts, sprintf("Final predictor stack: %d bands accepted, %d rejected after visual and statistical QC.",
                            length(accepted), sum(decision == "REJECT")))
}

pc <- prov |> dplyr::filter(data_type == "pointcloud")
if (nrow(pc) > 0L) facts <- c(facts, sprintf("%d point-cloud file(s) are processed directly by the analysis module.", nrow(pc)))

write_md(file.path(DIR$evidence, "05_preprocessing_pack.md"), "Preprocessing and quality control", list(
  "Processing facts (generated from code, not from AI)" = paste0("- ", facts, collapse = "\n"),
  "Reference vectors" = vector_qc,
  "Band-level QC" = if (nrow(qc_tbl) > 0L) read_csv_safe(PREP_CONFIG$qc_file) |>
    dplyr::select(band, resampling, na_pct_in_aoi, min, max, auto_flag, qc_decision, qc_note) else "No predictor rasters."
))

mark_stage("05_preprocesamiento", "PASSED")
cat("\nEtapa 05 completada. Siguiente: 06_ejecutar_metodo.R\n")
