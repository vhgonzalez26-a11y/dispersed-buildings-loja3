# ============================================================
# 04_inventario_datos.R
# Inventario automático de raster, vectores y nubes de puntos +
# registro de procedencia (fuente, fecha, licencia, rol).
# Salidas : outputs/tables/04_data_inventory.csv
#           outputs/decisions/04_data_provenance.csv  (lo completa usted)
# Puerta  : procedencia completa, CRS definido, un solo raster de referencia
# ============================================================

if (!exists("ROOT")) {
  ROOT <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), error = function(e) NULL)
  if (is.null(ROOT) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ROOT <- dirname(rstudioapi::getActiveDocumentContext()$path)
  }
  if (is.null(ROOT) || !nzchar(ROOT)) ROOT <- getwd()
}
source(file.path(ROOT, "00_setup.R"), encoding = "UTF-8")

cat("\n== ETAPA 04: INVENTARIO Y PROCEDENCIA DE DATOS ==\n\n")
require_stage("03_protocolo", "04_datos")
assert_protocol_unchanged(CARD)

INV_CONFIG <- list(
  hash_files_smaller_than_mb = 2048,   # MD5 de archivos originales para trazabilidad
  provenance_file = file.path(DIR$decisions, "04_data_provenance.csv")
)

EXT <- list(
  raster = "\\.(tif|tiff|img|vrt|jp2|nc)$",
  vector = "\\.(gpkg|shp|geojson|kml|fgb)$",
  pointcloud = "\\.(las|laz)$",
  tabular = "\\.(csv|xlsx)$"
)

files <- list.files(DIR$raw, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
files <- files[!grepl("LEEME_datos\\.txt$|\\.aux\\.xml$|\\.(dbf|shx|prj|cpg|sbn|sbx|ovr)$", files, ignore.case = TRUE)]
if (length(files) == 0L) human_gate("data/raw está vacío. Copie sus datos originales en raster/, vector/, pointcloud/ o tabular/.", "04_datos")

type_of <- function(f) {
  for (t in names(EXT)) if (grepl(EXT[[t]], f, ignore.case = TRUE)) return(t)
  "other"
}
types <- vapply(files, type_of, character(1))
if (any(types %in% c("raster", "vector"))) require_pkgs(c("terra", "sf"))
if (any(types == "pointcloud")) require_pkgs(c("lidR", "sf"))

rel <- function(f) sub(paste0("^", ROOT, "/?"), "", normalizePath(f, winslash = "/"))
project_crs <- card_get(CARD, "project.project_crs")

# ------------------------------------------------------------
# 1. Metadatos por tipo
# ------------------------------------------------------------

meta_raster <- function(f) {
  r <- terra::rast(f)
  e <- as.vector(terra::ext(r))
  list(crs = terra::crs(r, describe = TRUE)$code %||% "", crs_name = terra::crs(r, describe = TRUE)$name %||% "",
       n_layers = terra::nlyr(r), res_x = terra::res(r)[1], res_y = terra::res(r)[2],
       n_cells = terra::ncell(r), extent = paste(round(e, 3), collapse = ", "),
       layer_names = paste(names(r), collapse = "; "), datatype = paste(unique(terra::datatype(r)), collapse = "; "))
}

meta_vector <- function(f) {
  v <- sf::st_read(f, quiet = TRUE)
  crs <- sf::st_crs(v)
  bb <- sf::st_bbox(v)
  list(crs = if (is.na(crs$epsg)) "" else as.character(crs$epsg), crs_name = crs$Name %||% "",
       n_features = nrow(v), geometry = paste(unique(as.character(sf::st_geometry_type(v))), collapse = "; "),
       extent = paste(round(as.numeric(bb), 3), collapse = ", "),
       fields = paste(setdiff(names(v), attr(v, "sf_column")), collapse = "; "))
}

meta_pointcloud <- function(f) {
  h <- lidR::readLASheader(f)
  hd <- h@PHB
  n <- hd[["Number of point records"]]
  area <- (hd[["Max X"]] - hd[["Min X"]]) * (hd[["Max Y"]] - hd[["Min Y"]])
  epsg <- tryCatch(sf::st_crs(h)$epsg, error = function(e) NA)
  list(crs = if (isTRUE(epsg > 0)) as.character(epsg) else "", n_points = n,
       bbox_density_pts_m2 = if (area > 0) round(n / area, 2) else NA,
       extent = paste(round(c(hd[["Min X"]], hd[["Min Y"]], hd[["Max X"]], hd[["Max Y"]]), 3), collapse = ", "),
       las_version = paste0(hd[["Version Major"]], ".", hd[["Version Minor"]]))
}

inventory <- purrr::map_dfr(seq_along(files), function(i) {
  f <- files[i]; t <- types[i]
  info <- tryCatch(switch(t,
    raster = meta_raster(f), vector = meta_vector(f), pointcloud = meta_pointcloud(f),
    list()), error = function(e) list(read_error = conditionMessage(e)))
  size_mb <- file.info(f)$size / 1024^2
  md5 <- if (size_mb < INV_CONFIG$hash_files_smaller_than_mb) md5_file(f) else "not_hashed_large_file"
  tibble::as_tibble(c(list(file = rel(f), data_type = t, size_mb = round(size_mb, 2), md5 = md5),
                      lapply(info, function(x) clean_chr(paste(x, collapse = "; ")))))
})
inventory <- add_missing_columns(inventory, c("crs", "read_error"))
write_csv_safe(inventory, file.path(DIR$tables, "04_data_inventory.csv"))
register_table("tab_data_inventory", file.path(DIR$tables, "04_data_inventory.csv"), "Automatic inventory of raw datasets", "04")

# ------------------------------------------------------------
# 2. Plantilla de procedencia (se conservan sus respuestas)
# ------------------------------------------------------------

PROV_HUMAN <- c("role", "data_nature", "is_reference_grid", "source", "platform_sensor",
                "acquisition_date", "processing_level", "license", "citation_or_doi", "notes")

old_prov <- if (file.exists(INV_CONFIG$provenance_file)) read_csv_safe(INV_CONFIG$provenance_file) else NULL
prov <- merge_preserving_human(inventory |> dplyr::select(file, data_type, crs), old_prov, "file", PROV_HUMAN)
write_csv_safe(prov, INV_CONFIG$provenance_file)

# ------------------------------------------------------------
# 3. Puertas
# ------------------------------------------------------------

ROLES <- c("predictor", "reference_training", "reference_validation", "aoi", "pointcloud",
           "checkpoints", "field_measurements", "auxiliary", "exclude")
problems <- character()

read_err <- inventory |> dplyr::filter(nzchar(read_error))
if (nrow(read_err) > 0L) problems <- c(problems, paste0("No se pudo leer: ", read_err$file, " (", read_err$read_error, ")"))

active <- prov |> dplyr::filter(lower_clean(role) != "exclude")
bad_role <- prov |> dplyr::filter(!lower_clean(role) %in% ROLES)
if (nrow(bad_role) > 0L) problems <- c(problems, paste0("role vacío o inválido en: ", bad_role$file))

for (col in c("source", "acquisition_date", "license")) {
  miss <- active |> dplyr::filter(!nzchar(.data[[col]]))
  if (nrow(miss) > 0L) problems <- c(problems, paste0(col, " vacío en: ", miss$file))
}

spatial_active <- active |> dplyr::filter(data_type %in% c("raster", "vector", "pointcloud"))
no_crs <- spatial_active |> dplyr::filter(!nzchar(crs))
if (nrow(no_crs) > 0L) problems <- c(problems, paste0("Sin CRS definido (asígnelo en origen, no lo adivine): ", no_crs$file))

predictors <- active |> dplyr::filter(lower_clean(role) == "predictor")
bad_nature <- predictors |> dplyr::filter(!lower_clean(data_nature) %in% c("continuous", "categorical"))
if (nrow(bad_nature) > 0L) problems <- c(problems, paste0("data_nature debe ser continuous/categorical en: ", bad_nature$file))
if (sum(upper_clean(predictors$is_reference_grid) == "YES") > 1L) problems <- c(problems, "Solo un raster puede tener is_reference_grid = YES")

if (length(problems) > 0L) {
  human_gate(paste0(
    "Complete outputs/decisions/04_data_provenance.csv:\n", paste0("  - ", problems, collapse = "\n"),
    "\nRoles válidos: ", paste(ROLES, collapse = ", ")
  ), "04_datos")
}

# ------------------------------------------------------------
# 4. Diagnóstico: CRS y coincidencia temporal
# ------------------------------------------------------------

proj_code <- sub("^EPSG:", "", project_crs)
crs_report <- spatial_active |>
  dplyr::mutate(matches_project_crs = crs == proj_code, action = ifelse(matches_project_crs, "none", "reproject in stage 05"))
write_csv_safe(crs_report, file.path(DIR$tables, "04_crs_report.csv"))

# ADAPTACIÓN A1: as.Date() falla con error ante fechas parciales ("2000-02", "2020-2021") o textos; se normalizan año y año-mes y lo demás queda NA.
parse_date <- function(x) {
  x <- substr(clean_chr(x), 1, 10)
  x <- ifelse(grepl("^\\d{4}-\\d{2}$", x), paste0(x, "-01"), ifelse(grepl("^\\d{4}$", x), paste0(x, "-01-01"), x))
  suppressWarnings(as.Date(x, optional = TRUE))
}
dates <- parse_date(active$acquisition_date)
temporal <- active |>
  dplyr::mutate(date = dates) |>
  dplyr::group_by(role) |>
  dplyr::summarise(first = suppressWarnings(as.character(min(date, na.rm = TRUE))), last = suppressWarnings(as.character(max(date, na.rm = TRUE))), n = dplyr::n(), .groups = "drop")
write_csv_safe(temporal, file.path(DIR$tables, "04_temporal_coverage.csv"))

date_gap_days <- if (all(c("predictor", "reference_training") %in% lower_clean(active$role))) {
  p <- dates[lower_clean(active$role) == "predictor"]; r <- dates[lower_clean(active$role) == "reference_training"]
  if (all(is.na(p)) || all(is.na(r))) NA else as.numeric(abs(stats::median(p, na.rm = TRUE) - stats::median(r, na.rm = TRUE)))
} else NA
if (is.finite(date_gap_days) && date_gap_days > 90) {
  warning("Desfase mediano de ", date_gap_days, " días entre predictores y referencia. Justifíquelo en Métodos/Limitaciones.", call. = FALSE)
}

write_md(file.path(DIR$evidence, "04_data_pack.md"), "Data inventory and provenance", list(
  "Datasets and provenance" = prov |> dplyr::filter(lower_clean(role) != "exclude") |>
    dplyr::select(file, data_type, role, platform_sensor, acquisition_date, processing_level, source, license, citation_or_doi),
  "Technical metadata" = inventory |> dplyr::select(dplyr::any_of(c("file", "crs", "res_x", "res_y", "n_layers", "n_features", "n_points", "bbox_density_pts_m2"))),
  "CRS harmonisation" = crs_report |> dplyr::select(file, crs, matches_project_crs, action),
  "Temporal coverage by role" = temporal,
  "Median predictor-reference date gap (days)" = as.character(date_gap_days)
))

mark_stage("04_datos", "PASSED", paste0(nrow(active), " archivos activos"))
cat("\nEtapa 04 completada. Siguiente: 05_preprocesamiento_QC.R\n")
