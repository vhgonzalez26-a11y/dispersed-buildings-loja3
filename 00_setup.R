# ============================================================
# 00_setup.R
# Estructura del proyecto, paquetes base, ficha y bitácora
# Todos los demás scripts cargan este archivo.
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

detect_root <- function() {
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f) && nzchar(f)) return(dirname(normalizePath(f, winslash = "/")))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

if (!exists("ROOT") || !file.exists(file.path(ROOT, "00_setup.R"))) ROOT <- detect_root()
CARD_FILE <- file.path(ROOT, "project_card.yaml")

# ------------------------------------------------------------
# 1. Carpetas
# ------------------------------------------------------------

DIR <- list(
  raw               = file.path(ROOT, "data", "raw"),
  raw_raster        = file.path(ROOT, "data", "raw", "raster"),
  raw_vector        = file.path(ROOT, "data", "raw", "vector"),
  raw_pointcloud    = file.path(ROOT, "data", "raw", "pointcloud"),
  raw_tabular       = file.path(ROOT, "data", "raw", "tabular"),
  processed         = file.path(ROOT, "data", "processed"),
  processed_layers  = file.path(ROOT, "data", "processed", "layers"),
  processed_vectors = file.path(ROOT, "data", "processed", "vectors"),
  ideas             = file.path(ROOT, "ideas"),
  literature        = file.path(ROOT, "literature", "ris"),
  methods           = file.path(ROOT, "methods"),
  protocol          = file.path(ROOT, "outputs", "protocol"),
  decisions         = file.path(ROOT, "outputs", "decisions"),
  tables            = file.path(ROOT, "outputs", "tables"),
  figures           = file.path(ROOT, "outputs", "figures"),
  method            = file.path(ROOT, "outputs", "method"),
  evidence          = file.path(ROOT, "outputs", "evidence"),
  logs              = file.path(ROOT, "outputs", "logs")
)

invisible(lapply(DIR, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)))

# ------------------------------------------------------------
# 2. Paquetes base (los espaciales se instalan cuando se necesitan)
# ------------------------------------------------------------

CORE_PKGS <- c("dplyr", "tidyr", "purrr", "readr", "stringr", "stringi", "tibble",
               "openxlsx", "ggplot2", "svglite", "yaml", "jsonlite")

missing_core <- CORE_PKGS[!vapply(CORE_PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_core) > 0L) {
  message("Instalando paquetes base: ", paste(missing_core, collapse = ", "))
  install.packages(missing_core, repos = "https://cloud.r-project.org")
}

source(file.path(ROOT, "R_helpers.R"), encoding = "UTF-8")
# Ayudas genéricas para dejar listo un artículo: paleta segura para daltonismo, rastreo de cifras, suplemento en Word (ver LEEME_PLANTILLA.md)
if (file.exists(file.path(ROOT, "R_helpers_publicacion.R"))) source(file.path(ROOT, "R_helpers_publicacion.R"), encoding = "UTF-8")

# Paquetes por etapa (referencia)
STAGE_PKGS <- list(
  spatial   = c("terra", "sf"),
  M01       = c("terra", "sf", "ranger"),
  M02       = c("terra", "sf", "lidR"),
  maps      = c("terra", "sf", "ggspatial"),
  api       = c("httr2"),
  docx      = c("officer", "flextable")
)

# ------------------------------------------------------------
# 3. Ficha del proyecto y ajustes globales
# ------------------------------------------------------------

CARD <- read_card(CARD_FILE)
set.seed(as.integer(card_get(CARD, "project.random_seed", 2026)))


FIG_FONT <- card_get(CARD, "figures.font_family", "serif")
if (requireNamespace("systemfonts", quietly = TRUE)) {
  if (!any(grepl(FIG_FONT, systemfonts::system_fonts()$family, fixed = TRUE))) {
    warning("Fuente '", FIG_FONT, "' no detectada; se usará 'serif'.", call. = FALSE)
    FIG_FONT <- "serif"
  }
}

readme_raw <- file.path(DIR$raw, "LEEME_datos.txt")
if (!file.exists(readme_raw)) {
  writeLines(c(
    "DATOS ORIGINALES — NO EDITAR NI SOBRESCRIBIR",
    "raster/      GeoTIFF de satélite, ortomosaicos, DEM (.tif)",
    "vector/      AOI, puntos/polígonos de entrenamiento y validación, puntos GNSS (.gpkg/.shp)",
    "pointcloud/  nubes de puntos LiDAR o fotogramétricas (.las/.laz)",
    "tabular/     mediciones de campo (.csv)",
    "Toda transformación se hace con los scripts y queda en data/processed."
  ), readme_raw)
}

writeLines(capture.output(sessionInfo()), file.path(DIR$logs, "session_info.txt"))

cat("\n============================================================\n")
cat("PIPELINE DE INVESTIGACIÓN GEOESPACIAL — ", card_get(CARD, "project.short_name", ""), "\n", sep = "")
cat("============================================================\n")
cat("Raíz: ", ROOT, "\n\n", sep = "")
pipeline_status()
