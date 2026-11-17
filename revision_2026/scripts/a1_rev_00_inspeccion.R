# a1_rev_00_inspeccion.R: inspecciona campos y clases de las capas viales y el campo de uso de los edificios, y comprueba la cache M03 (insumos para los analisis de la revision de Land). No modifica resultados. Fecha: 2026-11-07
suppressMessages(library(sf))
ROOT <- "C:/Trabajos_y_Publicaciones/02_cambio_ciudades/000_A1_research"
v <- function(...) file.path(ROOT, "data/raw/vector", ...)
cat("== IGM via_l\n"); x <- st_read(v("igm_50k_2013_base.gpkg"), layer = "via_l", quiet = TRUE); print(names(x)); print(table(x$hct_desc, useNA = "ifany"))
for (l in c("sendero_l", "rodera_l")) { y <- st_read(v("igm_50k_2013_base.gpkg"), layer = l, quiet = TRUE); cat(l, nrow(y), "\n") }
cat("\n== OSM\n"); o <- st_read(v("osm_vias_snapshots.gpkg"), layer = "osm_2019_01", quiet = TRUE); print(names(o)); print(table(o$highway, useNA = "ifany"))
cat("\n== Censo ejes\n"); e <- st_read(v("censo_cartografia_aoi.gpkg"), layer = "ejes_2021", quiet = TRUE); print(names(e)); cat(nrow(e), "ejes\n")
cat("\n== Edificios\n"); b <- st_read(v("censo_cartografia_aoi.gpkg"), layer = "edificios_2021_puntos", quiet = TRUE); print(names(b)); cat(nrow(b), "\n")
for (nm in names(b)) { if (nm == "geom" || nm == "geometry") next; u <- unique(b[[nm]]); if (length(u) <= 25) { cat("--", nm, "\n"); print(table(b[[nm]], useNA = "ifany")) } }
cat("\n== Cache\n")
for (f in list.files(file.path(ROOT, "data/processed"), pattern = "^M03_dataset_.*[.]rds$", full.names = TRUE)) {
  d <- readRDS(f); cat(basename(f), format(file.mtime(f)), "| pres:", nrow(d$pres), "| pool:", nrow(d$pool), "\n")
  if (nrow(d$pres) == 17088 || TRUE) { cat("columnas pres:", paste(names(d$pres), collapse = ", "), "\n") }
}
