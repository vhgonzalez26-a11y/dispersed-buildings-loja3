# a1_rev_02_redes_viales.R
# Objetivo : (revisor 2, comentario 4) caracterizar las cuatro redes viales dentro del dominio principal: longitud, densidad, composicion por clase y cobertura espacial
#            (celdas de 1 km con al menos una via), con una correspondencia de clases equivalentes entre IGM y OSM; y calcular la distancia de cada punto a la red OSM 'vehicular'
#            (todo el OSM salvo sendas, caminos, escaleras y ciclovias) para poder comparar clases equivalentes en el modelo.
# Entradas : data/raw/vector/* (IGM 1:50 000, OSM 2019-01, ejes censales 2021, manzanas 2021, Podocarpus, AOI), cache M03_dataset_ea8ce5b426.rds
# Salidas  : outputs/revision/rev_redes_composicion.csv, rev_redes_resumen.csv, data/processed/REV_dist_osm_vehicular.rds
# Fecha    : 2026-11-07
suppressMessages({library(sf); library(terra); library(dplyr); library(purrr); library(tibble)})
ROOT <- "C:/Trabajos_y_Publicaciones/02_cambio_ciudades/000_A1_research"
OUT <- file.path(ROOT, "outputs/revision"); say <- function(...) cat(sprintf("[rev02 %s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
CRS <- 32717; v <- function(...) file.path(ROOT, "data/raw/vector", ...)
rd <- function(f, l = NULL) { x <- if (is.null(l)) st_read(v(f), quiet = TRUE) else st_read(v(f), layer = l, quiet = TRUE)
  if (any(as.character(st_geometry_type(x)) %in% c("POLYGON", "MULTIPOLYGON"))) x <- st_make_valid(x); st_transform(x, CRS) }
# dominio principal (misma construccion que M03: AOI amplio + 5 km, menos manzanas 2021 y Podocarpus, limitado al ráster de pendiente)
say("Construyendo el dominio principal...")
excl <- st_union(rd("areas_protegidas_snap_podocarpus_oficial.gpkg")); mz <- st_union(rd("censo_cartografia_aoi.gpkg", "manzanas_2021")); bad <- st_union(excl, mz)
aoi_amp <- st_union(rd("aoi_estudio.gpkg", "aoi_amplio")); ext_pre <- st_buffer(aoi_amp, 5000)
dom <- st_union(st_difference(ext_pre, bad), st_difference(aoi_amp, bad))
re <- as.vector(ext(rast(file.path(ROOT, "data/raw/raster/slope_srtm30_x2.tif"))))
dom <- st_intersection(dom, st_as_sfc(st_bbox(c(xmin = re[[1]], xmax = re[[2]], ymin = re[[3]], ymax = re[[4]]), crs = st_crs(CRS))))
area_km2 <- as.numeric(st_area(dom)) / 1e6; say("Area del dominio: ", round(area_km2, 1), " km2 (esperado 1639.0)")
stopifnot(abs(area_km2 - 1639.05) < 1)

# redes y clases
igm <- rd("igm_50k_2013_base.gpkg", "via_l"); igm$clase <- c("RUTA PRIMARIA" = "major", "RUTA SECUNDARIA" = "major", "RUTA LOCAL" = "local")[igm$hct_desc]
igm$detalle <- igm$hct_desc
sen <- rd("igm_50k_2013_base.gpkg", "sendero_l"); sen$clase <- "track_path"; sen$detalle <- "SENDERO"
rod <- rd("igm_50k_2013_base.gpkg", "rodera_l"); rod$clase <- "track_path"; rod$detalle <- "RODERA"
osm <- rd("osm_vias_snapshots.gpkg", "osm_2019_01")
maj <- c("trunk", "trunk_link", "primary", "primary_link", "secondary", "secondary_link")
loc <- c("tertiary", "tertiary_link", "residential", "unclassified", "living_street", "service", "road")
osm$clase <- ifelse(osm$highway %in% maj, "major", ifelse(osm$highway %in% loc, "local", "track_path")); osm$detalle <- osm$highway
cen <- rd("censo_cartografia_aoi.gpkg", "ejes_2021"); cen$clase <- "no_hierarchy"; cen$detalle <- as.character(cen$tipo_eje)
nets <- list(
  igm_vias = igm[, c("clase", "detalle")], igm_all = dplyr::bind_rows(igm[, c("clase", "detalle")], sen[, c("clase", "detalle")], rod[, c("clase", "detalle")]),
  osm2019 = osm[, c("clase", "detalle")], census2020 = cen[, c("clase", "detalle")])
# cuadricula de 1 km sobre el dominio
grid <- st_make_grid(dom, cellsize = 1000); grid <- grid[lengths(st_intersects(grid, dom)) > 0]
comp <- purrr::imap_dfr(nets, function(n, nm) {
  say("Red ", nm, " (", nrow(n), " tramos)...")
  n2 <- suppressWarnings(st_intersection(n, st_sf(geometry = dom))); n2 <- n2[!st_is_empty(n2), ]
  n2$len_km <- as.numeric(st_length(n2)) / 1000
  n2 |> st_drop_geometry() |> dplyr::group_by(clase, detalle) |> dplyr::summarise(length_km = sum(len_km), n_segments = dplyr::n(), .groups = "drop") |> dplyr::mutate(network = nm)
})
comp <- comp |> dplyr::group_by(network) |> dplyr::mutate(share_of_length = length_km / sum(length_km)) |> dplyr::ungroup() |> dplyr::mutate(density_km_per_km2 = length_km / area_km2)
write.csv(comp, file.path(OUT, "rev_redes_composicion.csv"), row.names = FALSE)
res <- purrr::imap_dfr(nets, function(n, nm) {
  nl <- suppressWarnings(st_intersection(n, st_sf(geometry = dom)))
  cov <- mean(lengths(st_intersects(grid, nl)) > 0)
  cm <- comp |> dplyr::filter(network == nm)
  tibble::tibble(network = nm, length_km = sum(cm$length_km), density_km_per_km2 = sum(cm$length_km) / area_km2, coverage_1km_cells = cov,
                 major_km = sum(cm$length_km[cm$clase == "major"]), local_km = sum(cm$length_km[cm$clase == "local"]), track_path_km = sum(cm$length_km[cm$clase == "track_path"]),
                 no_hierarchy_km = sum(cm$length_km[cm$clase == "no_hierarchy"]), n_cells = length(grid))
})
write.csv(res, file.path(OUT, "rev_redes_resumen.csv"), row.names = FALSE); print(as.data.frame(res)); print(as.data.frame(comp |> dplyr::select(network, clase, detalle, length_km, share_of_length)))

# distancia a la red OSM vehicular (OSM sin sendas, caminos, escaleras ni ciclovias) para presencias y fondo del dominio principal
say("Distancia a la red OSM vehicular...")
osm_v <- st_geometry(osm[osm$clase != "track_path", ]); osm_v <- osm_v[!st_is_empty(osm_v)]
cache <- readRDS(file.path(ROOT, "data/processed/M03_dataset_ea8ce5b426.rds"))
P0 <- cache$pres; D0 <- cache$pool; P0 <- P0[is.finite(P0$slope), ]; D0 <- D0[is.finite(D0$slope), ]
P <- P0[P0$in_extended, ]; D <- D0[D0$in_extended, ]; D <- D[seq_len(ceiling(20 * nrow(P))), ]
dist_to <- function(x, y) { g <- st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = CRS) |> st_geometry(); idx <- st_nearest_feature(g, osm_v)
  as.numeric(terra::distance(terra::vect(g), terra::vect(osm_v[idx]), pairwise = TRUE)) }
dP <- dist_to(P$x, P$y); dD <- dist_to(D$x, D$y)
# control: la red vehicular nunca queda mas lejos que la red completa
say("Control: media de d_osm2019 (completa) ", round(mean(P$d_osm2019)), " m; d vehicular ", round(mean(dP)), " m; vehicular >= completa: ", all(dP >= P$d_osm2019 - 1e-6))
saveRDS(list(P = tibble::tibble(x = P$x, y = P$y, d_osm_veh = dP, log_d_osm_veh = log1p(dP)), D = tibble::tibble(x = D$x, y = D$y, d_osm_veh = dD, log_d_osm_veh = log1p(dD))),
        file.path(ROOT, "data/processed/REV_dist_osm_vehicular.rds"))
say("Fin.")
