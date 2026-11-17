# a1_rev_00b_uso_edificios.R: cuenta las categorias de uso y 'lote vacio' entre los edificios fuera de manzanas (insumo de la sensibilidad por tipo de edificio pedida por el revisor 2). No modifica resultados. Fecha: 2026-11-07
suppressMessages(library(sf))
ROOT <- "C:/Trabajos_y_Publicaciones/02_cambio_ciudades/000_A1_research"
b <- st_read(file.path(ROOT, "data/raw/vector/censo_cartografia_aoi.gpkg"), layer = "edificios_2021_puntos", quiet = TRUE)
b <- st_drop_geometry(b); d <- b[b$en_manzana == 0, ]
cat("dispersos (todo el gpkg):", nrow(d), "\n")
print(table(cod_otros = d$cod_otros, otros = d$otros, useNA = "ifany"))
cat("edif_censo:\n"); print(summary(d$edif_censo)); print(head(table(d$edif_censo), 10))
cat("nom_edif no vacio:", sum(nzchar(trimws(d$nom_edif))), "\n"); print(head(sort(table(d$nom_edif[nzchar(trimws(d$nom_edif))]), decreasing = TRUE), 15))
