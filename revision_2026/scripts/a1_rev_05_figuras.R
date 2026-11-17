# a1_rev_05_figuras.R
# Objetivo : dos figuras nuevas para la revision de Land: (1) ubicacion de los pliegues de la validacion espacial por bloques (repeticion 1) y de los centros de la validacion dejando un centro fuera
#            (revisor 1: mostrar donde estan los pliegues); (2) sintesis de los efectos de pendiente, vias y nucleo urbano por centro frente al valor del dominio principal (revisor 2, comentario 6).
# Entradas : outputs/revision/rev_fold_blocks.csv, data/processed/M03_dataset_ea8ce5b426.rds, outputs/tables/tab_M03_centre_effects.csv, outputs/tables/tab_M03_coefficients.csv
# Salidas  : outputs/figures/fig_REV_folds.(png|svg), fig_REV_centre_synthesis.(png|svg)
# Estilo   : paleta Okabe-Ito con simbolos distintos, Arial 10 pt (project_card figures), etiquetas en ingles
# Fecha    : 2026-11-07
ROOT <- "C:/Trabajos_y_Publicaciones/02_cambio_ciudades/000_A1_research"; setwd(ROOT)
suppressMessages(source(file.path(ROOT, "00_setup.R"), encoding = "UTF-8")); suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
source(file.path(ROOT, "methods/M03_figuras.R"), encoding = "UTF-8")
blocks <- readr::read_csv(file.path(DIR$decisions, "..", "revision", "rev_fold_blocks.csv"), show_col_types = FALSE) |> dplyr::filter(repeat_id == 1)
cache <- readRDS(file.path(ROOT, "data/processed/M03_dataset_ea8ce5b426.rds")); P <- cache$pres; P <- P[is.finite(P$slope) & P$in_extended, ]
bs <- 10000
blocks <- blocks |> tidyr::separate(blk, c("bx", "by"), sep = " ", convert = TRUE) |> dplyr::mutate(xmin = bx * bs / 1000, xmax = (bx + 1) * bs / 1000, ymin = by * bs / 1000, ymax = (by + 1) * bs / 1000, Fold = factor(fold))
P$xk <- P$x / 1000; P$yk <- P$y / 1000
lab_panel <- c(a = "(a) Spatial block CV:
folds of 10 km blocks (repeat 1)", b = "(b) Leave-one-centre-out:
centre of each building")
pa <- P |> dplyr::mutate(panel = lab_panel[["a"]]); pb <- P |> dplyr::mutate(panel = lab_panel[["b"]], Centre = factor(unname(M03_CENTRE_LABEL[centre_census]), levels = unname(M03_CENTRE_LABEL[c("Catamayo", "Loja", "Malacatos", "San_Pedro", "Vilcabamba", "fuera_AOI")])))
levels(pb$Centre)[6] <- "Not assigned"
bl_a <- blocks |> dplyr::mutate(panel = lab_panel[["a"]]); bl_b <- blocks |> dplyr::mutate(panel = lab_panel[["b"]])
cols_fold <- stats::setNames(unname(OKABE_ITO[c("blue", "orange", "green", "pink", "yellow")]), 1:5)
cols_cen <- stats::setNames(c(unname(OKABE_ITO[c("blue", "orange", "green", "pink", "yellow")]), "grey55"), levels(pb$Centre))
shp_cen <- stats::setNames(c(16, 17, 15, 18, 8, 4), levels(pb$Centre))
fig_folds <- ggplot() +
  geom_rect(data = bl_a, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Fold), alpha = 0.55, colour = "grey30", linewidth = 0.2) +
  geom_rect(data = bl_b, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, colour = "grey30", linewidth = 0.2) +
  geom_point(data = pa, aes(xk, yk), size = 0.15, colour = "grey15", alpha = 0.6) +
  geom_point(data = pb, aes(xk, yk, colour = Centre, shape = Centre), size = 0.35, alpha = 0.8) +
  facet_wrap(~panel, ncol = 2) + coord_equal() +
  scale_fill_manual(name = "Fold (blocks)", values = cols_fold) + scale_colour_manual(name = "Centre (buildings)", values = cols_cen) + scale_shape_manual(name = "Centre (buildings)", values = shp_cen) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2, order = 2), shape = guide_legend(nrow = 2, order = 2), fill = guide_legend(nrow = 1, order = 1)) +
  labs(x = "Easting (km, UTM zone 17S)", y = "Northing (km)") + theme_paper() + theme(legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 9), legend.text = element_text(size = 9), axis.text = element_text(size = 8), strip.text = element_text(size = 9, hjust = 0), panel.spacing = unit(0.6, "cm"))
save_fig(fig_folds, "fig_REV_folds", 16, 12.5)

ce <- readr::read_csv(file.path(DIR$tables, "tab_M03_centre_effects.csv"), show_col_types = FALSE)
co <- readr::read_csv(file.path(DIR$tables, "tab_M03_coefficients.csv"), show_col_types = FALSE) |> dplyr::filter(feature_set == "igm_vias") |> dplyr::transmute(variable = c(slope = "slope", dist_road_igm_vias = "log_d_igm_vias", dist_nucleo = "log_d_nucleo")[variable], pooled = intensity_ratio)
ce <- ce |> dplyr::left_join(co, by = "variable") |> dplyr::mutate(Variable = factor(c(slope = "Slope", log_d_igm_vias = "Distance to IGM roads", log_d_nucleo = "Distance to urban core")[variable], levels = c("Slope", "Distance to IGM roads", "Distance to urban core")),
  Centre = factor(unname(M03_CENTRE_LABEL[centre]), levels = unname(M03_CENTRE_LABEL[c("Catamayo", "Loja", "Malacatos", "San_Pedro", "Vilcabamba")])))
if (any(is.na(ce$Variable))) { print(unique(ce$variable)); stop("variables inesperadas en tab_M03_centre_effects") }
fig_syn <- ggplot(ce, aes(intensity_ratio, Centre, colour = Centre, shape = Centre)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  geom_vline(aes(xintercept = pooled), colour = "grey20", linewidth = 0.5) +
  geom_point(size = 2.4) + scale_x_log10(breaks = c(0.4, 0.7, 1, 1.25), labels = c("0.4", "0.7", "1", "1.25")) +
  facet_wrap(~Variable, ncol = 3) + scale_colour_manual(values = cols_cen[levels(ce$Centre)], guide = "none") + scale_shape_manual(values = shp_cen[levels(ce$Centre)], guide = "none") +
  scale_y_discrete(limits = rev) +
  labs(x = "Intensity ratio per SD (log scale)", y = NULL) + theme_paper() + theme(axis.text = element_text(size = 9), strip.text = element_text(size = 9))
save_fig(fig_syn, "fig_REV_centre_synthesis", 17, 6.5)
cat("figuras guardadas en", DIR$figures, "\n")
