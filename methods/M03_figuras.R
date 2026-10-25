# ============================================================
# methods/M03_figuras.R
# Figuras del módulo M03 (presencia-fondo). Las dibuja el propio módulo al terminar (06) y también la etapa 08 a partir de las tablas ya guardadas,
# de modo que cambiar un estilo no obliga a repetir el modelo (43 min). Reglas de estilo de la cadena (lecciones de la auditoría de A1):
#   - paleta segura para daltonismo (Okabe-Ito, R_helpers_publicacion.R) con símbolos distintos, nunca rojo y verde juntos;
#   - etiquetas en inglés legible (no códigos internos como «dist_nucleo» ni nombres con guion bajo);
#   - la validación por centro muestra el modelo final (spline) y la métrica del criterio (AUC tras muestreo por pares), con su umbral.
# Entradas : data frames (tablas tab_M03_*) o los CSV de outputs/tables. Salidas: objetos ggplot y, con m03_render_figures(), PNG/SVG.
# ============================================================

M03_VAR_LABEL <- c(slope = "Slope", dist_road_igm_vias = "Distance to IGM roads (no tracks)", dist_road_igm_all = "Distance to IGM roads, tracks and unpaved roads",
                   dist_road_osm2019 = "Distance to OSM roads (2019)", dist_road_census2020 = "Distance to census road axes", dist_nucleo = "Distance to urban core",
                   dist_river = "Distance to rivers", log_d_igm_vias = "Distance to IGM roads (no tracks)", log_d_nucleo = "Distance to urban core")
M03_CENTRE_LABEL <- c(Catamayo = "Catamayo", Loja = "Loja", Malacatos = "Malacatos", San_Pedro = "San Pedro de la Bendita", Vilcabamba = "Vilcabamba", fuera_AOI = "Not assigned to a centre")

m03_fig_effects <- function(coef_tab) {
  eff <- coef_tab
  eff$var_f <- factor(M03_VAR_LABEL[eff$variable], levels = rev(unique(M03_VAR_LABEL[eff$variable])))
  fs <- sort(unique(eff$feature_set))
  ggplot2::ggplot(eff, ggplot2::aes(intensity_ratio, var_f, colour = feature_set, shape = feature_set)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = ir_boot_lo, xmax = ir_boot_hi), position = ggplot2::position_dodge(width = 0.6), size = 0.3) +
    ggplot2::scale_x_log10() + cb_scales(fs, "Feature set") +
    ggplot2::labs(x = "Intensity ratio per SD (log scale; block-bootstrap 95% interval)", y = NULL) + theme_paper()
}

m03_fig_response <- function(resp_tab) {
  rs <- resp_tab
  rs$var_f <- factor(M03_VAR_LABEL[rs$variable], levels = unique(unname(M03_VAR_LABEL[c("slope", "dist_road_igm_vias", "dist_nucleo")])))
  ggplot2::ggplot(rs, ggplot2::aes(x, relative_intensity)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) + ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey50") + ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~var_f, scales = "free_x") +
    ggplot2::labs(x = "Predictor value (slope in degrees; distances in metres, 1st-99th percentile of the domain)", y = "Relative intensity (log scale; reference = median)") + theme_paper()
}

m03_fig_validation <- function(cv_fold_tab) {
  cvl <- cv_fold_tab |> dplyr::select(feature_set, scheme, auc, cauc_pds, cauc, boyce) |>
    tidyr::pivot_longer(c(auc, cauc_pds, cauc, boyce), names_to = "metric", values_to = "value") |>
    dplyr::mutate(metric = factor(metric, levels = c("auc", "cauc_pds", "cauc", "boyce"),
                                  labels = c("AUC", "AUC after pairwise distance\nsampling, calibrated", "AUC calibrated with\ndistance null model", "Boyce index")),
                  scheme = factor(scheme, levels = c("random", "spatial_block"), labels = c("Random", "Spatial block")))
  ggplot2::ggplot(cvl, ggplot2::aes(feature_set, value, colour = scheme, shape = scheme)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), alpha = 0.75, size = 1.2) +
    ggplot2::facet_wrap(~metric, scales = "free_y") + ggplot2::coord_flip() +
    ggplot2::scale_colour_manual(name = "CV scheme", values = c("Random" = OKABE_ITO[["orange"]], "Spatial block" = OKABE_ITO[["blue"]])) +
    ggplot2::scale_shape_manual(name = "CV scheme", values = c("Random" = 17, "Spatial block" = 16)) +
    ggplot2::labs(x = NULL, y = NULL) + theme_paper()
}

# Validación dejando un centro fuera: modelo final (spline) y métrica del criterio (AUC tras muestreo por pares); la línea discontinua es el umbral
m03_fig_centres <- function(loco_spline, threshold = NA_real_) {
  lc <- tibble::tibble(centre = unname(M03_CENTRE_LABEL[as.character(loco_spline$centre)]), `AUC` = as.numeric(loco_spline$auc),
                       `AUC after pairwise distance sampling (cAUC_pds)` = as.numeric(loco_spline$cauc_pds), `Boyce index` = as.numeric(loco_spline$boyce)) |>
    tidyr::pivot_longer(-centre, names_to = "metric", values_to = "value")
  lc$metric <- factor(lc$metric, levels = c("AUC", "AUC after pairwise distance sampling (cAUC_pds)", "Boyce index"))
  ggplot2::ggplot(lc, ggplot2::aes(centre, value, colour = metric, shape = metric)) +
    list(if (is.finite(threshold)) ggplot2::geom_hline(yintercept = threshold, linetype = 2, colour = "grey40")) +
    ggplot2::geom_point(size = 2.4, position = ggplot2::position_dodge(width = 0.5)) + ggplot2::coord_flip() +
    cb_scales(levels(lc$metric)) +
    ggplot2::labs(x = NULL, y = if (is.finite(threshold)) paste0("Value (dashed line: criterion ", threshold, " for cAUC_pds)") else "Value", colour = NULL, shape = NULL) + theme_paper() +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
}

m03_fig_convergence <- function(conv_tab) {
  cv0 <- conv_tab
  cv0$var_f <- factor(M03_VAR_LABEL[cv0$variable], levels = unique(unname(M03_VAR_LABEL[c("slope", "log_d_igm_vias", "log_d_nucleo")])))
  ggplot2::ggplot(cv0, ggplot2::aes(ratio, intensity_ratio, colour = var_f, shape = var_f)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2) + ggplot2::scale_x_log10() + cb_scales(levels(cv0$var_f)) +
    ggplot2::labs(x = "Background points per presence (log scale)", y = "Intensity ratio per SD", colour = NULL, shape = NULL) + theme_paper()
}

# Dominio de modelación y edificios dispersos por centro (a partir de los vectores procesados; el dominio excluye manzanas censales y áreas protegidas)
m03_fig_domain <- function(DIR, root, card) {
  dom_f <- file.path(DIR$processed_vectors, "model_domain.gpkg"); tr_f <- file.path(DIR$processed_vectors, "training.gpkg")
  if (!file.exists(dom_f) || !file.exists(tr_f)) return(NULL)
  dom_main <- sf::st_read(dom_f, layer = "main", quiet = TRUE)
  cores <- sf::st_read(dom_f, layer = "urban_cores", quiet = TRUE)
  dom_geom <- sf::st_union(dom_main)
  par <- card_get(card, "method_params.M03_presencia_fondo.data", list())
  sub <- tryCatch({
    cen_f <- file.path(root, par$census_gpkg %||% "data/raw/vector/censo_cartografia_aoi.gpkg")
    mz <- sf::st_make_valid(sf::st_read(cen_f, layer = par$blocks_layer %||% "manzanas_2021", quiet = TRUE))
    mz <- sf::st_transform(mz, sf::st_crs(dom_geom)); mz <- mz[lengths(sf::st_intersects(mz, dom_geom)) > 0L, ]
    g <- sf::st_difference(dom_geom, sf::st_union(sf::st_geometry(mz)))
    ex_f <- file.path(root, par$exclude_vector %||% "")
    if (nzchar(par$exclude_vector %||% "") && file.exists(ex_f)) g <- sf::st_difference(g, sf::st_union(sf::st_transform(sf::st_make_valid(sf::st_read(ex_f, quiet = TRUE)), sf::st_crs(dom_geom))))
    g
  }, error = function(e) { message("Figura del dominio: no se pudieron restar manzanas ni áreas protegidas: ", conditionMessage(e)); dom_geom })
  tr <- sf::st_read(tr_f, quiet = TRUE)
  tr <- tr[!is.na(tr[[par$in_block_field %||% "en_manzana"]]) & tr[[par$in_block_field %||% "en_manzana"]] == 0, ]
  tr <- tr[lengths(sf::st_intersects(tr, sf::st_transform(dom_geom, sf::st_crs(tr)))) > 0L, ]
  cen <- as.character(tr$centro_aoi); cen[is.na(cen) | !nzchar(trimws(cen))] <- "fuera_AOI"
  tr$centre <- factor(unname(M03_CENTRE_LABEL[cen]), levels = unname(M03_CENTRE_LABEL))
  L_D <- "Main modelling domain"; L_C <- "Urban cores"
  lev <- levels(tr$centre)
  pal <- stats::setNames(c(cb_values(5), "#4D4D4D")[seq_along(lev)], lev)
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = sf::st_sf(geometry = sub), ggplot2::aes(fill = L_D), colour = NA) +
    ggplot2::geom_sf(data = cores, ggplot2::aes(fill = L_C), colour = NA, alpha = 0.7) +
    ggplot2::geom_sf(data = tr, ggplot2::aes(colour = centre), shape = 16, size = 0.35, alpha = 0.75) +
    ggplot2::scale_fill_manual(name = "Layer", breaks = c(L_D, L_C), values = stats::setNames(c("grey88", "firebrick"), c(L_D, L_C))) +
    ggplot2::scale_colour_manual(name = "Centre", drop = FALSE, values = pal) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1, override.aes = list(size = 2.2, alpha = 1)), fill = ggplot2::guide_legend(order = 2)) +
    ggplot2::labs(x = NULL, y = NULL) + theme_paper() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.2), legend.title = ggplot2::element_text(face = "bold", size = 10),
                   legend.text = ggplot2::element_text(size = 9), legend.key.size = grid::unit(0.5, "cm"))
  if (requireNamespace("ggspatial", quietly = TRUE)) {
    p <- p + ggspatial::annotation_scale(location = "bl", width_hint = 0.25, text_family = FIG_FONT) +
      ggspatial::annotation_north_arrow(location = "tr", height = grid::unit(0.8, "cm"), width = grid::unit(0.6, "cm"), style = ggspatial::north_arrow_minimal())
  }
  p
}

# Dibuja y guarda las figuras del módulo que salen de tablas (las demás, p. ej. correlograma o residuos, las dibuja el módulo)
m03_render_figures <- function(DIR, card, root, dpi = 600) {
  rdt <- function(n) readr::read_csv(file.path(DIR$tables, n), show_col_types = FALSE, progress = FALSE)
  out <- character()
  guardar <- function(p, id, w, h) { if (!is.null(p)) { save_fig(p, id, w, h, dpi = dpi); out <<- c(out, id) } }
  guardar(m03_fig_effects(rdt("tab_M03_coefficients.csv")), "fig_M03_effects", 13, 8)
  guardar(m03_fig_response(rdt("tab_M03_response_curves.csv")), "fig_M03_response_curves", 17, 7)
  guardar(m03_fig_validation(rdt("tab_M03_cv_folds.csv")), "fig_M03_validation", 17, 7)
  # umbral del criterio de validación por centro (H6), tomado del protocolo congelado
  thr <- tryCatch({ h <- Filter(function(x) identical(as.character(x$subset), "spl_center_min"), read_frozen_protocol()$payload$hypotheses); as.numeric(h[[1]]$threshold) }, error = function(e) NA_real_)
  guardar(m03_fig_centres(rdt("tab_M03_loco_spline.csv"), threshold = thr), "fig_M03_centres", 13, 8)
  guardar(m03_fig_convergence(rdt("tab_M03_background_convergence.csv")), "fig_M03_background_convergence", 12, 8)
  guardar(tryCatch(m03_fig_domain(DIR, root, card), error = function(e) { message("Figura del dominio no regenerada: ", conditionMessage(e)); NULL }), "fig_M03_domain", 15, 11)
  out
}

# Convención de la cadena: la etapa 08 busca render_method_figures() en methods/<prefijo>_figuras.R
render_method_figures <- m03_render_figures
