# ============================================================
# methods/M03_presencia_fondo.R
# Localización de edificios dispersos como proceso puntual de Poisson no
# homogéneo (IPP) con covariables, estimado con puntos de fondo (cuadratura),
# validado con bloques espaciales, modelo nulo de distancia y índice de Boyce.
#
# Documentación completa: methods/M03_presencia_fondo_README.md
# Cada decisión de diseño remite a una referencia verificada en texto completo
# (literature/Estado_metodologia_tras_textos_completos.md):
#   - Verosimilitud IPP y fondo grande ...... Warton & Shepherd 2010; Renner et al. 2015;
#                                              Fithian & Hastie 2013
#   - Comprobación de convergencia del fondo . Renner et al. 2015; Barbet-Massin et al. 2012;
#                                              Whitford et al. 2024
#   - Validación por bloques espaciales ...... Roberts et al. 2017; Valavi et al. 2019;
#                                              Pohjankukka et al. 2017
#   - Sesgo de ordenamiento y modelo nulo .... Hijmans 2012 (cAUC = AUC + 0.5 - max(0.5, AUC nulo))
#   - Límites del AUC; índice de Boyce ....... Lobo et al. 2008; Hirzel et al. 2006
#   - Colinealidad ........................... Dormann et al. 2013
#   - Autocorrelación espacial en cambio de uso Mets et al. 2017; Segurado et al. 2006
#   - Endogeneidad de las vías ............... Chomitz & Gray 1996; Duranton & Turner 2012
#   - Estructura de reporte .................. Zurell et al. 2020 (ODMAP)
#
# Contrato (methods/M00_PLANTILLA.R): devuelve metrics, tables, figures, facts,
# limitations, products. No escribe archivos finales salvo la caché de datos
# derivados en data/processed.
# ============================================================

MODULE_INFO <- list(
  id = "M03_presencia_fondo",
  description = paste(
    "Presence-background model of dispersed building locations estimated as an inhomogeneous Poisson point process",
    "(large uniform background, quadrature likelihood), with standardised intensity ratios and block-bootstrap intervals,",
    "spatial block cross-validation with distance-null calibration (cAUC) and continuous Boyce index,",
    "leave-one-centre-out validation, sensitivity to background size, study domain, block size and random seed,",
    "and descriptive shares of buildings near roads and above slope values."
  ),
  requires_stages = "04_datos",
  packages = c("sf", "terra"),
  # nombres exactos (la etapa 03 comprueba que cada hipótesis use una métrica de esta lista)
  metrics = c(paste0("intensity_ratio_", c("slope", "dist_road_igm_vias", "dist_road_igm_all", "dist_road_osm2019", "dist_road_census2020", "dist_nucleo", "dist_river")),
              paste0("intensity_ratio_iqr_", c("slope", "dist_road_igm_vias", "dist_road_igm_all", "dist_road_osm2019", "dist_road_census2020", "dist_nucleo", "dist_river")),
              paste0("odds_ratio_", c("slope", "dist_road_igm_vias", "dist_nucleo")),
              "auc_spatial_cv", "auc_random_cv", "auc_pds_spatial_cv", "auc_null_distance", "cauc_spatial_cv", "cauc_pds_spatial_cv", "boyce_spatial_cv",
              "cv_log_likelihood", "spatial_sorting_bias", "spatial_sorting_bias_after_pds", "log_likelihood", "aic",
              "share_within_250m_road", "share_within_25m_state_road", "share_within_30m_state_road", "share_above_slope_threshold", "share_inside_urban_limit", "building_density_per_km2",
              "n_presences", "n_background", "block_size_m", "residual_range_m"),
  # patrones: <feature_set> = igm_vias, igm_all, osm2019, census2020, full_igm_amenity; spl_<feature_set> = modelo con splines
  subsets = c("<feature_set>", "spl_<feature_set>", "center_<name>", "center_min", "spl_center_<name>", "spl_center_min", "ring_<a>_<b>km",
              "sens_domain_<name>", "sens_bgratio_<k>", "sens_seed_<s>", "sens_blocksize_<m>", "sens_spline_df_<k>", "sens_logit_3to1")
)

M03_VAR_LABEL <- c(slope = "slope", log_d_igm_vias = "dist_road_igm_vias", log_d_igm_all = "dist_road_igm_all",
                   log_d_osm2019 = "dist_road_osm2019", log_d_census2020 = "dist_road_census2020",
                   log_d_nucleo = "dist_nucleo", log_d_river = "dist_river")

validate_params <- function(p) {
  msg <- character()
  need <- c("data", "domains", "background", "cv", "bootstrap", "feature_sets", "final_feature_set")
  for (k in need) if (is.null(p[[k]])) msg <- c(msg, paste0("falta ", k))
  if (length(msg) > 0L) return(msg)
  for (k in c("census_gpkg", "presence_layer", "blocks_layer", "aoi_gpkg", "exclude_vector", "slope_raster", "igm_gpkg", "osm_gpkg", "osm_layer",
              "census_axes_layer", "core_layer")) {
    if (!nzchar(clean_chr(p$data[[k]] %||% ""))) msg <- c(msg, paste0("falta data$", k))
  }
  if (length(p$thresholds$road_distance_sensitivity_m) > 0L && length(p$state_networks) == 0L) msg <- c(msg, "falta state_networks para los umbrales de sensibilidad de distancia a vías")
  vars <- unique(unlist(p$feature_sets))
  bad <- setdiff(vars, names(M03_VAR_LABEL))
  if (length(bad) > 0L) msg <- c(msg, paste0("variables desconocidas en feature_sets: ", paste(bad, collapse = ", ")))
  if (!clean_chr(p$final_feature_set) %in% names(p$feature_sets)) msg <- c(msg, "final_feature_set debe ser un nombre de feature_sets")
  if (!is.numeric(p$background$ratio_main) || p$background$ratio_main < 1) msg <- c(msg, "background$ratio_main debe ser >= 1")
  if (!all(p$background$ratio_main %in% p$background$ladder)) msg <- c(msg, "background$ladder debe incluir ratio_main")
  if (!is.numeric(p$cv$folds) || p$cv$folds < 2) msg <- c(msg, "cv$folds debe ser >= 2")
  msg
}

# ------------------------------------------------------------
# Funciones estadísticas (sin dependencias adicionales)
# ------------------------------------------------------------

auc_rank <- function(pos, neg) {
  n1 <- as.numeric(length(pos)); n0 <- as.numeric(length(neg))   # numeric: n1 * n0 desborda enteros con fondos grandes
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Índice continuo de Boyce (Hirzel et al. 2006) sobre la escala de rangos del fondo de evaluación
boyce_cont <- function(score_p, score_d, nclass = 101L, window = 0.1) {
  Fd <- stats::ecdf(score_d)
  p <- Fd(score_p); b <- Fd(score_d)
  mids <- seq(window / 2, 1 - window / 2, length.out = nclass)
  P <- vapply(mids, function(m) mean(p >= m - window / 2 & p <= m + window / 2), numeric(1))
  E <- vapply(mids, function(m) mean(b >= m - window / 2 & b <= m + window / 2), numeric(1))
  pe <- ifelse(E > 0, P / E, NA_real_)
  ok <- is.finite(pe)
  idx <- if (sum(ok) > 5L) suppressWarnings(stats::cor(mids[ok], pe[ok], method = "spearman")) else NA_real_
  list(index = idx, curve = tibble::tibble(mid = mids, pe = pe))
}

# Muestreo por pares de distancias (Hijmans 2012): cada presencia de prueba se empareja con el fondo de prueba
# cuya distancia al conjunto de entrenamiento es la más parecida
pds_idx <- function(d_p, d_d, seed = 1L, tol = 0.33) {
  # Hijmans (2012): si las dos distancias difieren más de 33 % (aquí relativo a la mayor de ambas) la presencia no se usa
  n <- length(d_d); ob <- order(d_d); srt <- d_d[ob]; avail <- rep(TRUE, n)
  set.seed(seed); ip <- integer(length(d_p)); idd <- integer(length(d_p)); k <- 0L
  for (i in sample.int(length(d_p))) {
    pos <- findInterval(d_p[i], srt); lo <- pos; hi <- pos + 1L
    while (lo >= 1L && !avail[lo]) lo <- lo - 1L
    while (hi <= n && !avail[hi]) hi <- hi + 1L
    pick <- if (lo < 1L) hi else if (hi > n) lo else if ((d_p[i] - srt[lo]) <= (srt[hi] - d_p[i])) lo else hi
    if (pick < 1L || pick > n) next
    if (abs(d_p[i] - srt[pick]) > tol * max(d_p[i], srt[pick])) next
    avail[pick] <- FALSE; k <- k + 1L; ip[k] <- i; idd[k] <- ob[pick]
  }
  list(p = ip[seq_len(k)], d = idd[seq_len(k)])
}

wilson_ci <- function(k, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  ph <- k / n; den <- 1 + z^2 / n
  ctr <- (ph + z^2 / (2 * n)) / den
  hw <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / den
  c(max(0, ctr - hw), min(1, ctr + hw))
}

# Verosimilitud IPP: sum_i x_i' b - c * sum_j exp(x_j' b); Newton-Raphson (cóncava). c = área por punto de fondo (km2)
ipp_fit <- function(Xp, Xd, cc, beta = NULL, maxit = 60L, tol = 1e-9) {
  k <- ncol(Xp); sp <- colSums(Xp)
  if (is.null(beta)) { beta <- stats::setNames(rep(0, k), colnames(Xp)); beta[1] <- log(nrow(Xp) / (cc * nrow(Xd))) }
  ll_fun <- function(b) sum(Xp %*% b) - cc * sum(exp(Xd %*% b))
  ll <- ll_fun(beta); conv <- FALSE; it <- 0L
  for (it in seq_len(maxit)) {
    mu <- cc * exp(drop(Xd %*% beta))
    g <- sp - drop(crossprod(Xd, mu))
    H <- crossprod(Xd * mu, Xd)
    step <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(step)) break
    t <- 1; ok <- FALSE
    for (h in 1:25) {
      bn <- beta + t * step; lln <- ll_fun(bn)
      if (is.finite(lln) && lln >= ll - 1e-12) { ok <- TRUE; break }
      t <- t / 2
    }
    if (!ok) break
    d <- abs(lln - ll); beta <- bn; ll <- lln
    if (d < tol * (1 + abs(ll)) || max(abs(t * step)) < 1e-8) { conv <- TRUE; break }
  }
  mu <- cc * exp(drop(Xd %*% beta))
  V <- tryCatch(solve(crossprod(Xd * mu, Xd)), error = function(e) matrix(NA_real_, k, k))
  list(beta = beta, ll = ll, vcov = V, se = sqrt(diag(V)), converged = conv, iter = it)
}

morans_correlogram <- function(xy, r, band_m = 1000, nband = 20L, nperm = 99L, seed = 1L) {
  n <- length(r); z <- r - mean(r); s2 <- sum(z^2)
  D <- as.matrix(stats::dist(xy))
  set.seed(seed)
  Zp <- vapply(seq_len(nperm), function(i) sample(z), numeric(n))
  out <- lapply(seq_len(nband), function(b) {
    W <- (D > (b - 1) * band_m & D <= b * band_m) * 1
    S0 <- sum(W); if (S0 == 0) return(tibble::tibble(band = b, dist_km = b * band_m / 1000, moran_i = NA_real_, p_value = NA_real_, n_pairs = 0))
    I <- (n / S0) * drop(crossprod(z, W %*% z)) / s2
    Ip <- (n / S0) * colSums((W %*% Zp) * Zp) / s2
    tibble::tibble(band = b, dist_km = b * band_m / 1000, moran_i = I, p_value = (1 + sum(Ip >= I)) / (nperm + 1), n_pairs = S0)
  })
  dplyr::bind_rows(out)
}

# ------------------------------------------------------------
# Módulo
# ------------------------------------------------------------

run_method <- function(ctx) {
  p <- ctx$params; seed <- ctx$seed; crs <- ctx$crs
  t_start <- Sys.time()
  say <- function(...) cat(sprintf("[M03 %s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  path <- function(x) if (grepl("^(/|[A-Za-z]:)", x)) x else file.path(ctx$root, x)
  rd <- function(file, layer = NULL) {
    v <- if (is.null(layer)) sf::st_read(path(file), quiet = TRUE) else sf::st_read(path(file), layer = layer, quiet = TRUE)
    # solo los poligonos se reparan: st_make_valid() puede convertir lineas en GEOMETRYCOLLECTION y terra las descartaria
    if (any(as.character(sf::st_geometry_type(v)) %in% c("POLYGON", "MULTIPOLYGON"))) v <- sf::st_make_valid(v)
    sf::st_transform(v, crs)
  }
  dat <- p$data; scale_test <- p$testing_scale %||% 1
  ratio_main <- p$background$ratio_main * scale_test
  ladder <- sort(unique(p$background$ladder * scale_test))
  seeds <- as.integer(p$background$seeds)
  fsets <- p$feature_sets; final_fs <- clean_chr(p$final_feature_set)
  allvars <- unique(unlist(fsets))
  n_rep <- as.integer(p$cv$repeats %||% 3L); k_cv <- as.integer(p$cv$folds)
  thr_road <- as.numeric(p$thresholds$road_distance_m %||% 250)
  # umbrales de sensibilidad normativos (derecho de via 25 m desde el eje; 30 m con el retiro adicional de 5 m): Reglamento de la Ley de Caminos, Art. 4
  thr_road_sens <- as.numeric(unlist(p$thresholds$road_distance_sensitivity_m))
  # proxies de la red vial estatal (el derecho de vía de 25 m + 5 m del reglamento de la Ley de Caminos aplica a la red estatal, no a todas las vías)
  state_nets <- p$state_networks
  slope_thr <- p$thresholds$slope_deg
  min_n <- as.integer(p$min_presences_subset %||% 200L)
  dom_main <- clean_chr(p$domains$main)
  dom_sens <- as.character(unlist(p$domains$sensitivity))

  # ----------------------------------------------------------
  # 1. Dominios, presencias y fondo (con caché de covariables)
  # ----------------------------------------------------------
  say("Leyendo capas y construyendo dominios...")
  excl <- sf::st_union(rd(dat$exclude_vector))
  mz <- sf::st_union(rd(dat$census_gpkg, dat$blocks_layer))
  bad <- sf::st_union(excl, mz)
  aoi_lyr <- function(nm) sf::st_union(rd(dat$aoi_gpkg, nm))
  aoi_amp <- aoi_lyr("aoi_amplio")
  ext_pre <- sf::st_buffer(aoi_amp, as.numeric(p$domains$extended_buffer_m %||% 5000))
  zu <- rd(dat$igm_gpkg, dat$core_layer)
  cores <- sf::st_union(zu[zu$nam %in% as.character(unlist(dat$core_names)), ])
  core_pre <- sf::st_intersection(sf::st_buffer(cores, as.numeric(p$domains$core_buffer_m %||% 10000)), ext_pre)
  dom_def <- list(
    aoi_estricto = sf::st_difference(aoi_lyr("aoi_estricto"), bad),
    aoi_base = sf::st_difference(aoi_lyr("aoi_base"), bad),
    aoi_amplio = sf::st_difference(aoi_amp, bad),
    extended = sf::st_difference(ext_pre, bad),
    core_buffer = sf::st_difference(core_pre, bad)
  )
  dom_def$extended <- sf::st_union(dom_def$extended, dom_def$aoi_amplio)
  # todos los dominios se limitan a la extensión del ráster de pendiente (sin pendiente no hay covariables)
  re <- as.vector(terra::ext(terra::rast(path(dat$slope_raster))))
  raster_box <- sf::st_as_sfc(sf::st_bbox(c(xmin = re[[1]], xmax = re[[2]], ymin = re[[3]], ymax = re[[4]]), crs = sf::st_crs(crs)))
  dom_def <- lapply(dom_def, function(g) sf::st_intersection(g, raster_box))
  dom_area <- vapply(dom_def, function(g) as.numeric(sf::st_area(g)) / 1e6, numeric(1))
  if (!dom_main %in% names(dom_def)) stop("domains$main desconocido: ", dom_main, call. = FALSE)
  say("Áreas (km2): ", paste(names(dom_area), round(dom_area), collapse = "; "))

  cache_key <- substr(hash_object(list(dat = dat, dom = p$domains, bg = list(p$background, scale_test), min_n = min_n, allvars = allvars, st = state_nets, code_version = "M03-2026-10-21c")), 1, 10)
  cache_file <- file.path(ctx$DIR$processed, paste0("M03_dataset_", cache_key, ".rds"))

  if (file.exists(cache_file)) {
    say("Cargando caché de covariables: ", basename(cache_file))
    cache <- readRDS(cache_file)
  } else {
    pts_all <- rd(dat$census_gpkg, dat$presence_layer)
    pts <- pts_all[pts_all[[dat$in_block_field %||% "en_manzana"]] == 0, ]
    n_raw <- nrow(pts)
    in_ext <- lengths(sf::st_intersects(pts, dom_def$extended)) > 0
    pts <- pts[in_ext, ]
    say("Presencias (edificios fuera de manzanas, fuera de Podocarpus): ", nrow(pts), " de ", n_raw)

    # fondo: muestreo uniforme en el dominio más amplio; las bolsas se anidan por dominio
    n_main_p <- sum(lengths(sf::st_intersects(pts, dom_def[[dom_main]])) > 0)
    pool_n <- ceiling(max(ladder) * n_main_p * dom_area[["extended"]] / dom_area[[dom_main]] * 1.03)
    say("Generando fondo uniforme: ", pool_n, " puntos en el dominio extendido (semilla ", seeds[1], ")...")
    set.seed(seeds[1]); pool <- sf::st_sample(dom_def$extended, size = pool_n, type = "random", exact = TRUE)
    pool <- sf::st_sf(geometry = pool[sample.int(length(pool))])
    seed_pools <- list()
    for (s in seeds[-1]) {
      set.seed(s); n_s <- ceiling(ratio_main * n_main_p)
      sp <- sf::st_sample(dom_def[[dom_main]], size = n_s, type = "random", exact = TRUE)
      seed_pools[[as.character(s)]] <- sf::st_sf(geometry = sp[sample.int(length(sp))])
    }

    say("Cargando predictores...")
    slope_r <- terra::rast(path(dat$slope_raster))
    ex_r <- as.vector(terra::ext(slope_r))   # xmin, xmax, ymin, ymax
    lines_of <- function(gpkg, layers) do.call(c, lapply(layers, function(l) sf::st_geometry(rd(gpkg, l))))
    L <- list(
      igm_vias = lines_of(dat$igm_gpkg, as.character(unlist(dat$igm_roads_layers))),
      igm_all = lines_of(dat$igm_gpkg, as.character(unlist(dat$igm_all_layers))),
      osm2019 = lines_of(dat$osm_gpkg, dat$osm_layer),
      census2020 = lines_of(dat$census_gpkg, dat$census_axes_layer),
      river = lines_of(dat$igm_gpkg, as.character(unlist(dat$river_layers))),
      nucleo = sf::st_geometry(zu[zu$nam %in% as.character(unlist(dat$core_names)), ])
    )
    lines_where <- function(gpkg, layer, field, values) { x <- rd(gpkg, layer); sf::st_geometry(x[x[[field]] %in% values, ]) }
    for (nm in names(state_nets)) {
      sn <- state_nets[[nm]]; vals <- as.character(unlist(sn$values))
      L[[paste0("state_", nm)]] <- if (identical(sn$source, "igm")) lines_where(dat$igm_gpkg, as.character(sn$layer), sn$field, vals) else lines_where(dat$osm_gpkg, dat$osm_layer, sn$field, vals)
    }
    L <- lapply(L, function(l) l[!sf::st_is_empty(l)])   # terra descarta geometrías vacías y el emparejamiento por pares fallaría
    centres_poly <- rd(dat$aoi_gpkg, "aoi_amplio")
    # distancia exacta al elemento mas cercano (misma cifra que sf::st_distance(by_element = TRUE), ~7 veces mas rapido; verificado con 20 000 puntos)
    dist_to <- function(g, lines) {
      idx <- sf::st_nearest_feature(g, lines)
      as.numeric(terra::distance(terra::vect(g), terra::vect(lines[idx]), pairwise = TRUE))
    }
    covariates <- function(sfobj, label) {
      g <- sf::st_geometry(sfobj)
      say("Covariables de ", label, " (", length(g), " puntos)...")
      xy <- sf::st_coordinates(g)
      inside <- xy[, 1] >= ex_r[1] & xy[, 1] <= ex_r[2] & xy[, 2] >= ex_r[3] & xy[, 2] <= ex_r[4]
      sl <- rep(NA_real_, length(g))
      sl[inside] <- terra::extract(slope_r, xy[inside, , drop = FALSE])[, 1] / (dat$slope_raster_divisor %||% 2)
      n_na_slope <- sum(inside & is.na(sl))
      sl[inside & is.na(sl)] <- 0     # el ráster usa 0 como NoData: pendiente < 0,25 grados
      out <- tibble::tibble(x = xy[, 1], y = xy[, 2], slope = sl,
                            d_igm_vias = dist_to(g, L$igm_vias), d_igm_all = dist_to(g, L$igm_all),
                            d_osm2019 = dist_to(g, L$osm2019), d_census2020 = dist_to(g, L$census2020),
                            d_nucleo = dist_to(g, L$nucleo), d_river = dist_to(g, L$river))
      for (nm in names(state_nets)) out[[paste0("d_state_", nm)]] <- dist_to(g, L[[paste0("state_", nm)]])
      ci <- sf::st_intersects(g, centres_poly)
      out$centre <- vapply(ci, function(i) if (length(i) == 0L) "fuera_AOI" else as.character(centres_poly$centro[i[1]]), character(1))
      for (nm in grep("^d_", names(out), value = TRUE)) out[[paste0("log_", nm)]] <- log1p(out[[nm]])
      attr(out, "n_na_slope_inside") <- n_na_slope
      attr(out, "n_outside_slope_raster") <- sum(!inside)
      out
    }
    flags <- function(g) {
      out <- lapply(names(dom_def), function(d) lengths(sf::st_intersects(g, dom_def[[d]])) > 0)
      names(out) <- paste0("in_", names(dom_def)); tibble::as_tibble(out)
    }
    cp <- covariates(pts, "presencias"); cp <- dplyr::bind_cols(cp, flags(sf::st_geometry(pts)))
    cp$centre_census <- as.character(pts$centro_aoi); cp$centre_census[is.na(cp$centre_census) | !nzchar(cp$centre_census)] <- "fuera_AOI"
    cd <- covariates(pool, "fondo principal")
    n_na_slope <- attr(cd, "n_na_slope_inside"); n_outside_slope <- attr(cd, "n_outside_slope_raster")
    cd <- dplyr::bind_cols(cd, flags(sf::st_geometry(pool)))
    cs <- lapply(seed_pools, function(sp) covariates(sp, "fondo de semilla"))
    cache <- list(pres = cp, pool = cd, seed_pools = cs, n_raw_presence = n_raw,
                  n_na_slope = n_na_slope, n_outside_slope = n_outside_slope)
    saveRDS(cache, cache_file)
    say("Caché guardada: ", basename(cache_file))
  }

  P0 <- cache$pres; D0 <- cache$pool
  bad_slope <- sum(!is.finite(P0$slope)) + sum(!is.finite(D0$slope))
  if (bad_slope > 0L) {
    warning("Se descartan ", bad_slope, " puntos sin pendiente (fuera del ráster).", call. = FALSE)
    P0 <- P0[is.finite(P0$slope), ]; D0 <- D0[is.finite(D0$slope), ]
  }

  dom_data <- function(dm, pool = D0, n_back = NULL) {
    fl <- paste0("in_", dm)
    P <- P0[P0[[fl]], ]; D <- pool[pool[[fl]], ]
    if (!is.null(n_back)) {
      if (nrow(D) < n_back) { warning("Fondo insuficiente en ", dm, ": ", nrow(D), " < ", n_back, call. = FALSE) } else D <- D[seq_len(n_back), ]
    }
    list(P = P, D = D)
  }

  main <- dom_data(dom_main, n_back = ceiling(ratio_main * sum(P0[[paste0("in_", dom_main)]])))
  nP <- nrow(main$P); nD <- nrow(main$D)
  cc <- dom_area[[dom_main]] / nD
  say("Modelo principal: ", nP, " presencias, ", nD, " puntos de fondo; c = ", signif(cc, 4), " km2 por punto")
  zs_of <- function(D) list(mean = vapply(c("slope", paste0("log_d_", c("igm_vias", "igm_all", "osm2019", "census2020", "nucleo", "river"))), function(v) mean(D[[v]]), numeric(1)),
                            sd = vapply(c("slope", paste0("log_d_", c("igm_vias", "igm_all", "osm2019", "census2020", "nucleo", "river"))), function(v) stats::sd(D[[v]]), numeric(1)))
  zs <- zs_of(main$D)
  add_z <- function(T, zsx) { for (v in names(zsx$mean)) T[[paste0("z_", v)]] <- (T[[v]] - zsx$mean[[v]]) / zsx$sd[[v]]; T }
  mk_X <- function(T, vars, ns_df = 0L, bases = NULL) {
    cols <- lapply(vars, function(v) {
      z <- T[[paste0("z_", v)]]
      if (ns_df > 0L) { m <- stats::predict(bases[[v]], z); colnames(m) <- paste0(v, "_ns", seq_len(ncol(m))); m }
      else matrix(z, ncol = 1L, dimnames = list(NULL, v))
    })
    cbind(`(Intercept)` = 1, do.call(cbind, cols))
  }
  fit_set <- function(P, D, vars, cc, zsx, ns_df = 0L, bases = NULL, beta = NULL) {
    Pz <- add_z(P, zsx); Dz <- add_z(D, zsx)
    ipp_fit(mk_X(Pz, vars, ns_df, bases), mk_X(Dz, vars, ns_df, bases), cc, beta)
  }

  # ----------------------------------------------------------
  # 2. Colinealidad (Dormann et al. 2013)
  # ----------------------------------------------------------
  Dz <- add_z(main$D, zs)
  # bases de splines naturales (modelo principal de ajuste y validación; el lineal resume el efecto promedio para las hipótesis)
  ns_df <- as.integer(p$nonlinear_df %||% 3L)
  bases <- if (ns_df > 0L) lapply(stats::setNames(allvars, allvars), function(v) splines::ns(Dz[[paste0("z_", v)]], df = ns_df)) else NULL
  coll <- purrr::imap_dfr(fsets, function(vars, nm) {
    vars <- as.character(unlist(vars)); cm <- stats::cor(Dz[, paste0("z_", vars), drop = FALSE])
    vif <- if (length(vars) > 1L) diag(solve(cm)) else 1
    tibble::tibble(feature_set = nm, variable = vars, vif = as.numeric(vif),
                   max_abs_r_with_others = vapply(seq_along(vars), function(i) if (length(vars) > 1L) max(abs(cm[i, -i])) else 0, numeric(1)))
  })
  coll_pairs <- {
    cm <- stats::cor(Dz[, paste0("z_", allvars), drop = FALSE]); dimnames(cm) <- list(M03_VAR_LABEL[allvars], M03_VAR_LABEL[allvars])
    as.data.frame(as.table(cm)) |> stats::setNames(c("var1", "var2", "pearson_r")) |> dplyr::filter(as.character(var1) < as.character(var2)) |>
      dplyr::mutate(exceeds_0.7 = abs(pearson_r) > 0.7)
  }

  # ----------------------------------------------------------
  # 3. Ajuste del modelo principal por conjunto de variables y convergencia del fondo
  # ----------------------------------------------------------
  say("Ajustando conjuntos de variables en el dominio principal...")
  fits <- lapply(fsets, function(vars) fit_set(main$P, main$D, as.character(unlist(vars)), cc, zs))
  fits_spl <- if (ns_df > 0L) lapply(fsets, function(vars) fit_set(main$P, main$D, as.character(unlist(vars)), cc, zs, ns_df, bases)) else list()
  not_conv <- c(names(fits)[!vapply(fits, function(f) isTRUE(f$converged), logical(1))],
                names(fits_spl)[!vapply(fits_spl, function(f) isTRUE(f$converged), logical(1))])
  if (length(not_conv) > 0L) warning("Newton-Raphson no convergió en: ", paste(not_conv, collapse = ", "), call. = FALSE)

  coef_rows <- function(f, vars, subset, lo = NULL, hi = NULL, prefix = "intensity_ratio_", note = "") {
    b <- f$beta[vars]; se <- f$se[vars]
    lo <- lo %||% exp(b - 1.96 * se); hi <- hi %||% exp(b + 1.96 * se)
    dplyr::bind_rows(lapply(seq_along(vars), function(i) metric_rows(paste0(prefix, M03_VAR_LABEL[[vars[i]]]), exp(b[[i]]), subset, lo[[i]], hi[[i]], "ratio per SD", note)))
  }
  wald_note <- "95% interval from model-based (Poisson) standard errors; underestimates uncertainty under spatial clustering"

  # convergencia del fondo (Renner et al. 2015)
  say("Convergencia del tamaño del fondo...")
  vars_fin <- as.character(unlist(fsets[[final_fs]]))
  conv_tab <- purrr::map_dfr(ladder, function(r) {
    d <- dom_data(dom_main, n_back = ceiling(r * nP))
    cc_r <- dom_area[[dom_main]] / nrow(d$D)
    f <- fit_set(d$P, d$D, vars_fin, cc_r, zs)
    # error de Monte Carlo de la integral del fondo (= error de la log-verosimilitud): con puntos aleatorios la log-verosimilitud fluctúa sin tendencia
    mu_r <- cc_r * exp(drop(mk_X(add_z(d$D, zs), vars_fin) %*% f$beta))
    tibble::tibble(ratio = r, n_background = nrow(d$D), log_likelihood = f$ll, ll_mc_se = sqrt(nrow(d$D)) * stats::sd(mu_r), converged = f$converged,
                   variable = vars_fin, beta = f$beta[vars_fin], se = f$se[vars_fin], intensity_ratio = exp(f$beta[vars_fin]))
  })
  ll_by_ratio <- conv_tab |> dplyr::distinct(ratio, n_background, log_likelihood, ll_mc_se) |> dplyr::arrange(ratio) |>
    dplyr::mutate(delta_ll_vs_previous = log_likelihood - dplyr::lag(log_likelihood),
                  delta_mc_se = sqrt(ll_mc_se^2 + dplyr::lag(ll_mc_se)^2),
                  delta_within_2_mc_se = abs(delta_ll_vs_previous) < 2 * delta_mc_se)
  big <- max(ladder); prev <- sort(ladder)[length(ladder) - 1L]
  max_rel_change_ir <- conv_tab |> dplyr::filter(ratio %in% c(big, prev)) |> dplyr::group_by(variable) |>
    dplyr::summarise(rel = abs(diff(intensity_ratio)) / dplyr::first(intensity_ratio[order(ratio)]), .groups = "drop")
  # criterio: razones estables (< 2 %) y cambio de log-verosimilitud compatible con el error de Monte Carlo (Renner et al. 2015 usan ΔLL < 2 con malla regular;
  # con puntos aleatorios la log-verosimilitud tiene ruido de muestreo, por eso se compara con 2 errores estándar de la diferencia)
  convergence_ok <- all(max_rel_change_ir$rel < 0.02) && isTRUE(utils::tail(ll_by_ratio$delta_within_2_mc_se, 1))
  say("Convergencia del fondo (último paso): cambio máx. en IR = ", round(100 * max(max_rel_change_ir$rel), 2), "%; ΔLL = ",
      round(utils::tail(ll_by_ratio$delta_ll_vs_previous, 1), 2), " (2 EE de Monte Carlo = ", round(2 * utils::tail(ll_by_ratio$delta_mc_se, 1), 1), ") -> ", if (convergence_ok) "OK" else "NO CONVERGE")

  # logit ingenuo 1:3 sin ponderar (para mostrar el sesgo señalado por Fithian & Hastie 2013)
  r_naive <- as.numeric(p$background$naive_logit_ratio %||% 3) * scale_test
  dn <- dom_data(dom_main, n_back = ceiling(r_naive * nP))
  naive_d <- rbind(add_z(dn$P, zs)[, paste0("z_", vars_fin)] |> dplyr::mutate(y = 1), add_z(dn$D, zs)[, paste0("z_", vars_fin)] |> dplyr::mutate(y = 0))
  naive <- stats::glm(stats::reformulate(paste0("z_", vars_fin), "y"), data = naive_d, family = stats::binomial())
  naive_ci <- suppressMessages(stats::confint.default(naive))

  # ----------------------------------------------------------
  # 4. Autocorrelación residual y tamaño de bloque (Roberts 2017; Valavi 2019)
  # ----------------------------------------------------------
  say("Correlograma de residuos y tamaño de bloque...")
  f_fin <- fits[[final_fs]]
  # residuos del modelo principal de ajuste (splines si están activos; si no, el lineal)
  eta_d <- if (ns_df > 0L) drop(mk_X(Dz, vars_fin, ns_df, bases) %*% fits_spl[[final_fs]]$beta) else drop(mk_X(Dz, vars_fin) %*% f_fin$beta)
  cell_m <- 1000
  cellid <- function(x, y) paste(floor(x / cell_m), floor(y / cell_m))
  E_cell <- tapply(cc * exp(eta_d), cellid(main$D$x, main$D$y), sum)
  O_cell <- table(factor(cellid(main$P$x, main$P$y), levels = names(E_cell)))
  cells <- tibble::tibble(cell = names(E_cell), E = as.numeric(E_cell), O = as.integer(O_cell[names(E_cell)])) |>
    dplyr::mutate(cx = (as.numeric(sub(" .*", "", cell)) + 0.5) * cell_m, cy = (as.numeric(sub(".* ", "", cell)) + 0.5) * cell_m,
                  resid = (O - E) / sqrt(E)) |> dplyr::filter(E >= 0.5)
  corr <- morans_correlogram(cbind(cells$cx, cells$cy), cells$resid, band_m = cell_m, nband = as.integer(p$cv$correlogram_bands %||% 20L), seed = seed)
  first_ns <- which(!(corr$p_value < 0.05 & corr$moran_i > 0))[1]
  range_m <- if (is.na(first_ns)) max(corr$dist_km) * 1000 else (first_ns - 1L) * cell_m
  bnd <- as.numeric(unlist(p$cv$block_auto_bounds_m %||% c(2000, 16000)))
  block_m <- if (identical(clean_chr(p$cv$block_size_m), "auto")) min(max(ceiling(range_m / 1000) * 1000, bnd[1]), bnd[2]) else as.numeric(p$cv$block_size_m)
  say("Rango de autocorrelación residual ≈ ", range_m, " m -> bloque = ", block_m, " m")

  # ----------------------------------------------------------
  # 5. Validación cruzada espacial (+ aleatoria como referencia)
  # ----------------------------------------------------------
  AL_of <- function(d) {
    A <- dplyr::bind_rows(d$P |> dplyr::mutate(yy = 1L), d$D |> dplyr::mutate(yy = 0L))
    add_z(A, zs)
  }
  run_split <- function(A, tr, te, setlist, ns_df = 0L, bases = NULL, seed_i = 1L) {
    trp <- tr & A$yy == 1L; trd <- tr & A$yy == 0L; tep <- te & A$yy == 1L; ted <- te & A$yy == 0L
    if (sum(tep) < 5L || sum(ted) < 5L || sum(trp) < 20L) return(NULL)
    g_tr <- sf::st_as_sf(A[trp, c("x", "y")], coords = c("x", "y"), crs = crs)
    near <- function(idx) { g <- sf::st_as_sf(A[idx, c("x", "y")], coords = c("x", "y"), crs = crs)
      as.numeric(sf::st_distance(g, g_tr[sf::st_nearest_feature(g, g_tr), ], by_element = TRUE)) }
    d_p <- near(tep); d_d <- near(ted)
    pr <- pds_idx(d_p, d_d, seed_i)
    nauc <- auc_rank(-d_p, -d_d); nauc_pds <- auc_rank(-d_p[pr$p], -d_d[pr$d])
    ssb <- mean(d_p) / mean(d_d); ssb_pds <- mean(d_p[pr$p]) / mean(d_d[pr$d])   # Dp / Da (Hijmans 2012): 1 = sin sesgo de ordenamiento
    purrr::imap_dfr(setlist, function(vars, nm) {
      vars <- as.character(unlist(vars))
      f <- ipp_fit(mk_X(A[trp, ], vars, ns_df, bases), mk_X(A[trd, ], vars, ns_df, bases), cc)
      bt <- f$beta; ep <- drop(mk_X(A[tep, ], vars, ns_df, bases) %*% bt); ed <- drop(mk_X(A[ted, ], vars, ns_df, bases) %*% bt)
      a <- auc_rank(ep, ed); bc <- boyce_cont(ep, ed)
      a_pds <- auc_rank(ep[pr$p], ed[pr$d])
      tibble::tibble(feature_set = nm, auc = a, auc_pds = a_pds, auc_null = nauc, auc_null_pds = nauc_pds,
                     cauc = a + 0.5 - max(0.5, nauc), cauc_pds = a_pds + 0.5 - max(0.5, nauc_pds), boyce = bc$index,
                     ssb = ssb, ssb_pds = ssb_pds, n_pds_pairs = length(pr$p),
                     cv_ll = sum(ep) - cc * sum(exp(ed)),   # log-verosimilitud IPP fuera de muestra (no depende de la independencia de los datos de prueba)
                     n_test_p = sum(tep), n_test_d = sum(ted), converged = f$converged)
    })
  }
  cv_run <- function(A, setlist, S, reps, ns_df = 0L, bases = NULL, random = FALSE) {
    purrr::map_dfr(seq_len(reps), function(r) {
      fo <- if (random) { set.seed(seed + 100L + r); sample(rep_len(seq_len(k_cv), nrow(A))) } else spatial_block_folds(A$x, A$y, S, k_cv, seed + r)
      purrr::map_dfr(seq_len(k_cv), function(f) {
        res <- run_split(A, fo != f, fo == f, setlist, ns_df, bases, seed_i = seed + r * 10L + f)
        if (is.null(res)) NULL else dplyr::mutate(res, repeat_id = r, fold = f)
      })
    })
  }
  A_main <- AL_of(main)
  say("Validación cruzada espacial (bloque ", block_m, " m, ", k_cv, " pliegues x ", n_rep, " repeticiones)...")
  cv_sp <- cv_run(A_main, fsets, block_m, n_rep)
  say("Validación cruzada aleatoria (referencia)...")
  cv_rnd <- cv_run(A_main, fsets, block_m, 1L, random = TRUE)
  summ_cv <- function(cv, prefix_note) {
    cv |> dplyr::group_by(feature_set) |> dplyr::summarise(
      dplyr::across(c(auc, auc_pds, auc_null, cauc, cauc_pds, boyce, ssb, ssb_pds, cv_ll), list(mean = ~mean(.x, na.rm = TRUE), lo = ~min(.x, na.rm = TRUE), hi = ~max(.x, na.rm = TRUE)), .names = "{.col}__{.fn}"),
      n_splits = dplyr::n(), .groups = "drop")
  }
  s_sp <- summ_cv(cv_sp); s_rnd <- summ_cv(cv_rnd)
  cv_rows <- function(s, sub_fun, which, nm_map, note) {
    dplyr::bind_rows(lapply(seq_len(nrow(s)), function(i) dplyr::bind_rows(lapply(names(nm_map), function(col)
      metric_rows(nm_map[[col]], s[[paste0(col, "__mean")]][i], sub_fun(s$feature_set[i]), s[[paste0(col, "__lo")]][i], s[[paste0(col, "__hi")]][i], "", note)))))
  }
  note_cv <- paste0("Mean over folds; interval = min-max across ", k_cv, "-fold x ", n_rep, " spatial block splits (not a confidence interval)")
  m_cv <- cv_rows(s_sp, identity, "sp", c(auc = "auc_spatial_cv", auc_pds = "auc_pds_spatial_cv", auc_null = "auc_null_distance", cauc = "cauc_spatial_cv", cauc_pds = "cauc_pds_spatial_cv", boyce = "boyce_spatial_cv", ssb = "spatial_sorting_bias", ssb_pds = "spatial_sorting_bias_after_pds", cv_ll = "cv_log_likelihood"), note_cv)
  # modelo principal de ajuste: splines naturales (subconjuntos spl_<conjunto>)
  if (ns_df > 0L) {
    say("Validación cruzada espacial del modelo con splines...")
    cv_spl <- cv_run(A_main, fsets, block_m, n_rep, ns_df, bases) |> dplyr::mutate(feature_set = paste0("spl_", feature_set))
    s_spl <- summ_cv(cv_spl)
    m_cv <- dplyr::bind_rows(m_cv, cv_rows(s_spl, identity, "sp", c(auc = "auc_spatial_cv", auc_pds = "auc_pds_spatial_cv", auc_null = "auc_null_distance", cauc = "cauc_spatial_cv", cauc_pds = "cauc_pds_spatial_cv", boyce = "boyce_spatial_cv", ssb = "spatial_sorting_bias", ssb_pds = "spatial_sorting_bias_after_pds", cv_ll = "cv_log_likelihood"), note_cv))
    s_sp <- dplyr::bind_rows(s_sp, s_spl); cv_sp <- dplyr::bind_rows(cv_sp, cv_spl)
  }
  m_rnd <- cv_rows(s_rnd, identity, "rnd", c(auc = "auc_random_cv"), "Mean over folds of non-spatial random k-fold CV; shown to quantify optimism (Roberts et al. 2017; Pohjankukka et al. 2017)")

  # curva P/E agrupada del modelo final (todos los pliegues, una repetición)
  fo1 <- spatial_block_folds(A_main$x, A_main$y, block_m, k_cv, seed + 1L)
  pe_curve <- purrr::map_dfr(seq_len(k_cv), function(f) {
    tr <- fo1 != f; te <- fo1 == f
    trp <- tr & A_main$yy == 1L; trd <- tr & A_main$yy == 0L
    fx <- ipp_fit(mk_X(A_main[trp, ], vars_fin, ns_df, bases), mk_X(A_main[trd, ], vars_fin, ns_df, bases), cc)
    ep <- drop(mk_X(A_main[te & A_main$yy == 1L, ], vars_fin, ns_df, bases) %*% fx$beta); ed <- drop(mk_X(A_main[te & A_main$yy == 0L, ], vars_fin, ns_df, bases) %*% fx$beta)
    boyce_cont(ep, ed)$curve |> dplyr::mutate(fold = f)
  })

  # ----------------------------------------------------------
  # 6. Efectos: razones de intensidad con bootstrap por bloques
  # ----------------------------------------------------------
  B <- as.integer(p$bootstrap$B %||% 300L)
  boot_sets <- as.character(unlist(p$bootstrap$sets %||% names(fsets)))
  blk_p <- cellid_block <- NULL
  blk <- function(T) paste(floor(T$x / block_m), floor(T$y / block_m))
  boot_fit <- function(Xp, Xd, bp, bd, beta0, B, sd0) {
    ub <- union(unique(bp), unique(bd)); spp <- split(seq_len(nrow(Xp)), factor(bp, levels = ub)); sdd <- split(seq_len(nrow(Xd)), factor(bd, levels = ub))
    out <- matrix(NA_real_, B, ncol(Xp), dimnames = list(NULL, colnames(Xp))); set.seed(sd0)
    for (b in seq_len(B)) {
      s <- sample(ub, length(ub), replace = TRUE)
      ip <- unlist(spp[s], use.names = FALSE); id <- unlist(sdd[s], use.names = FALSE)
      f <- tryCatch(ipp_fit(Xp[ip, , drop = FALSE], Xd[id, , drop = FALSE], cc, beta0, maxit = 30L), error = function(e) NULL)
      if (!is.null(f) && isTRUE(f$converged)) out[b, ] <- f$beta
    }
    out
  }
  say("Bootstrap por bloques (B = ", B, ")...")
  Pz <- add_z(main$P, zs)
  boot_res <- lapply(stats::setNames(boot_sets, boot_sets), function(nm) {
    vars <- as.character(unlist(fsets[[nm]]))
    boot_fit(mk_X(Pz, vars), mk_X(Dz, vars), blk(main$P), blk(main$D), fits[[nm]]$beta, B, seed + 7L)
  })
  # bootstrap por bloques del modelo con splines (conjunto final): curvas de respuesta y razones entre cuartiles
  resp_tab <- tibble::tibble(); m_spl_eff <- NULL
  if (ns_df > 0L) {
    say("Bootstrap por bloques del modelo con splines y curvas de respuesta...")
    f_s <- fits_spl[[final_fs]]
    boot_spl <- boot_fit(mk_X(Pz, vars_fin, ns_df, bases), mk_X(Dz, vars_fin, ns_df, bases), blk(main$P), blk(main$D), f_s$beta, B, seed + 9L)
    ok_b <- stats::complete.cases(boot_spl)
    contrib <- function(v, z, beta) { m <- stats::predict(bases[[v]], z); drop(m %*% beta[paste0(v, "_ns", seq_len(ncol(m)))]) }
    resp_tab <- purrr::map_dfr(vars_fin, function(v) {
      zv <- Dz[[paste0("z_", v)]]; zg <- stats::quantile(zv, seq(0.01, 0.99, length.out = 60)); zr <- stats::median(zv)
      est <- contrib(v, zg, f_s$beta) - contrib(v, zr, f_s$beta)
      bb <- t(vapply(which(ok_b), function(i) contrib(v, zg, boot_spl[i, ]) - contrib(v, zr, boot_spl[i, ]), numeric(length(zg))))
      xg <- zg * zs$sd[[v]] + zs$mean[[v]]
      tibble::tibble(variable = M03_VAR_LABEL[[v]], x = if (v == "slope") xg else expm1(xg), x_unit = if (v == "slope") "degrees" else "metres",
                     z_percentile = seq(1, 99, length.out = 60), relative_intensity = exp(est),
                     ci_low = exp(apply(bb, 2, stats::quantile, 0.025)), ci_high = exp(apply(bb, 2, stats::quantile, 0.975)))
    })
    iqr <- purrr::map_dfr(vars_fin, function(v) {
      zv <- Dz[[paste0("z_", v)]]; z25 <- stats::quantile(zv, 0.25); z75 <- stats::quantile(zv, 0.75)
      e <- exp(contrib(v, z75, f_s$beta) - contrib(v, z25, f_s$beta))
      bb <- vapply(which(ok_b), function(i) exp(contrib(v, z75, boot_spl[i, ]) - contrib(v, z25, boot_spl[i, ])), numeric(1))
      tibble::tibble(variable = v, ratio = e, lo = stats::quantile(bb, 0.025), hi = stats::quantile(bb, 0.975))
    })
    m_spl_eff <- dplyr::bind_rows(lapply(seq_len(nrow(iqr)), function(i) metric_rows(paste0("intensity_ratio_iqr_", M03_VAR_LABEL[[iqr$variable[i]]]), iqr$ratio[i], paste0("spl_", final_fs),
                                                                                       iqr$lo[i], iqr$hi[i], "ratio P25 to P75", paste0("Spline model; ratio of intensity between the 75th and 25th percentile of the predictor in the domain; block-bootstrap 95% interval (B = ", sum(ok_b), ")"))))
  }
  m_ir <- dplyr::bind_rows(lapply(names(fsets), function(nm) {
    vars <- as.character(unlist(fsets[[nm]]))
    if (nm %in% names(boot_res)) {
      bm <- boot_res[[nm]]; q <- apply(exp(bm[, vars, drop = FALSE]), 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      coef_rows(fits[[nm]], vars, nm, lo = q[1, ], hi = q[2, ], note = paste0("95% block-bootstrap percentile interval (B = ", B, ", blocks of ", block_m, " m)"))
    } else coef_rows(fits[[nm]], vars, nm, note = wald_note)
  }))
  coef_tab <- purrr::imap_dfr(fits, function(f, nm) {
    vars <- as.character(unlist(fsets[[nm]]))
    tibble::tibble(feature_set = nm, variable = M03_VAR_LABEL[vars], beta = f$beta[vars], se_model = f$se[vars], intensity_ratio = exp(f$beta[vars]),
                   ir_boot_lo = if (nm %in% names(boot_res)) apply(exp(boot_res[[nm]][, vars, drop = FALSE]), 2, stats::quantile, probs = 0.025, na.rm = TRUE) else NA_real_,
                   ir_boot_hi = if (nm %in% names(boot_res)) apply(exp(boot_res[[nm]][, vars, drop = FALSE]), 2, stats::quantile, probs = 0.975, na.rm = TRUE) else NA_real_,
                   log_likelihood = f$ll, converged = f$converged, n_boot_ok = if (nm %in% names(boot_res)) sum(stats::complete.cases(boot_res[[nm]])) else NA_integer_)
  })

  m_extra <- dplyr::bind_rows(
    dplyr::bind_rows(lapply(names(fits), function(nm) metric_rows("log_likelihood", fits[[nm]]$ll, nm, unit = "", note = "IPP log-likelihood up to a constant; comparable only for the same data"))),
    dplyr::bind_rows(lapply(names(fits), function(nm) metric_rows("aic", -2 * fits[[nm]]$ll + 2 * length(fits[[nm]]$beta), nm, note = "AIC of the IPP likelihood; spatial clustering makes likelihood differences overstate evidence, so rely on cv_log_likelihood"))),
    if (ns_df > 0L) dplyr::bind_rows(lapply(names(fits_spl), function(nm) metric_rows("log_likelihood", fits_spl[[nm]]$ll, paste0("spl_", nm), note = paste0("Natural cubic splines, df = ", ns_df, " per predictor")))),
    if (ns_df > 0L) dplyr::bind_rows(lapply(names(fits_spl), function(nm) metric_rows("aic", -2 * fits_spl[[nm]]$ll + 2 * length(fits_spl[[nm]]$beta), paste0("spl_", nm), note = "AIC of the IPP likelihood; spatial clustering makes likelihood differences overstate evidence, so rely on cv_log_likelihood"))),
    dplyr::bind_rows(lapply(names(fits), function(nm) metric_rows("n_presences", nP, nm, unit = "count"))),
    dplyr::bind_rows(lapply(names(fits), function(nm) metric_rows("n_background", nD, nm, unit = "count"))),
    metric_rows("block_size_m", block_m, final_fs, unit = "m", note = "Spatial block size used in cross-validation and bootstrap"),
    metric_rows("residual_range_m", range_m, final_fs, unit = "m", note = "Distance up to which Moran's I of 1-km Pearson residuals is significantly positive (p < 0.05, 99 permutations)")
  )

  # ----------------------------------------------------------
  # 7. Heterogeneidad: un centro fuera, y anillos de distancia al núcleo
  # ----------------------------------------------------------
  say("Validación dejando un centro fuera...")
  cen_lab <- A_main$centre; cen_lab[A_main$yy == 1L] <- A_main$centre_census[A_main$yy == 1L]
  centres <- names(which(table(cen_lab[A_main$yy == 1L]) >= 30L)); centres <- setdiff(centres, "fuera_AOI")
  A_main$blk <- blk(A_main)
  loco <- purrr::map_dfr(centres, function(cn) {
    te <- cen_lab == cn
    tr <- !te & !(A_main$blk %in% unique(A_main$blk[te]))
    r <- run_split(A_main, tr, te, fsets[final_fs], seed_i = seed)
    if (is.null(r)) NULL else dplyr::mutate(r, centre = cn, n_presences = sum(te & A_main$yy == 1L))
  })
  m_loco <- dplyr::bind_rows(
    if (nrow(loco) > 0L) dplyr::bind_rows(lapply(seq_len(nrow(loco)), function(i) dplyr::bind_rows(
      metric_rows("auc_spatial_cv", loco$auc[i], paste0("center_", loco$centre[i]), unit = "", note = "Trained on other centres' blocks, tested on this centre; includes distance-null calibration in tab_M03_loco"),
      metric_rows("cauc_spatial_cv", loco$cauc[i], paste0("center_", loco$centre[i]), unit = "", note = "Hijmans (2012) calibrated AUC"),
      metric_rows("cauc_pds_spatial_cv", loco$cauc_pds[i], paste0("center_", loco$centre[i]), unit = "", note = "Hijmans (2012): AUC after pairwise distance sampling, calibrated"),
      metric_rows("boyce_spatial_cv", loco$boyce[i], paste0("center_", loco$centre[i]), unit = "")))),
    if (nrow(loco) > 0L) metric_rows("auc_spatial_cv", min(loco$auc), "center_min", unit = "", note = "Minimum AUC across held-out centres"),
    if (nrow(loco) > 0L) metric_rows("cauc_spatial_cv", min(loco$cauc), "center_min", unit = "", note = "Minimum calibrated AUC across held-out centres"),
    if (nrow(loco) > 0L) metric_rows("cauc_pds_spatial_cv", min(loco$cauc_pds), "center_min", unit = "", note = "Minimum PDS-calibrated AUC across held-out centres")
  )
  # un centro fuera con el modelo de splines (modelo principal de validación)
  loco_spl <- tibble::tibble(); m_loco_spl <- NULL
  if (ns_df > 0L) {
    loco_spl <- purrr::map_dfr(centres, function(cn) {
      te <- cen_lab == cn
      tr <- !te & !(A_main$blk %in% unique(A_main$blk[te]))
      r <- run_split(A_main, tr, te, fsets[final_fs], ns_df, bases, seed_i = seed)
      if (is.null(r)) NULL else dplyr::mutate(r, centre = cn, n_presences = sum(te & A_main$yy == 1L))
    })
    if (nrow(loco_spl) > 0L) m_loco_spl <- dplyr::bind_rows(
      dplyr::bind_rows(lapply(seq_len(nrow(loco_spl)), function(i) dplyr::bind_rows(
        metric_rows("auc_spatial_cv", loco_spl$auc[i], paste0("spl_center_", loco_spl$centre[i]), note = "Spline model; trained on other centres' blocks, tested on this centre"),
        metric_rows("cauc_pds_spatial_cv", loco_spl$cauc_pds[i], paste0("spl_center_", loco_spl$centre[i]), note = "Spline model; AUC after pairwise distance sampling (Hijmans 2012), calibrated"),
        metric_rows("boyce_spatial_cv", loco_spl$boyce[i], paste0("spl_center_", loco_spl$centre[i]))))),
      metric_rows("auc_spatial_cv", min(loco_spl$auc), "spl_center_min", note = "Spline model; minimum AUC across held-out centres"),
      metric_rows("cauc_pds_spatial_cv", min(loco_spl$cauc_pds), "spl_center_min", note = "Spline model; minimum PDS-calibrated AUC across held-out centres"))
  }
  # efectos por centro (modelo propio en cada centro)
  centre_fit <- purrr::map_dfr(centres, function(cn) {
    ip <- which(main$P$centre_census == cn); id <- which(main$D$centre == cn)
    if (length(ip) < min_n || length(id) < 200L) return(NULL)
    f <- fit_set(main$P[ip, ], main$D[id, ], vars_fin, cc, zs)
    tibble::tibble(centre = cn, n_presences = length(ip), variable = vars_fin, intensity_ratio = exp(f$beta[vars_fin]),
                   ci_low = exp(f$beta[vars_fin] - 1.96 * f$se[vars_fin]), ci_high = exp(f$beta[vars_fin] + 1.96 * f$se[vars_fin]))
  })
  m_centre_ir <- if (nrow(centre_fit) > 0L) dplyr::bind_rows(lapply(seq_len(nrow(centre_fit)), function(i)
    metric_rows(paste0("intensity_ratio_", M03_VAR_LABEL[[centre_fit$variable[i]]]), centre_fit$intensity_ratio[i], paste0("center_", centre_fit$centre[i]),
                centre_fit$ci_low[i], centre_fit$ci_high[i], "ratio per SD", paste0("Model fitted within the centre (n = ", centre_fit$n_presences[i], "); ", wald_note)))) else NULL

  brk <- as.numeric(unlist(p$rings_km %||% c(0, 2, 5, 10))) * 1000
  ring_of <- function(d) cut(d, c(brk, Inf), right = FALSE, labels = paste0(brk / 1000, "_", c(brk[-1] / 1000, "inf"), "km"))
  ring_tab <- tibble::tibble(ring = levels(ring_of(0)), n_presences = as.integer(table(ring_of(main$P$d_nucleo))[levels(ring_of(0))]),
                             area_km2 = as.numeric(table(ring_of(main$D$d_nucleo))[levels(ring_of(0))]) * cc) |>
    dplyr::mutate(density_per_km2 = n_presences / area_km2) |> dplyr::filter(n_presences > 0L, area_km2 > 0)
  ring_fit <- purrr::map_dfr(ring_tab$ring[ring_tab$n_presences >= min_n], function(rg) {
    ip <- which(ring_of(main$P$d_nucleo) == rg); id <- which(ring_of(main$D$d_nucleo) == rg)
    v2 <- setdiff(vars_fin, "log_d_nucleo")
    zr <- zs
    f <- fit_set(main$P[ip, ], main$D[id, ], v2, cc, zr)
    tibble::tibble(ring = rg, variable = v2, intensity_ratio = exp(f$beta[v2]), ci_low = exp(f$beta[v2] - 1.96 * f$se[v2]), ci_high = exp(f$beta[v2] + 1.96 * f$se[v2]))
  })
  m_ring <- dplyr::bind_rows(
    dplyr::bind_rows(lapply(seq_len(nrow(ring_tab)), function(i) metric_rows("building_density_per_km2", ring_tab$density_per_km2[i], paste0("ring_", ring_tab$ring[i]), unit = "buildings/km2", note = paste0("n = ", ring_tab$n_presences[i])))),
    if (nrow(ring_fit) > 0L) dplyr::bind_rows(lapply(seq_len(nrow(ring_fit)), function(i) metric_rows(paste0("intensity_ratio_", M03_VAR_LABEL[[ring_fit$variable[i]]]), ring_fit$intensity_ratio[i], paste0("ring_", ring_fit$ring[i]),
                                                                                                    ring_fit$ci_low[i], ring_fit$ci_high[i], "ratio per SD", wald_note))))

  # ----------------------------------------------------------
  # 8. Sensibilidades
  # ----------------------------------------------------------
  say("Sensibilidades...")
  m_sens <- list(); sens_tab <- list()
  # fondo
  for (r in ladder) {
    f <- conv_tab[conv_tab$ratio == r, ]
    sub <- paste0("sens_bgratio_", r / scale_test)
    m_sens[[length(m_sens) + 1L]] <- dplyr::bind_rows(lapply(seq_len(nrow(f)), function(i) metric_rows(paste0("intensity_ratio_", M03_VAR_LABEL[[f$variable[i]]]), f$intensity_ratio[i], sub,
                                                                                                    exp(f$beta[i] - 1.96 * f$se[i]), exp(f$beta[i] + 1.96 * f$se[i]), "ratio per SD", wald_note)))
    m_sens[[length(m_sens) + 1L]] <- metric_rows("log_likelihood", f$log_likelihood[1], sub, note = "IPP log-likelihood with this background size")
  }
  # logit ingenuo
  nb <- stats::coef(naive)[paste0("z_", vars_fin)]
  m_sens[[length(m_sens) + 1L]] <- dplyr::bind_rows(lapply(seq_along(vars_fin), function(i) metric_rows(paste0("odds_ratio_", M03_VAR_LABEL[[vars_fin[i]]]), exp(nb[[i]]), "sens_logit_3to1",
                                                                                                         exp(naive_ci[paste0("z_", vars_fin[i]), 1]), exp(naive_ci[paste0("z_", vars_fin[i]), 2]), "odds ratio per SD",
                                                                                                         paste0("Unweighted logistic regression, background:presence = ", r_naive / scale_test, ":1; biased relative to the IPP when the model is misspecified (Fithian & Hastie 2013)"))))
  # semillas de fondo
  for (s in names(cache$seed_pools)) {
    Ds <- cache$seed_pools[[s]]; Ds <- Ds[is.finite(Ds$slope), ]; Ds <- Ds[seq_len(min(nrow(Ds), nD)), ]
    f <- fit_set(main$P, Ds, vars_fin, dom_area[[dom_main]] / nrow(Ds), zs)
    m_sens[[length(m_sens) + 1L]] <- coef_rows(f, vars_fin, paste0("sens_seed_", s), note = wald_note)
  }
  # dominio
  dom_summary <- list()
  for (dm in dom_sens) {
    d <- dom_data(dm, n_back = NULL); ccd <- dom_area[[dm]] / nrow(d$D)
    if (nrow(d$P) < min_n) next
    zsd <- zs_of(d$D)
    f <- fit_set(d$P, d$D, vars_fin, ccd, zsd)
    sub <- paste0("sens_domain_", dm)
    m_sens[[length(m_sens) + 1L]] <- coef_rows(f, vars_fin, sub, note = wald_note)
    zsave <- zs; Ad <- add_z(dplyr::bind_rows(d$P |> dplyr::mutate(yy = 1L), d$D |> dplyr::mutate(yy = 0L)), zsd)
    cv_d <- {
      fo <- spatial_block_folds(Ad$x, Ad$y, block_m, k_cv, seed + 1L)
      purrr::map_dfr(seq_len(k_cv), function(fk) {
        tr <- fo != fk; te <- fo == fk; trp <- tr & Ad$yy == 1L; trd <- tr & Ad$yy == 0L; tep <- te & Ad$yy == 1L; ted <- te & Ad$yy == 0L
        if (sum(tep) < 5L) return(NULL)
        g_tr <- sf::st_as_sf(Ad[trp, c("x", "y")], coords = c("x", "y"), crs = crs)
        nr <- function(ix) { g <- sf::st_as_sf(Ad[ix, c("x", "y")], coords = c("x", "y"), crs = crs); as.numeric(sf::st_distance(g, g_tr[sf::st_nearest_feature(g, g_tr), ], by_element = TRUE)) }
        dpp <- nr(tep); ddd <- nr(ted); nauc <- auc_rank(-dpp, -ddd)
        fx <- ipp_fit(mk_X(Ad[trp, ], vars_fin), mk_X(Ad[trd, ], vars_fin), ccd)
        ep <- drop(mk_X(Ad[tep, ], vars_fin) %*% fx$beta); ed <- drop(mk_X(Ad[ted, ], vars_fin) %*% fx$beta)
        a <- auc_rank(ep, ed); tibble::tibble(auc = a, cauc = a + 0.5 - max(0.5, nauc), boyce = boyce_cont(ep, ed)$index)
      })
    }
    m_sens[[length(m_sens) + 1L]] <- dplyr::bind_rows(
      metric_rows("auc_spatial_cv", mean(cv_d$auc), sub, min(cv_d$auc), max(cv_d$auc), note = "Mean over spatial block folds (1 repeat). AUC and Boyce depend on the extent of the domain and are NOT comparable across domains (Lobo et al. 2008; Hirzel et al. 2006)"),
      metric_rows("cauc_spatial_cv", mean(cv_d$cauc), sub, min(cv_d$cauc), max(cv_d$cauc), note = "Hijmans (2012) calibrated AUC"),
      metric_rows("boyce_spatial_cv", mean(cv_d$boyce, na.rm = TRUE), sub),
      metric_rows("n_presences", nrow(d$P), sub, unit = "count"), metric_rows("n_background", nrow(d$D), sub, unit = "count"))
    dom_summary[[dm]] <- tibble::tibble(domain = dm, area_km2 = dom_area[[dm]], n_presences = nrow(d$P), n_background = nrow(d$D))
  }
  # tamaño de bloque
  bs_list <- as.numeric(unlist(p$cv$block_sizes_sensitivity %||% numeric()))
  bs_tab <- purrr::map_dfr(bs_list, function(S) {
    r <- tryCatch(cv_run(A_main, fsets[final_fs], S, 1L), error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0L) return(NULL)
    r |> dplyr::summarise(block_size_m = S, auc_lo = min(auc), auc_hi = max(auc), auc = mean(auc), cauc = mean(cauc), boyce = mean(boyce, na.rm = TRUE))
  })
  if (nrow(bs_tab) > 0L) for (i in seq_len(nrow(bs_tab))) m_sens[[length(m_sens) + 1L]] <- dplyr::bind_rows(
    metric_rows("auc_spatial_cv", bs_tab$auc[i], paste0("sens_blocksize_", bs_tab$block_size_m[i]), bs_tab$auc_lo[i], bs_tab$auc_hi[i], note = "Mean over folds (1 repeat)"),
    metric_rows("cauc_spatial_cv", bs_tab$cauc[i], paste0("sens_blocksize_", bs_tab$block_size_m[i])),
    metric_rows("boyce_spatial_cv", bs_tab$boyce[i], paste0("sens_blocksize_", bs_tab$block_size_m[i])))
  # sensibilidad a los grados de libertad de los splines (el modelo principal usa nonlinear_df)
  if (ns_df > 0L) for (dfk in setdiff(as.integer(unlist(p$nonlinear_df_sensitivity %||% c(2L, 5L))), ns_df)) {
    bases_k <- lapply(stats::setNames(vars_fin, vars_fin), function(v) splines::ns(Dz[[paste0("z_", v)]], df = dfk))
    cv_k <- cv_run(A_main, fsets[final_fs], block_m, 1L, ns_df = dfk, bases = bases_k)
    m_sens[[length(m_sens) + 1L]] <- dplyr::bind_rows(
      metric_rows("cauc_pds_spatial_cv", mean(cv_k$cauc_pds), paste0("sens_spline_df_", dfk), min(cv_k$cauc_pds), max(cv_k$cauc_pds), note = "Mean over folds (1 repeat)"),
      metric_rows("cv_log_likelihood", mean(cv_k$cv_ll), paste0("sens_spline_df_", dfk), note = "Mean held-out IPP log-likelihood per fold"))
  }
  m_sens <- dplyr::bind_rows(m_sens)

  # ----------------------------------------------------------
  # 9. Participaciones descriptivas (umbrales de planificación)
  # ----------------------------------------------------------
  say("Participaciones descriptivas...")
  share_rows <- list(); share_tab <- list()
  for (thr in thr_road) for (src in c("igm_vias", "igm_all", "osm2019", "census2020")) {
    dcol <- paste0("d_", src)
    k <- sum(main$P[[dcol]] <= thr); ci <- wilson_ci(k, nP); e <- mean(main$D[[dcol]] <= thr)
    share_rows[[length(share_rows) + 1L]] <- metric_rows(paste0("share_within_", format(thr, scientific = FALSE), "m_road"), k / nP, src, ci[1], ci[2], "proportion",
                                                          paste0("Share of dispersed buildings within ", thr, " m of the road source; background (area) share = ", round(e, 3), "; 95% Wilson interval"))
    share_tab[[length(share_tab) + 1L]] <- tibble::tibble(road_source = src, threshold_m = thr, share_buildings = k / nP, ci_low = ci[1], ci_high = ci[2], share_area = e, enrichment = (k / nP) / e)
  }
  # umbrales normativos de sensibilidad (25 y 30 m) solo sobre proxies de la red vial estatal
  for (thr in thr_road_sens) for (nm in names(state_nets)) {
    dcol <- paste0("d_state_", nm)
    k <- sum(main$P[[dcol]] <= thr); ci <- wilson_ci(k, nP); e <- mean(main$D[[dcol]] <= thr)
    share_rows[[length(share_rows) + 1L]] <- metric_rows(paste0("share_within_", format(thr, scientific = FALSE), "m_state_road"), k / nP, nm, ci[1], ci[2], "proportion",
                                                          paste0("Share of dispersed buildings within ", thr, " m of the state-network proxy ", nm, "; background (area) share = ", round(e, 3), "; 95% Wilson interval"))
    share_tab[[length(share_tab) + 1L]] <- tibble::tibble(road_source = paste0("state:", nm), threshold_m = thr, share_buildings = k / nP, ci_low = ci[1], ci_high = ci[2], share_area = e, enrichment = (k / nP) / e)
  }
  slope_ref <- tibble::tibble(slope_deg = c(10, 15, 20, 25, 30, 40)) |>
    dplyr::mutate(share_buildings_above = vapply(slope_deg, function(s) mean(main$P$slope > s), numeric(1)), share_area_above = vapply(slope_deg, function(s) mean(main$D$slope > s), numeric(1)))
  if (!is.null(slope_thr) && is.numeric(slope_thr)) {
    k <- sum(main$P$slope > slope_thr); ci <- wilson_ci(k, nP)
    share_rows[[length(share_rows) + 1L]] <- metric_rows("share_above_slope_threshold", k / nP, final_fs, ci[1], ci[2], "proportion", paste0("Share of buildings on slopes > ", slope_thr, " degrees (threshold from the protocol)"))
  }
  # límite urbano oficial 2021 (cubre Loja, Malacatos y Vilcabamba)
  ul_tab <- NULL
  if (nzchar(clean_chr(dat$urban_limit_vector %||% ""))) {
    ul <- sf::st_union(rd(dat$urban_limit_vector))
    insideU <- lengths(sf::st_intersects(sf::st_as_sf(main$P, coords = c("x", "y"), crs = crs), ul)) > 0
    ul_tab <- main$P |> dplyr::mutate(inside = insideU) |> dplyr::filter(centre_census %in% c("Loja", "Malacatos", "Vilcabamba")) |>
      dplyr::group_by(centre = centre_census) |> dplyr::summarise(n = dplyr::n(), n_inside = sum(inside), share_inside = n_inside / n, .groups = "drop")
    for (i in seq_len(nrow(ul_tab))) {
      ci <- wilson_ci(ul_tab$n_inside[i], ul_tab$n[i])
      share_rows[[length(share_rows) + 1L]] <- metric_rows("share_inside_urban_limit", ul_tab$share_inside[i], paste0("center_", ul_tab$centre[i]), ci[1], ci[2], "proportion",
                                                           "Dispersed buildings (outside census blocks) inside the official 2021 urban limit; Catamayo and San Pedro de la Bendita have no limit in the data")
    }
  }
  m_share <- dplyr::bind_rows(share_rows)

  # ----------------------------------------------------------
  # 10. Métricas, tablas y figuras
  # ----------------------------------------------------------
  metrics_all <- dplyr::bind_rows(m_ir, m_spl_eff, m_cv, m_rnd, m_extra, m_loco, m_loco_spl, m_centre_ir, m_ring, m_sens, m_share) |>
    dplyr::mutate(value = as.numeric(value)) |> dplyr::distinct(metric, subset, .keep_all = TRUE)
  # Main metrics (hypotheses, effects, validation of the main models, centres, rings, shares) go to
  # `metrics`; las demás (sensibilidades, duplicados por conjunto de variables, métricas de apoyo) van íntegras a tab_M03_metrics_supplementary,
  # con las mismas columnas. nrow(metrics) - filas añadidas + nrow(supplementary) = nrow(metrics_all) (se comprueba abajo).
  all_sets <- c(names(fsets), paste0("spl_", names(fsets)))
  report_sets <- c(final_fs, paste0("spl_", final_fs), "osm2019", "spl_osm2019")      # feature sets whose validation is reported (H4, H5)
  cv_reported <- c("auc_spatial_cv", "auc_pds_spatial_cv", "cauc_pds_spatial_cv", "cauc_spatial_cv", "boyce_spatial_cv", "cv_log_likelihood", "auc_random_cv")
  cv_design <- c("auc_null_distance", "spatial_sorting_bias", "spatial_sorting_bias_after_pds")   # idénticas en todos los conjuntos: una sola fila
  to_supp <- startsWith(metrics_all$subset, "sens_") |
    (metrics_all$subset %in% all_sets & metrics_all$metric %in% cv_reported & !metrics_all$subset %in% report_sets) |
    (metrics_all$subset %in% all_sets & metrics_all$metric %in% cv_design & metrics_all$subset != final_fs) |
    (metrics_all$subset %in% all_sets & metrics_all$metric %in% c("log_likelihood", "aic", "n_presences", "n_background")) |
    (grepl("^(spl_)?center_", metrics_all$subset) & metrics_all$metric %in% c("boyce_spatial_cv", "cauc_spatial_cv"))
  metrics_supp <- metrics_all[to_supp, ]
  metrics_main_rows <- dplyr::bind_rows(metric_rows("n_presences", nP, "main_domain", unit = "count", note = paste0("Dispersed buildings in the main domain '", dom_main, "'")),
                                        metric_rows("n_background", nD, "main_domain", unit = "count", note = "Uniform background points of the main model"))
  metrics <- dplyr::bind_rows(metrics_all[!to_supp, ], metrics_main_rows)
  stopifnot(nrow(metrics) - nrow(metrics_main_rows) + nrow(metrics_supp) == nrow(metrics_all))
  say("Métricas: ", nrow(metrics), " en `metrics` y ", nrow(metrics_supp), " en la tabla suplementaria (total ", nrow(metrics_all), ")")

  cv_fold_tab <- dplyr::bind_rows(dplyr::mutate(cv_sp, scheme = "spatial_block"), dplyr::mutate(cv_rnd, scheme = "random"))
  cv_tab <- dplyr::bind_rows(dplyr::mutate(s_sp, scheme = "spatial_block"), dplyr::mutate(s_rnd, scheme = "random"))
  tables <- list(
    tab_M03_coefficients = list(data = coef_tab, description = "Standardised intensity ratios (exp of coefficient per SD) with model-based SE and block-bootstrap percentile intervals, by feature set"),
    tab_M03_cv_summary = list(data = cv_tab, description = "Cross-validation summary by feature set and scheme (mean and min-max over folds): AUC, pairwise-distance AUC, distance-null AUC, calibrated AUC, Boyce"),
    tab_M03_cv_folds = list(data = cv_fold_tab, description = "Fold-level cross-validation results"),
    tab_M03_loco = list(data = loco, description = "Leave-one-centre-out validation of the final feature set (linear model)"),
    tab_M03_loco_spline = list(data = loco_spl, description = "Leave-one-centre-out validation of the final feature set (spline model)"),
    tab_M03_response_curves = list(data = resp_tab, description = "Spline model: intensity relative to the median predictor value, with block-bootstrap 95% band, in the original predictor units"),
    tab_M03_model_comparison = list(data = {
      lin <- purrr::imap_dfr(fits, function(f, nm) tibble::tibble(model = "linear", feature_set = nm, n_parameters = length(f$beta), log_likelihood = f$ll, aic = -2 * f$ll + 2 * length(f$beta)))
      spl <- purrr::imap_dfr(fits_spl, function(f, nm) tibble::tibble(model = paste0("splines_df", ns_df), feature_set = nm, n_parameters = length(f$beta), log_likelihood = f$ll, aic = -2 * f$ll + 2 * length(f$beta)))
      cvs <- cv_tab |> dplyr::filter(scheme == "spatial_block") |> dplyr::transmute(key = feature_set, cv_loglik_mean = cv_ll__mean, auc_pds = auc_pds__mean, cauc_pds = cauc_pds__mean, boyce = boyce__mean)
      dplyr::bind_rows(lin, spl) |> dplyr::mutate(key = ifelse(model == "linear", feature_set, paste0("spl_", feature_set))) |> dplyr::left_join(cvs, by = "key") |> dplyr::select(-key)
    }, description = "Linear versus spline model: log-likelihood, AIC and spatially cross-validated held-out log-likelihood, AUC after pairwise distance sampling and Boyce"),
    tab_M03_background_convergence = list(data = conv_tab, description = "Intensity ratios and log-likelihood as the background size increases (Renner et al. 2015 convergence check)"),
    tab_M03_loglik_by_background = list(data = ll_by_ratio, description = "Change in log-likelihood between successive background sizes"),
    tab_M03_naive_vs_ipp = list(data = tibble::tibble(variable = M03_VAR_LABEL[vars_fin], odds_ratio_logit_naive = exp(nb), intensity_ratio_ipp = exp(f_fin$beta[vars_fin])),
                                description = "Unweighted logistic regression (3:1) versus IPP estimates for the final feature set"),
    tab_M03_collinearity_vif = list(data = coll, description = "VIF and maximum absolute correlation within each feature set (Dormann et al. 2013)"),
    tab_M03_collinearity_pairs = list(data = coll_pairs, description = "Pearson correlations among all predictors on the background sample"),
    tab_M03_residual_correlogram = list(data = corr, description = "Moran's I of 1-km Pearson residuals by distance band"),
    tab_M03_domains = list(data = dplyr::bind_rows(tibble::tibble(domain = dom_main, area_km2 = dom_area[[dom_main]], n_presences = nP, n_background = nD), dplyr::bind_rows(dom_summary)), description = "Study domains used in the main model and in the domain sensitivity"),
    tab_M03_block_size_sensitivity = list(data = bs_tab, description = "Spatial CV results for alternative block sizes"),
    tab_M03_rings = list(data = ring_tab, description = "Building counts, area and density by distance ring from the urban core"),
    tab_M03_ring_effects = list(data = ring_fit, description = "Standardised intensity ratios within distance rings"),
    tab_M03_centre_effects = list(data = centre_fit, description = "Standardised intensity ratios fitted separately within each centre"),
    tab_M03_shares = list(data = dplyr::bind_rows(share_tab), description = "Share of dispersed buildings within each road-distance threshold (descriptive 250 m; regulatory sensitivity 25 and 30 m) versus share of area"),
    tab_M03_slope_reference = list(data = slope_ref, description = "Descriptive share of buildings and area above reference slope values (not official planning thresholds)"),
    tab_M03_urban_limit = list(data = ul_tab %||% tibble::tibble(), description = "Dispersed buildings inside the official 2021 urban limit by centre")
  )

  tables <- Filter(function(t) is.data.frame(t$data) && nrow(t$data) > 0L && ncol(t$data) > 0L, tables)
  # Figuras que salen de tablas: se dibujan con methods/M03_figuras.R (el mismo código que usa la etapa 08 para redibujarlas sin repetir el modelo).
  figs_env <- new.env(parent = environment()); sys.source(file.path(ctx$root, "methods", "M03_figuras.R"), envir = figs_env)
  fig_eff <- figs_env$m03_fig_effects(coef_tab)
  fig_cv <- figs_env$m03_fig_validation(cv_fold_tab)
  fig_conv <- figs_env$m03_fig_convergence(conv_tab)
  fig_cor <- figs_env$m03_fig_correlogram(corr, block_m / 1000)
  pe <- pe_curve |> dplyr::group_by(mid) |> dplyr::summarise(lo = min(pe, na.rm = TRUE), hi = max(pe, na.rm = TRUE), pe = stats::median(pe, na.rm = TRUE), .groups = "drop")
  fig_pe <- ggplot2::ggplot(pe, ggplot2::aes(mid, pe)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) + ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    ggplot2::labs(x = "Predicted intensity (percentile of background)", y = "Predicted-to-expected ratio") + theme_paper()
  figures <- list(
    list(id = "fig_M03_effects", plot = fig_eff, width_cm = 13, height_cm = 8, caption_draft = "Standardised intensity ratios (per SD of each predictor) for each feature set with block-bootstrap 95% intervals."),
    list(id = "fig_M03_validation", plot = fig_cv, width_cm = 17, height_cm = 7, caption_draft = "Cross-validated AUC, distance-null-calibrated AUC and continuous Boyce index under random and spatial block k-fold schemes."),
    list(id = "fig_M03_background_convergence", plot = fig_conv, width_cm = 12, height_cm = 8, caption_draft = "Intensity ratios of the final model as the number of background points per presence increases."),
    list(id = "fig_M03_correlogram", plot = fig_cor, width_cm = 12, height_cm = 8, caption_draft = "Moran's I of 1-km Pearson residuals of the final model by distance band; the dashed line marks the block size used for validation."),
    list(id = "fig_M03_boyce_curve", plot = fig_pe, width_cm = 11, height_cm = 8, caption_draft = "Predicted-to-expected frequency of held-out dispersed buildings along the predicted-intensity gradient (median and range across spatial folds).")
  )
  if (nrow(loco_spl) > 0L) {   # modelo final (spline) y métrica del criterio de H6 (AUC tras muestreo por pares)
    thr_h6 <- tryCatch({ h <- Filter(function(x) identical(as.character(x$subset), "spl_center_min"), read_frozen_protocol()$payload$hypotheses); as.numeric(h[[1]]$threshold) }, error = function(e) NA_real_)
    figures[[length(figures) + 1L]] <- list(id = "fig_M03_centres", width_cm = 13, height_cm = 8,
      plot = figs_env$m03_fig_centres(loco_spl, threshold = thr_h6),
      caption_draft = "Leave-one-centre-out validation of the final spline model: AUC, AUC after pairwise distance sampling (cAUC_pds) and continuous Boyce index for each held-out centre.")
  }
  if (nrow(resp_tab) > 0L) {
    fig_resp <- figs_env$m03_fig_response(resp_tab)
    figures[[length(figures) + 1L]] <- list(id = "fig_M03_response_curves", plot = fig_resp, width_cm = 17, height_cm = 7,
      caption_draft = "Response curves of the spline model: intensity of dispersed buildings relative to the median value of each predictor, with block-bootstrap 95% bands.")
  }
  # mapa de residuos por celda de 1 km (Lobo et al. 2008: los índices de un solo número no muestran dónde falla el modelo)
  fig_res <- figs_env$m03_fig_residual_map(cells, cell_m)
  figures[[length(figures) + 1L]] <- list(id = "fig_M03_residual_map", plot = fig_res, width_cm = 14, height_cm = 12,
    caption_draft = "Pearson residuals (observed minus expected dispersed buildings, scaled) of the final model in 1-km cells; orange cells have more buildings than the model expects and blue cells fewer.")
  tables$tab_M03_metrics_supplementary <- list(data = metrics_supp,
    description = "Supplementary metrics (sensitivity analyses, per-feature-set duplicates, supporting validation metrics); same columns as metrics_long.csv, nothing omitted")
  # modelling domains for the study-area map (layer "main" = main domain)
  dom_file <- file.path(ctx$DIR$processed_vectors, "model_domain.gpkg")
  ensure_dir(dirname(dom_file)); if (file.exists(dom_file)) unlink(dom_file)
  # contornos sin huecos (las manzanas y el parque se restan del dominio de modelación, pero el mapa solo necesita el contorno)
  outline_of <- function(g, fill = TRUE) tryCatch({
    pl <- sf::st_collection_extract(sf::st_make_valid(g), "POLYGON")
    pl <- sf::st_cast(sf::st_sfc(sf::st_union(pl), crs = sf::st_crs(pl)), "POLYGON")
    if (fill) pl <- sf::st_sfc(lapply(pl, function(q) sf::st_polygon(list(q[[1]]))), crs = sf::st_crs(pl))
    sf::st_sf(geometry = pl)
  }, error = function(e) NULL)
  for (nm in c("main", "aoi_amplio", "urban_cores")) {
    pg <- outline_of(switch(nm, main = dom_def[[dom_main]], aoi_amplio = dom_def$aoi_amplio, urban_cores = cores), fill = nm != "urban_cores")
    if (!is.null(pg)) { pg$layer_name <- nm; suppressWarnings(sf::st_write(pg, dom_file, layer = nm, quiet = TRUE, append = FALSE)) }
  }
  tables$tab_M03_residual_cells <- list(data = cells, description = "Observed and expected dispersed buildings and Pearson residuals per 1-km cell (final model)")

  # El mapa del dominio y los edificios dispersos por centro (fig_M03_domain) lo dibuja la etapa 08 con methods/M03_figuras.R a partir de los vectores procesados.

  # ----------------------------------------------------------
  # 11. Hechos y limitaciones
  # ----------------------------------------------------------
  runtime <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1)
  facts <- list(
    algorithm = "Inhomogeneous Poisson point process (log-linear intensity) fitted by maximising the point-process likelihood with uniform background (quadrature) points; Newton-Raphson, implemented in base R",
    equivalence = "Equivalent in the limit to infinitely weighted logistic regression (Fithian & Hastie 2013) and to a weighted Poisson regression (Warton & Shepherd 2010; Renner et al. 2015)",
    presences = paste0(nP, " census buildings (2021) located outside census blocks and outside Podocarpus National Park, within the domain '", dom_main, "'"),
    domain_main = dom_main, domain_area_km2 = round(dom_area[[dom_main]], 1),
    background = paste0(nD, " points sampled uniformly at random in the domain (", ratio_main / scale_test, " per presence; seed ", seeds[1], "); area per point = ", signif(cc, 4), " km2"),
    background_convergence = paste0("ladder (points per presence): ", paste(ladder / scale_test, collapse = ", "), "; convergence criterion: relative change in each intensity ratio < 2% and change in log-likelihood within 2 Monte Carlo standard errors of the difference between the two largest sizes (random quadrature points make the log-likelihood fluctuate; Renner et al. 2015 use a change < 2 with a regular mesh) ->", if (convergence_ok) "met" else "NOT met"),
    predictors = paste(M03_VAR_LABEL[allvars], collapse = ", "),
    predictor_transform = "distances as log(1 + d in metres); slope in degrees (SRTM 30 m); all predictors standardised with the mean and SD of the background sample, so ratios are per SD",
    effects_interval = paste0("block bootstrap, B = ", B, ", blocks of ", block_m, " m, percentile 95% interval"),
    residual_autocorrelation = paste0("Moran's I of 1-km Pearson residuals (99 permutations): significant positive autocorrelation up to ", range_m, " m"),
    cross_validation = paste0(k_cv, "-fold spatial block CV, blocks of ", block_m, " m (rule: ", clean_chr(p$cv$block_size_m), "; bounds ", paste(bnd, collapse = "-"), " m), ", n_rep, " repeats; random k-fold shown for reference"),
    validation_metrics = "AUC (Mann-Whitney), AUC after pairwise distance sampling, distance-null AUC, calibrated AUC = AUC + 0.5 - max(0.5, null AUC) (Hijmans 2012), continuous Boyce index on the rank scale of the held-out background (Hirzel et al. 2006)",
    leave_one_centre_out = paste(centres, collapse = ", "),
    sensitivity = paste0("background size, naive logistic 3:1, background seeds (", paste(seeds, collapse = ", "), "), domain (", paste(dom_sens, collapse = ", "), "), block sizes (", paste(bs_list, collapse = ", "), " m), natural splines (df = ", ns_df, ")"),
    collinearity = paste0("max |r| among all predictors = ", round(max(abs(coll_pairs$pearson_r)), 2), "; max VIF within feature sets = ", round(max(coll$vif), 2)),
    shares = paste0("share of buildings within ", thr_road, " m of each road source versus share of area (descriptive); share within ", paste(thr_road_sens, collapse = " and "), " m only for proxies of the state road network (", paste(names(state_nets), collapse = ", "), "), as sensitivity thresholds from the Ley de Caminos regulation Art. 4; slope shares are descriptive (no official threshold)"),
    random_seed = seed,
    metrics_policy = paste0(nrow(metrics), " metrics are reported in metrics_long.csv (hypotheses, effects, validation of the main models, centres, rings, shares); the other ", nrow(metrics_supp), " (sensitivity analyses and supporting metrics) are in tab_M03_metrics_supplementary"),
    model_roles = if (ns_df > 0L) paste0("Primary model for fit, cross-validation and response curves: natural cubic splines (df = ", ns_df, " per predictor), subsets 'spl_<feature set>'. Log-linear model: summary of the average effect as intensity ratios per SD, used for the directional hypotheses on slope, roads, core and river (subsets '<feature set>')") else "Log-linear model only",
    model_comparison = "Linear versus spline model compared by held-out (spatial block) IPP log-likelihood, AUC after pairwise distance sampling and Boyce; in-sample AIC is reported but overstates evidence under clustering",
    co_principal_domain = "aoi_amplio is reported as a sensitivity domain (subset sens_domain_aoi_amplio) next to the main domain",
    odmap = "Structure follows ODMAP (Zurell et al. 2020): see methods/M03_presencia_fondo_README.md"
  )
  limitations <- c(
    "Presences are census buildings, not dwellings or residences; seasonal use cannot be separated at building level.",
    "The model describes where dispersed buildings are located relative to the available environment (relative intensity), not the probability of conversion or the processes that created them; effects are associations.",
    "The study domain (AOI) was built with slope, road-corridor and map-house rules, so part of the association may be built into the domain; sensitivity to the AOI scenarios and to a wider domain is reported but a domain fully independent of the covariates is not available.",
    "Roads from the 1:50,000 base (2009-2013), OSM 2019 and census axes 2020-2021 are not contemporaneous with 2021 buildings; recent roads and tracks are partly endogenous to settlement (Chomitz & Gray 1996; Duranton & Turner 2012), so intensity ratios of road distance are not causal effects.",
    "Intensity ratios per SD from the log-linear model summarise an average effect of relationships that are not linear in the log-intensity; response curves and ratios between the 25th and 75th percentile come from the spline model.",
    "Buildings are spatially clustered; the Poisson model assumes conditional independence, so model-based standard errors are too small; intervals reported for the main effects use a block bootstrap.",
    "Presence-background AUC is inflated by spatial sorting bias and depends on the background extent (Lobo et al. 2008; Hijmans 2012); AUC after pairwise distance sampling, calibrated AUC, the Boyce index and the effects are reported alongside it, and AUC or Boyce values are not compared across domains of different extent.",
    "The Boyce index is computed on the rank scale of the held-out background (equal-area windows of width 0.1), a variant of the continuous index of Hirzel et al. (2006) that avoids empty windows when intensities have a long tail; it is rank-based and cannot distinguish the shape of the predicted-to-expected curve.",
    "The census cartography may map dispersed buildings less completely away from roads (survey-effort bias, Phillips et al. 2009); in that case the association with road distance is overstated. A census of buildings has no opportunistic sampling, but the completeness of the 2021 update by distance to roads is not documented.",
    "The Boyce-based and AUC-based metrics describe discrimination of building locations within the domain and say nothing about the accuracy of predicted densities; the spatial distribution of errors is shown in the residual map (Lobo et al. 2008).",
    "No independent validation data (a different building inventory) are available; validation is internal (spatial block and leave-one-centre-out).",
    "Official planning thresholds for slope and road distance are not available; reported shares are descriptive.",
    "The state road network is proxied by IGM route class (primary; primary plus secondary) and OSM highway class; the official MTOP inventory is not available, the IGM classes date from 2013 and census axes have no hierarchy, so the 25 and 30 m shares are sensitivity values, not compliance.",
    if (!convergence_ok) "The background-size convergence criterion was not met at the largest size evaluated; coefficients may still change with a larger background."
  )
  limitations <- limitations[!vapply(limitations, is.null, logical(1))]

  list(metrics = metrics, tables = tables, figures = figures, facts = facts, limitations = limitations, products = dom_file)
}
