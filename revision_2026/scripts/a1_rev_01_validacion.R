# a1_rev_01_validacion.R
# Objetivo : analisis de validacion pedidos por el revisor 1 de Land (correcciones mayores, nov-2026), SIN tocar el protocolo congelado ni los resultados enviados:
#            V0 = validacion original (estandarizacion y nudos de splines con todo el fondo) para comprobar que se reproduce el manuscrito;
#            V1 = preprocesamiento (media, DE y nudos de splines) calculado solo con el fondo de entrenamiento de cada particion;
#            B5 / B10 = V1 mas una zona de exclusion (buffer) de 5 o 10 km entre el conjunto de prueba y el de entrenamiento.
#            Se aplican a la validacion espacial por bloques (5 pliegues x 3 repeticiones) y a la de dejar un centro fuera.
# Entradas : data/processed/M03_dataset_ea8ce5b426.rds (cache de covariables de la corrida completa del modulo M03: 17 088 presencias, fondo de 704 026 puntos)
#            parametros del protocolo (semilla 2026, bloque 10 km, 5 pliegues, 3 repeticiones, splines con 3 g.l., 20 puntos de fondo por presencia)
# Salidas  : outputs/revision/rev_cv_splits.csv (todas las particiones), rev_cv_summary.csv, rev_loco.csv, rev_fold_blocks.csv, rev_fold_separation.csv, rev_reproduccion.txt
# Uso      : Rscript a1_rev_01_validacion.R [quick]   ('quick' = una repeticion y dos pliegues/centros, solo para probar)
# Fecha    : 2026-11-07
suppressMessages({library(sf); library(terra); library(dplyr); library(purrr); library(tibble); library(splines)})
args <- commandArgs(trailingOnly = TRUE); QUICK <- "quick" %in% args
ROOT <- "C:/Trabajos_y_Publicaciones/02_cambio_ciudades/000_A1_research"
OUT <- file.path(ROOT, "outputs/revision"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf("[rev01 %s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
CRS <- "EPSG:32717"
seed <- 2026L; k_cv <- 5L; n_rep <- if (QUICK) 1L else 3L; block_m <- 10000; ns_df <- 3L; ratio_main <- 20
AREA_MAIN <- 1639.0470387658202   # km2, tab_M03_domains.csv

# ---------------- funciones (copiadas sin cambios de methods/M03_presencia_fondo.R y R_helpers.R) ----------------
auc_rank <- function(pos, neg) { n1 <- as.numeric(length(pos)); n0 <- as.numeric(length(neg)); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(c(pos, neg)); (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
boyce_cont <- function(score_p, score_d, nclass = 101L, window = 0.1) {
  Fd <- stats::ecdf(score_d); p <- Fd(score_p); b <- Fd(score_d)
  mids <- seq(window / 2, 1 - window / 2, length.out = nclass)
  P <- vapply(mids, function(m) mean(p >= m - window / 2 & p <= m + window / 2), numeric(1))
  E <- vapply(mids, function(m) mean(b >= m - window / 2 & b <= m + window / 2), numeric(1))
  pe <- ifelse(E > 0, P / E, NA_real_); ok <- is.finite(pe)
  idx <- if (sum(ok) > 5L) suppressWarnings(stats::cor(mids[ok], pe[ok], method = "spearman")) else NA_real_
  list(index = idx) }
pds_idx <- function(d_p, d_d, seed = 1L, tol = 0.33) {
  n <- length(d_d); ob <- order(d_d); srt <- d_d[ob]; avail <- rep(TRUE, n)
  set.seed(seed); ip <- integer(length(d_p)); idd <- integer(length(d_p)); k <- 0L
  for (i in sample.int(length(d_p))) {
    pos <- findInterval(d_p[i], srt); lo <- pos; hi <- pos + 1L
    while (lo >= 1L && !avail[lo]) lo <- lo - 1L
    while (hi <= n && !avail[hi]) hi <- hi + 1L
    pick <- if (lo < 1L) hi else if (hi > n) lo else if ((d_p[i] - srt[lo]) <= (srt[hi] - d_p[i])) lo else hi
    if (pick < 1L || pick > n) next
    if (abs(d_p[i] - srt[pick]) > tol * max(d_p[i], srt[pick])) next
    avail[pick] <- FALSE; k <- k + 1L; ip[k] <- i; idd[k] <- ob[pick] }
  list(p = ip[seq_len(k)], d = idd[seq_len(k)]) }
ipp_fit <- function(Xp, Xd, cc, beta = NULL, maxit = 60L, tol = 1e-9) {
  k <- ncol(Xp); sp <- colSums(Xp)
  if (is.null(beta)) { beta <- stats::setNames(rep(0, k), colnames(Xp)); beta[1] <- log(nrow(Xp) / (cc * nrow(Xd))) }
  ll_fun <- function(b) sum(Xp %*% b) - cc * sum(exp(Xd %*% b))
  ll <- ll_fun(beta); conv <- FALSE
  for (it in seq_len(maxit)) {
    mu <- cc * exp(drop(Xd %*% beta)); g <- sp - drop(crossprod(Xd, mu)); H <- crossprod(Xd * mu, Xd)
    step <- tryCatch(solve(H, g), error = function(e) NULL); if (is.null(step)) break
    t <- 1; ok <- FALSE
    for (h in 1:25) { bn <- beta + t * step; lln <- ll_fun(bn); if (is.finite(lln) && lln >= ll - 1e-12) { ok <- TRUE; break }; t <- t / 2 }
    if (!ok) break
    d <- abs(lln - ll); beta <- bn; ll <- lln
    if (d < tol * (1 + abs(ll)) || max(abs(t * step)) < 1e-8) { conv <- TRUE; break } }
  list(beta = beta, ll = ll, converged = conv) }
spatial_block_folds <- function(x, y, block_size, k, seed) {
  block_id <- paste(floor(x / block_size), floor(y / block_size), sep = "_"); blocks <- unique(block_id)
  set.seed(seed); fold_of_block <- stats::setNames(sample(rep_len(seq_len(k), length(blocks))), blocks); unname(fold_of_block[block_id]) }

# ---------------- datos: misma construccion que el modelo principal ----------------
cache <- readRDS(file.path(ROOT, "data/processed/M03_dataset_ea8ce5b426.rds"))
P0 <- cache$pres; D0 <- cache$pool
P0 <- P0[is.finite(P0$slope), ]; D0 <- D0[is.finite(D0$slope), ]
P <- P0[P0$in_extended, ]; D <- D0[D0$in_extended, ]; D <- D[seq_len(ceiling(ratio_main * nrow(P))), ]
cc <- AREA_MAIN / nrow(D)
say("Presencias: ", nrow(P), "; fondo: ", nrow(D), "; c = ", signif(cc, 4), " km2/punto (esperado 17 088, 341 760, 0.0048)")
stopifnot(nrow(P) == 17088L, nrow(D) == 341760L)
VARS <- c("slope", paste0("log_d_", c("igm_vias", "osm2019")), "log_d_nucleo")
zs_of <- function(Dd) list(mean = vapply(VARS, function(v) mean(Dd[[v]]), numeric(1)), sd = vapply(VARS, function(v) stats::sd(Dd[[v]]), numeric(1)))
add_z <- function(T, zsx) { for (v in names(zsx$mean)) T[[paste0("z_", v)]] <- (T[[v]] - zsx$mean[[v]]) / zsx$sd[[v]]; T }
mk_X <- function(T, vars, ns_df = 0L, bases = NULL) {
  cols <- lapply(vars, function(v) { z <- T[[paste0("z_", v)]]
    if (ns_df > 0L) { m <- stats::predict(bases[[v]], z); colnames(m) <- paste0(v, "_ns", seq_len(ncol(m))); m } else matrix(z, ncol = 1L, dimnames = list(NULL, v)) })
  cbind(`(Intercept)` = 1, do.call(cbind, cols)) }
keep <- c("x", "y", VARS)
A <- dplyr::bind_rows(P[, keep] |> dplyr::mutate(yy = 1L, centre = P$centre_census), D[, keep] |> dplyr::mutate(yy = 0L, centre = D$centre))
A <- as.data.frame(A)
zs_g <- zs_of(D); A_g <- add_z(A, zs_g)
bases_g <- lapply(stats::setNames(VARS, VARS), function(v) splines::ns(A_g[[paste0("z_", v)]][A_g$yy == 0L], df = ns_df))
A$blk <- paste(floor(A$x / block_m), floor(A$y / block_m))

MODELS <- list(igm_vias = list(vars = c("slope", "log_d_igm_vias", "log_d_nucleo"), ns = 0L),
               spl_igm_vias = list(vars = c("slope", "log_d_igm_vias", "log_d_nucleo"), ns = ns_df),
               osm2019 = list(vars = c("slope", "log_d_osm2019", "log_d_nucleo"), ns = 0L),
               spl_osm2019 = list(vars = c("slope", "log_d_osm2019", "log_d_nucleo"), ns = ns_df))
if (QUICK) MODELS <- MODELS[c("igm_vias", "spl_igm_vias")]
VARIANTS <- list(V0 = list(prep = "global", buffer_km = 0), V1 = list(prep = "train", buffer_km = 0),
                 B5 = list(prep = "train", buffer_km = 5), B10 = list(prep = "train", buffer_km = 10))

# plantilla ráster para la zona de exclusión (distancia al conjunto de prueba, resolución 250 m)
r_tpl <- rast(ext(min(A$x) - 1000, max(A$x) + 1000, min(A$y) - 1000, max(A$y) + 1000), resolution = 250, crs = CRS)
far_from_test <- function(te, buffer_m) {
  if (buffer_m <= 0) return(rep(TRUE, nrow(A)))
  rt <- terra::rasterize(terra::vect(cbind(A$x[te], A$y[te]), crs = CRS), r_tpl, field = 1)
  dd <- terra::extract(terra::distance(rt), cbind(A$x, A$y))[, 1]
  is.na(dd) | dd > buffer_m
}

run_split <- function(tr, te, vname, seed_i) {
  vr <- VARIANTS[[vname]]
  far <- far_from_test(te, vr$buffer_km * 1000)
  tr_eff <- tr & far
  trp <- tr_eff & A$yy == 1L; trd <- tr_eff & A$yy == 0L; tep <- te & A$yy == 1L; ted <- te & A$yy == 0L
  if (sum(tep) < 5L || sum(ted) < 5L || sum(trp) < 20L) return(NULL)
  if (vr$prep == "global") { Az <- A_g; bases <- bases_g } else {
    zs_tr <- zs_of(A[trd, ]); Az <- add_z(A, zs_tr)
    bases <- lapply(stats::setNames(VARS, VARS), function(v) splines::ns(Az[[paste0("z_", v)]][trd], df = ns_df)) }
  g_tr <- sf::st_as_sf(A[trp, c("x", "y")], coords = c("x", "y"), crs = 32717)
  near <- function(idx) { g <- sf::st_as_sf(A[idx, c("x", "y")], coords = c("x", "y"), crs = 32717)
    as.numeric(sf::st_distance(g, g_tr[sf::st_nearest_feature(g, g_tr), ], by_element = TRUE)) }
  d_p <- near(tep); d_d <- near(ted)
  pr <- pds_idx(d_p, d_d, seed_i)
  nauc <- auc_rank(-d_p, -d_d); nauc_pds <- auc_rank(-d_p[pr$p], -d_d[pr$d])
  purrr::imap_dfr(MODELS, function(m, nm) {
    f <- ipp_fit(mk_X(Az[trp, ], m$vars, m$ns, bases), mk_X(Az[trd, ], m$vars, m$ns, bases), cc)
    ep <- drop(mk_X(Az[tep, ], m$vars, m$ns, bases) %*% f$beta); ed <- drop(mk_X(Az[ted, ], m$vars, m$ns, bases) %*% f$beta)
    a <- auc_rank(ep, ed); a_pds <- auc_rank(ep[pr$p], ed[pr$d])
    tibble::tibble(variant = vname, model = nm, auc = a, auc_pds = a_pds, auc_null = nauc, auc_null_pds = nauc_pds,
                   cauc = a + 0.5 - max(0.5, nauc), cauc_pds = a_pds + 0.5 - max(0.5, nauc_pds), boyce = boyce_cont(ep, ed)$index,
                   ssb = mean(d_p) / mean(d_d), ssb_pds = mean(d_p[pr$p]) / mean(d_d[pr$d]),
                   n_train_p = sum(trp), n_train_d = sum(trd), n_excl_buffer_p = sum(tr & !far & A$yy == 1L), n_excl_buffer_d = sum(tr & !far & A$yy == 0L),
                   n_test_p = sum(tep), n_test_d = sum(ted), n_pds_pairs = length(pr$p), pds_share_retained = length(pr$p) / sum(tep),
                   dist_min_m = min(d_p), dist_p05_m = unname(stats::quantile(d_p, 0.05)), dist_median_m = stats::median(d_p), converged = f$converged)
  })
}

cen_lab <- A$centre
centres <- setdiff(names(which(table(cen_lab[A$yy == 1L]) >= 30L)), "fuera_AOI")
if (QUICK) centres <- centres[1:2]
cv_all <- list(); loco_all <- list(); blocks_all <- list()
# bloques por repeticion (para el mapa de pliegues)
for (r in seq_len(n_rep)) {
  fo <- spatial_block_folds(A$x, A$y, block_m, k_cv, seed + r)
  b <- A |> dplyr::mutate(fold = fo) |> dplyr::group_by(blk, fold) |> dplyr::summarise(n_presences = sum(yy == 1L), n_background = sum(yy == 0L), .groups = "drop") |> dplyr::mutate(repeat_id = r)
  blocks_all[[r]] <- b
}
for (vn in names(VARIANTS)) {
  say("Variante ", vn, ": validacion espacial por bloques...")
  res <- purrr::map_dfr(seq_len(n_rep), function(r) {
    fo <- spatial_block_folds(A$x, A$y, block_m, k_cv, seed + r)
    purrr::map_dfr(seq_len(if (QUICK) 2L else k_cv), function(f) {
      x <- run_split(fo != f, fo == f, vn, seed + r * 10L + f); if (is.null(x)) NULL else dplyr::mutate(x, repeat_id = r, fold = f) })
  })
  cv_all[[vn]] <- res; write.csv(dplyr::bind_rows(cv_all), file.path(OUT, if (QUICK) "rev_cv_splits_quick.csv" else "rev_cv_splits.csv"), row.names = FALSE)
  say("Variante ", vn, ": dejando un centro fuera...")
  lo <- purrr::map_dfr(centres, function(cn) {
    te <- cen_lab == cn; tr <- !te & !(A$blk %in% unique(A$blk[te]))
    x <- run_split(tr, te, vn, seed); if (is.null(x)) NULL else dplyr::mutate(x, centre = cn, n_blocks_removed = length(unique(A$blk[te & A$yy == 1L])))
  })
  loco_all[[vn]] <- lo; write.csv(dplyr::bind_rows(loco_all), file.path(OUT, if (QUICK) "rev_loco_quick.csv" else "rev_loco.csv"), row.names = FALSE)
  say("Variante ", vn, " lista. Mejor cAUC_pds spl_igm_vias (media pliegues): ", round(mean(res$cauc_pds[res$model == "spl_igm_vias"]), 3))
}
cv <- dplyr::bind_rows(cv_all); lo <- dplyr::bind_rows(loco_all)
summ <- cv |> dplyr::group_by(variant, model) |> dplyr::summarise(dplyr::across(c(auc, auc_pds, cauc_pds, auc_null, boyce, ssb, n_train_p, n_train_d, n_pds_pairs, pds_share_retained, dist_min_m, dist_p05_m),
          list(mean = ~mean(.x), lo = ~min(.x), hi = ~max(.x)), .names = "{.col}__{.fn}"), n_splits = dplyr::n(), .groups = "drop")
write.csv(summ, file.path(OUT, if (QUICK) "rev_cv_summary_quick.csv" else "rev_cv_summary.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(blocks_all), file.path(OUT, "rev_fold_blocks.csv"), row.names = FALSE)
sep <- dplyr::bind_rows(cv |> dplyr::mutate(scheme = "spatial_block", unit = paste0("r", repeat_id, "_f", fold)), lo |> dplyr::mutate(scheme = "leave_one_centre_out", unit = centre)) |>
  dplyr::filter(model == "spl_igm_vias") |> dplyr::select(scheme, unit, variant, n_train_p, n_train_d, n_excl_buffer_p, n_excl_buffer_d, n_test_p, n_test_d, n_pds_pairs, pds_share_retained, dist_min_m, dist_p05_m, dist_median_m, ssb, ssb_pds)
write.csv(sep, file.path(OUT, if (QUICK) "rev_fold_separation_quick.csv" else "rev_fold_separation.csv"), row.names = FALSE)
# comprobacion de reproduccion (V0 frente a lo publicado)
v0 <- summ |> dplyr::filter(variant == "V0", model == "spl_igm_vias"); l0 <- lo |> dplyr::filter(variant == "V0", model == "spl_igm_vias")
txt <- c(sprintf("V0 spl_igm_vias cAUC_pds media = %.3f (manuscrito: 0.811)", v0$cauc_pds__mean), sprintf("V0 spl_igm_vias AUC media = %.3f", v0$auc__mean),
         paste("V0 LOCO cAUC_pds:", paste(l0$centre, round(l0$cauc_pds, 3), collapse = "; "), "(manuscrito: Catamayo 0.709, Loja 0.811, Malacatos 0.795, San Pedro 0.811, Vilcabamba 0.882)"))
writeLines(txt, file.path(OUT, if (QUICK) "rev_reproduccion_quick.txt" else "rev_reproduccion.txt")); cat(txt, sep = "\n")
say("Fin.")
