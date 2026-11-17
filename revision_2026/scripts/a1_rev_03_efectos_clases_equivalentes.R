# a1_rev_03_efectos_clases_equivalentes.R
# Objetivo : atender a los revisores 1 (unidades coherentes, comparacion del buffer geometrico) y 2 (clases viales equivalentes; la diferencia entre fuentes no puede atribuirse solo a la fuente):
#   A. razones de intensidad en unidades FIJAS (por 10 grados de pendiente; por duplicacion de la distancia a la via y al nucleo) para el dominio principal, con intervalo por bootstrap de bloques de 10 km (B = 300),
#      para las cuatro fuentes originales y para pares de clases equivalentes IGM/OSM: 'mayor' (IGM primaria+secundaria vs OSM trunk/primary/secondary), 'vehicular' (IGM via_l vs OSM sin sendas/caminos), 'todo' (IGM con sendas y roderas vs OSM completo);
#   B. las mismas unidades fijas en los cinco dominios (intervalos del modelo, que subestiman la incertidumbre, como en el manuscrito);
#   C. validacion espacial por bloques (como V0 del manuscrito) de cada conjunto, con la diferencia pareada OSM - IGM dentro de cada clase equivalente;
#   D. influencia de los bloques: ajuste dejando fuera cada uno de los 29 bloques ocupados (jackknife por bloques) y bootstrap con B = 1000 del conjunto principal (estabilidad del intervalo).
# Entradas : cache M03_dataset_ea8ce5b426.rds, data/processed/REV_dist_osm_vehicular.rds (a1_rev_02), outputs/tables/tab_M03_coefficients.csv (comprobacion)
# Salidas  : outputs/revision/rev_efectos_unidades_fijas.csv, rev_efectos_dominios.csv, rev_cv_clases_equivalentes.csv, rev_cv_pareada.csv, rev_jackknife_bloques.csv, rev_bootstrap_B1000.csv
# Fecha    : 2026-11-07
suppressMessages({library(dplyr); library(purrr); library(tibble); library(sf); library(splines)})
args <- commandArgs(trailingOnly = TRUE); QUICK <- "quick" %in% args
ROOT <- "C:/Trabajos_y_Publicaciones/02_cambio_ciudades/000_A1_research"; OUT <- file.path(ROOT, "outputs/revision")
say <- function(...) cat(sprintf("[rev03 %s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
seed <- 2026L; k_cv <- 5L; n_rep <- if (QUICK) 1L else 3L; block_m <- 10000; B <- if (QUICK) 5L else 300L; B_big <- if (QUICK) 10L else 1000L
AREAS <- c(extended = 1639.0470387658202, aoi_amplio = 756.3560183225472, aoi_estricto = 169.39033944050968, aoi_base = 432.5138803817142, core_buffer = 1309.5335926617217)

auc_rank <- function(pos, neg) { n1 <- as.numeric(length(pos)); n0 <- as.numeric(length(neg)); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(c(pos, neg)); (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
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
  ll_fun <- function(b) sum(Xp %*% b) - cc * sum(exp(Xd %*% b)); ll <- ll_fun(beta); conv <- FALSE
  for (it in seq_len(maxit)) {
    mu <- cc * exp(drop(Xd %*% beta)); g <- sp - drop(crossprod(Xd, mu)); H <- crossprod(Xd * mu, Xd)
    step <- tryCatch(solve(H, g), error = function(e) NULL); if (is.null(step)) break
    t <- 1; ok <- FALSE
    for (h in 1:25) { bn <- beta + t * step; lln <- ll_fun(bn); if (is.finite(lln) && lln >= ll - 1e-12) { ok <- TRUE; break }; t <- t / 2 }
    if (!ok) break
    d <- abs(lln - ll); beta <- bn; ll <- lln
    if (d < tol * (1 + abs(ll)) || max(abs(t * step)) < 1e-8) { conv <- TRUE; break } }
  mu <- cc * exp(drop(Xd %*% beta)); V <- tryCatch(solve(crossprod(Xd * mu, Xd)), error = function(e) matrix(NA_real_, k, k))
  list(beta = beta, ll = ll, se = sqrt(diag(V)), converged = conv) }
spatial_block_folds <- function(x, y, block_size, k, seed) {
  block_id <- paste(floor(x / block_size), floor(y / block_size), sep = "_"); blocks <- unique(block_id)
  set.seed(seed); fold_of_block <- stats::setNames(sample(rep_len(seq_len(k), length(blocks))), blocks); unname(fold_of_block[block_id]) }

# ---------------- datos ----------------
cache <- readRDS(file.path(ROOT, "data/processed/M03_dataset_ea8ce5b426.rds")); veh <- readRDS(file.path(ROOT, "data/processed/REV_dist_osm_vehicular.rds"))
P0 <- cache$pres; D0 <- cache$pool; P0 <- P0[is.finite(P0$slope), ]; D0 <- D0[is.finite(D0$slope), ]
P <- P0[P0$in_extended, ]; D <- D0[D0$in_extended, ]; D <- D[seq_len(ceiling(20 * nrow(P))), ]
stopifnot(nrow(P) == 17088L, nrow(D) == 341760L, isTRUE(all.equal(P$x, veh$P$x)), isTRUE(all.equal(D$x, veh$D$x)))
P$log_d_osm_veh <- veh$P$log_d_osm_veh; D$log_d_osm_veh <- veh$D$log_d_osm_veh
cc <- AREAS[["extended"]] / nrow(D)
# conjuntos: variable vial de cada uno (todos incluyen pendiente y distancia al nucleo)
SETS <- list(
  igm_vias = "log_d_igm_vias", igm_all = "log_d_igm_all", osm2019 = "log_d_osm2019", census2020 = "log_d_census2020",
  major_igm = "log_d_state_igm_primaria_secundaria", major_osm = "log_d_state_osm_primaria_secundaria",
  vehicular_igm = "log_d_igm_vias", vehicular_osm = "log_d_osm_veh", all_igm = "log_d_igm_all", all_osm = "log_d_osm2019")
CLASS_OF <- c(igm_vias = "fuente original", igm_all = "fuente original", osm2019 = "fuente original", census2020 = "fuente original",
              major_igm = "mayor", major_osm = "mayor", vehicular_igm = "vehicular", vehicular_osm = "vehicular", all_igm = "todo", all_osm = "todo")
if (QUICK) SETS <- SETS[c("igm_vias", "major_igm", "major_osm")]
allv <- unique(c("slope", "log_d_nucleo", unlist(SETS)))
cen <- vapply(allv, function(v) mean(D[[v]]), numeric(1)); sdv <- vapply(allv, function(v) stats::sd(D[[v]]), numeric(1))
# matriz en unidades fijas: pendiente por 10 grados; distancias por duplicacion (log1p(d)/log 2)
rawX <- function(T, road) cbind(`(Intercept)` = 1, slope_10deg = (T$slope - cen[["slope"]]) / 10, road_doubling = (T[[road]] - cen[[road]]) / log(2), core_doubling = (T$log_d_nucleo - cen[["log_d_nucleo"]]) / log(2))
fit_raw <- function(Pd, Dd, road, cc_) ipp_fit(rawX(Pd, road), rawX(Dd, road), cc_)
blk <- function(T) paste(floor(T$x / block_m), floor(T$y / block_m))
boot_fit <- function(Xp, Xd, bp, bd, beta0, B, sd0, cc_) {
  ub <- union(unique(bp), unique(bd)); spp <- split(seq_len(nrow(Xp)), factor(bp, levels = ub)); sdd <- split(seq_len(nrow(Xd)), factor(bd, levels = ub))
  out <- matrix(NA_real_, B, ncol(Xp), dimnames = list(NULL, colnames(Xp))); set.seed(sd0)
  for (b in seq_len(B)) { s <- sample(ub, length(ub), replace = TRUE); ip <- unlist(spp[s], use.names = FALSE); id <- unlist(sdd[s], use.names = FALSE)
    f <- tryCatch(ipp_fit(Xp[ip, , drop = FALSE], Xd[id, , drop = FALSE], cc_, beta0, maxit = 30L), error = function(e) NULL)
    if (!is.null(f) && isTRUE(f$converged)) out[b, ] <- f$beta }
  out }
bp <- blk(P); bd <- blk(D)

# ---------------- A. unidades fijas en el dominio principal con bootstrap de bloques ----------------
say("A. Efectos en unidades fijas (B = ", B, ")...")
eff <- purrr::imap_dfr(SETS, function(road, nm) {
  f <- fit_raw(P, D, road, cc); bt <- boot_fit(rawX(P, road), rawX(D, road), bp, bd, f$beta, B, seed + 7L, cc); ok <- stats::complete.cases(bt)
  vv <- c(slope_10deg = "slope_10deg", road_doubling = "road_doubling", core_doubling = "core_doubling")
  purrr::imap_dfr(vv, function(col, term) tibble::tibble(feature_set = nm, class = CLASS_OF[[nm]], term = term, ratio = exp(f$beta[[col]]),
    boot_lo = unname(stats::quantile(exp(bt[ok, col]), 0.025)), boot_hi = unname(stats::quantile(exp(bt[ok, col]), 0.975)),
    model_lo = exp(f$beta[[col]] - 1.96 * f$se[[col]]), model_hi = exp(f$beta[[col]] + 1.96 * f$se[[col]]), n_boot_ok = sum(ok),
    ratio_per_SD_check = exp(f$beta[[col]] * c(slope_10deg = sdv[["slope"]] / 10, road_doubling = sdv[[road]] / log(2), core_doubling = sdv[["log_d_nucleo"]] / log(2))[[term]])))
})
write.csv(eff, file.path(OUT, if (QUICK) "rev_efectos_unidades_fijas_quick.csv" else "rev_efectos_unidades_fijas.csv"), row.names = FALSE)
pub <- read.csv(file.path(ROOT, "outputs/tables/tab_M03_coefficients.csv")); pubm <- pub |> dplyr::filter(feature_set %in% names(SETS)[1:min(4, length(SETS))])
chk <- eff |> dplyr::filter(feature_set %in% pubm$feature_set) |> dplyr::mutate(variable = dplyr::case_when(term == "slope_10deg" ~ "slope", term == "core_doubling" ~ "dist_nucleo", TRUE ~ NA_character_)) |>
  dplyr::filter(!is.na(variable)) |> dplyr::left_join(pubm |> dplyr::select(feature_set, variable, pub = intensity_ratio, pub_lo = ir_boot_lo, pub_hi = ir_boot_hi), by = c("feature_set", "variable"))
say("Comprobacion con el manuscrito (razon por DE recalculada vs publicada; intervalo bootstrap recalculado por 10 grados / duplicacion, no comparable directamente):")
print(as.data.frame(chk |> dplyr::select(feature_set, term, ratio_per_SD_check, pub) |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))))
stopifnot(all(abs(chk$ratio_per_SD_check - chk$pub) < 0.0015))

# ---------------- B. unidades fijas por dominio (conjunto principal igm_vias; intervalos del modelo) ----------------
say("B. Dominios...")
dom_eff <- purrr::imap_dfr(AREAS, function(a, dn) {
  fl <- paste0("in_", dn); Pd <- P0[P0[[fl]], ]; Dd <- D0[D0[[fl]], ]
  if (dn == "extended") { Pd <- P; Dd <- D }
  cc_d <- a / nrow(Dd); f <- fit_raw(Pd, Dd, "log_d_igm_vias", cc_d)
  purrr::imap_dfr(c(slope_10deg = "slope_10deg", road_doubling = "road_doubling", core_doubling = "core_doubling"), function(col, term)
    tibble::tibble(domain = dn, n_presences = nrow(Pd), n_background = nrow(Dd), area_km2 = a, term = term, ratio = exp(f$beta[[col]]), model_lo = exp(f$beta[[col]] - 1.96 * f$se[[col]]), model_hi = exp(f$beta[[col]] + 1.96 * f$se[[col]])))
})
write.csv(dom_eff, file.path(OUT, if (QUICK) "rev_efectos_dominios_quick.csv" else "rev_efectos_dominios.csv"), row.names = FALSE)
print(as.data.frame(dom_eff |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))))

# ---------------- C. validacion espacial por bloques (V0) de cada conjunto, lineal y con splines ----------------
say("C. Validacion espacial por conjunto...")
A <- dplyr::bind_rows(P |> dplyr::mutate(yy = 1L), D |> dplyr::mutate(yy = 0L)); A <- as.data.frame(A[, c("x", "y", "yy", allv)])
Zs <- A; for (v in allv) Zs[[v]] <- (A[[v]] - cen[[v]]) / sdv[[v]]
bases <- lapply(stats::setNames(allv, allv), function(v) splines::ns(Zs[[v]][Zs$yy == 0L], df = 3L))
Xmat <- function(idx, road, ns) { if (ns) { m <- function(v) { z <- stats::predict(bases[[v]], Zs[[v]][idx]); colnames(z) <- paste0(v, seq_len(ncol(z))); z }; cbind(`(Intercept)` = 1, m("slope"), m(road), m("log_d_nucleo")) }
  else cbind(`(Intercept)` = 1, slope = Zs$slope[idx], road = Zs[[road]][idx], core = Zs$log_d_nucleo[idx]) }
cvres <- purrr::map_dfr(seq_len(n_rep), function(r) {
  fo <- spatial_block_folds(A$x, A$y, block_m, k_cv, seed + r)
  purrr::map_dfr(seq_len(k_cv), function(f) {
    te <- fo == f; tr <- !te; trp <- which(tr & A$yy == 1L); trd <- which(tr & A$yy == 0L); tep <- which(te & A$yy == 1L); ted <- which(te & A$yy == 0L)
    g_tr <- sf::st_as_sf(A[trp, c("x", "y")], coords = c("x", "y"), crs = 32717)
    near <- function(idx) { g <- sf::st_as_sf(A[idx, c("x", "y")], coords = c("x", "y"), crs = 32717); as.numeric(sf::st_distance(g, g_tr[sf::st_nearest_feature(g, g_tr), ], by_element = TRUE)) }
    d_p <- near(tep); d_d <- near(ted); pr <- pds_idx(d_p, d_d, seed + r * 10L + f)
    purrr::imap_dfr(SETS, function(road, nm) purrr::map_dfr(c(FALSE, TRUE), function(ns) {
      fi <- ipp_fit(Xmat(trp, road, ns), Xmat(trd, road, ns), cc); ep <- drop(Xmat(tep, road, ns) %*% fi$beta); ed <- drop(Xmat(ted, road, ns) %*% fi$beta)
      tibble::tibble(feature_set = nm, class = CLASS_OF[[nm]], spline = ns, repeat_id = r, fold = f, auc = auc_rank(ep, ed), auc_pds = auc_rank(ep[pr$p], ed[pr$d])) }))
  }) })
write.csv(cvres, file.path(OUT, if (QUICK) "rev_cv_clases_equivalentes_quick.csv" else "rev_cv_clases_equivalentes.csv"), row.names = FALSE)
cvs <- cvres |> dplyr::group_by(feature_set, class, spline) |> dplyr::summarise(auc_mean = mean(auc), auc_pds_mean = mean(auc_pds), auc_pds_min = min(auc_pds), auc_pds_max = max(auc_pds), n_splits = dplyr::n(), .groups = "drop")
print(as.data.frame(cvs |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))))
pair <- purrr::map_dfr(c("mayor", "vehicular", "todo"), function(cl) purrr::map_dfr(c(FALSE, TRUE), function(sp) {
  nm <- switch(cl, mayor = c("major_igm", "major_osm"), vehicular = c("vehicular_igm", "vehicular_osm"), todo = c("all_igm", "all_osm"))
  a <- cvres |> dplyr::filter(feature_set == nm[1], spline == sp) |> dplyr::arrange(repeat_id, fold); b <- cvres |> dplyr::filter(feature_set == nm[2], spline == sp) |> dplyr::arrange(repeat_id, fold)
  if (nrow(a) == 0 || nrow(b) == 0) return(NULL); d <- b$auc_pds - a$auc_pds
  tibble::tibble(class = cl, spline = sp, igm = nm[1], osm = nm[2], mean_igm = mean(a$auc_pds), mean_osm = mean(b$auc_pds), mean_diff_osm_minus_igm = mean(d), diff_min = min(d), diff_max = max(d), n_folds_osm_higher = sum(d > 0), n_folds = length(d)) }))
write.csv(pair, file.path(OUT, if (QUICK) "rev_cv_pareada_quick.csv" else "rev_cv_pareada.csv"), row.names = FALSE); print(as.data.frame(pair |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))))

# ---------------- D. jackknife por bloques y bootstrap con B grande (conjunto principal) ----------------
say("D. Jackknife por bloques (conjunto principal)...")
f0 <- fit_raw(P, D, "log_d_igm_vias", cc); ub <- sort(unique(bp)); npb <- table(factor(bp, levels = ub))
jk <- purrr::map_dfr(ub, function(b) { ip <- bp != b; id <- bd != b; f <- ipp_fit(rawX(P[ip, ], "log_d_igm_vias"), rawX(D[id, ], "log_d_igm_vias"), cc, f0$beta)
  tibble::tibble(block_dropped = b, n_presences_dropped = as.integer(npb[[b]]), slope_10deg = exp(f$beta[["slope_10deg"]]), road_doubling = exp(f$beta[["road_doubling"]]), core_doubling = exp(f$beta[["core_doubling"]])) })
write.csv(jk, file.path(OUT, if (QUICK) "rev_jackknife_bloques_quick.csv" else "rev_jackknife_bloques.csv"), row.names = FALSE)
say("Rango jackknife por 10 grados: ", paste(round(range(jk$slope_10deg), 3), collapse = "-"), "; via: ", paste(round(range(jk$road_doubling), 3), collapse = "-"), "; nucleo: ", paste(round(range(jk$core_doubling), 3), collapse = "-"))
say("Bootstrap B = ", B_big, " del conjunto principal...")
bt <- boot_fit(rawX(P, "log_d_igm_vias"), rawX(D, "log_d_igm_vias"), bp, bd, f0$beta, B_big, seed + 7L, cc); ok <- stats::complete.cases(bt)
bb <- purrr::map_dfr(c("slope_10deg", "road_doubling", "core_doubling"), function(col) tibble::tibble(term = col, B = B_big, n_ok = sum(ok), ratio = exp(f0$beta[[col]]), lo = unname(stats::quantile(exp(bt[ok, col]), 0.025)), hi = unname(stats::quantile(exp(bt[ok, col]), 0.975)),
  lo_first300 = unname(stats::quantile(exp(bt[ok, col][seq_len(min(300, sum(ok)))]), 0.025)), hi_first300 = unname(stats::quantile(exp(bt[ok, col][seq_len(min(300, sum(ok)))]), 0.975))))
write.csv(bb, file.path(OUT, if (QUICK) "rev_bootstrap_B1000_quick.csv" else "rev_bootstrap_B1000.csv"), row.names = FALSE); print(as.data.frame(bb |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))))
say("Fin.")
