# a1_rev_04_sensibilidad_tipo_edificio.R
# Objetivo : (revisor 2, comentario 3) valorar si tratar como 'dispersos' todos los edificios fuera de manzanas condiciona los resultados. Se reajusta el modelo principal (IGM sin sendas) en submuestras:
#   S0 todos los edificios (referencia, dominio principal);
#   S1 sin los edificios con uso no residencial rotulado en el censo (campo cod_otros);
#   S2 sin los edificios dentro de los nucleos urbanos IGM ni del limite urbano oficial 2021 (Loja, Malacatos, Vilcabamba), y sin el fondo de esas zonas (dominio 'rural');
#   S3a / S3b sin los edificios a menos de 1 km / 2 km del nucleo urbano y sin el fondo en esas franjas (excluye la franja periurbana).
#   Los efectos se dan en unidades fijas (por 10 grados de pendiente; por duplicacion de distancia) con bootstrap de bloques de 10 km (B = 300) y se valida con la validacion espacial por bloques (5 x 3).
# Entradas : cache M03_dataset_ea8ce5b426.rds; data/raw/vector/censo_cartografia_aoi.gpkg (cod_otros); limite_urbano_2021.shp; igm_50k_2013_base.gpkg (zona_urbana_a)
# Salidas  : outputs/revision/rev_sensibilidad_tipo_edificio.csv, rev_sensibilidad_conteos.csv
# Fecha    : 2026-11-07
suppressMessages({library(dplyr); library(purrr); library(tibble); library(sf)})
args <- commandArgs(trailingOnly = TRUE); QUICK <- "quick" %in% args
ROOT <- "C:/Trabajos_y_Publicaciones/02_cambio_ciudades/000_A1_research"; OUT <- file.path(ROOT, "outputs/revision")
say <- function(...) cat(sprintf("[rev04 %s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
seed <- 2026L; k_cv <- 5L; n_rep <- if (QUICK) 1L else 3L; block_m <- 10000; B <- if (QUICK) 5L else 300L; CRS <- 32717
AREA <- 1639.0470387658202

auc_rank <- function(pos, neg) { n1 <- as.numeric(length(pos)); n0 <- as.numeric(length(neg)); if (n1 == 0 || n0 == 0) return(NA_real_); r <- rank(c(pos, neg)); (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
pds_idx <- function(d_p, d_d, seed = 1L, tol = 0.33) {
  n <- length(d_d); ob <- order(d_d); srt <- d_d[ob]; avail <- rep(TRUE, n); set.seed(seed); ip <- integer(length(d_p)); idd <- integer(length(d_p)); k <- 0L
  for (i in sample.int(length(d_p))) { pos <- findInterval(d_p[i], srt); lo <- pos; hi <- pos + 1L
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
  for (it in seq_len(maxit)) { mu <- cc * exp(drop(Xd %*% beta)); g <- sp - drop(crossprod(Xd, mu)); H <- crossprod(Xd * mu, Xd)
    step <- tryCatch(solve(H, g), error = function(e) NULL); if (is.null(step)) break
    t <- 1; ok <- FALSE
    for (h in 1:25) { bn <- beta + t * step; lln <- ll_fun(bn); if (is.finite(lln) && lln >= ll - 1e-12) { ok <- TRUE; break }; t <- t / 2 }
    if (!ok) break
    d <- abs(lln - ll); beta <- bn; ll <- lln; if (d < tol * (1 + abs(ll)) || max(abs(t * step)) < 1e-8) { conv <- TRUE; break } }
  mu <- cc * exp(drop(Xd %*% beta)); V <- tryCatch(solve(crossprod(Xd * mu, Xd)), error = function(e) matrix(NA_real_, k, k))
  list(beta = beta, ll = ll, se = sqrt(diag(V)), converged = conv) }
spatial_block_folds <- function(x, y, block_size, k, seed) { block_id <- paste(floor(x / block_size), floor(y / block_size), sep = "_"); blocks <- unique(block_id)
  set.seed(seed); fo <- stats::setNames(sample(rep_len(seq_len(k), length(blocks))), blocks); unname(fo[block_id]) }
boot_fit <- function(Xp, Xd, bp, bd, beta0, B, sd0, cc_) {
  ub <- union(unique(bp), unique(bd)); spp <- split(seq_len(nrow(Xp)), factor(bp, levels = ub)); sdd <- split(seq_len(nrow(Xd)), factor(bd, levels = ub))
  out <- matrix(NA_real_, B, ncol(Xp), dimnames = list(NULL, colnames(Xp))); set.seed(sd0)
  for (b in seq_len(B)) { s <- sample(ub, length(ub), replace = TRUE); ip <- unlist(spp[s], use.names = FALSE); id <- unlist(sdd[s], use.names = FALSE)
    f <- tryCatch(ipp_fit(Xp[ip, , drop = FALSE], Xd[id, , drop = FALSE], cc_, beta0, maxit = 30L), error = function(e) NULL); if (!is.null(f) && isTRUE(f$converged)) out[b, ] <- f$beta }
  out }

cache <- readRDS(file.path(ROOT, "data/processed/M03_dataset_ea8ce5b426.rds"))
P0 <- cache$pres; D0 <- cache$pool; P0 <- P0[is.finite(P0$slope), ]; D0 <- D0[is.finite(D0$slope), ]
P <- P0[P0$in_extended, ]; D <- D0[D0$in_extended, ]; D <- D[seq_len(ceiling(20 * nrow(P))), ]; stopifnot(nrow(P) == 17088L, nrow(D) == 341760L)
cc <- AREA / nrow(D)

# uso rotulado en el censo (uniendo por coordenadas) y zonas urbanas
say("Uniendo el uso rotulado de los edificios...")
v <- function(...) file.path(ROOT, "data/raw/vector", ...)
b <- st_read(v("censo_cartografia_aoi.gpkg"), layer = "edificios_2021_puntos", quiet = TRUE); b <- b[b$en_manzana == 0, ]; b <- st_transform(b, CRS)
bx <- st_coordinates(b); key <- function(x, y) paste(round(x, 1), round(y, 1))
lab <- data.frame(key = key(bx[, 1], bx[, 2]), nonres = nzchar(trimws(as.character(b$cod_otros))), stringsAsFactors = FALSE) |> dplyr::distinct(key, .keep_all = TRUE)
P$nonres <- lab$nonres[match(key(P$x, P$y), lab$key)]
say("Edificios con rotulo no residencial en el dominio principal: ", sum(P$nonres, na.rm = TRUE), " (sin cruce: ", sum(is.na(P$nonres)), ")")
P$nonres[is.na(P$nonres)] <- FALSE
ul <- st_union(st_transform(st_make_valid(st_read(v("limite_urbano_2021.shp"), quiet = TRUE)), CRS))
in_ul <- function(x, y) lengths(st_intersects(st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = CRS), ul)) > 0
P$in_ul <- in_ul(P$x, P$y); D$in_ul <- in_ul(D$x, D$y)
P$in_core <- P$d_nucleo == 0; D$in_core <- D$d_nucleo == 0

SUBS <- list(
  S0_todos = list(p = rep(TRUE, nrow(P)), d = rep(TRUE, nrow(D))),
  S1_sin_no_residencial = list(p = !P$nonres, d = rep(TRUE, nrow(D))),
  S2_sin_zonas_urbanas = list(p = !(P$in_core | P$in_ul), d = !(D$in_core | D$in_ul)),
  S3a_mas_de_1km_del_nucleo = list(p = P$d_nucleo > 1000, d = D$d_nucleo > 1000),
  S3b_mas_de_2km_del_nucleo = list(p = P$d_nucleo > 2000, d = D$d_nucleo > 2000))
cnt <- purrr::imap_dfr(SUBS, function(s, nm) tibble::tibble(subset = nm, n_presences = sum(s$p), n_removed_presences = nrow(P) - sum(s$p), n_background = sum(s$d), area_km2 = sum(s$d) * cc,
  n_non_residential_label = sum(P$nonres & s$p), n_in_core = sum(P$in_core & s$p), n_in_urban_limit = sum(P$in_ul & s$p)))
write.csv(cnt, file.path(OUT, if (QUICK) "rev_sensibilidad_conteos_quick.csv" else "rev_sensibilidad_conteos.csv"), row.names = FALSE); print(as.data.frame(cnt))

cen <- c(slope = mean(D$slope), road = mean(D$log_d_igm_vias), core = mean(D$log_d_nucleo))
rawX <- function(T, core = TRUE) { m <- cbind(`(Intercept)` = 1, slope_10deg = (T$slope - cen[["slope"]]) / 10, road_doubling = (T$log_d_igm_vias - cen[["road"]]) / log(2)); if (core) cbind(m, core_doubling = (T$log_d_nucleo - cen[["core"]]) / log(2)) else m }
blk <- function(T) paste(floor(T$x / block_m), floor(T$y / block_m))
res <- list(); cvl <- list()
for (nm in names(SUBS)) {
  s <- SUBS[[nm]]; Pp <- P[s$p, ]; Dd <- D[s$d, ]; say(nm, ": ", nrow(Pp), " presencias, ", nrow(Dd), " fondo")
  # en S3 la distancia al nucleo casi no varia dentro de la franja excluida; se mantiene en el modelo para comparar con el principal
  f <- ipp_fit(rawX(Pp), rawX(Dd), cc); bt <- boot_fit(rawX(Pp), rawX(Dd), blk(Pp), blk(Dd), f$beta, B, seed + 7L, cc); ok <- stats::complete.cases(bt)
  res[[nm]] <- purrr::map_dfr(c("slope_10deg", "road_doubling", "core_doubling"), function(col) tibble::tibble(subset = nm, term = col, ratio = exp(f$beta[[col]]), boot_lo = unname(stats::quantile(exp(bt[ok, col]), 0.025)), boot_hi = unname(stats::quantile(exp(bt[ok, col]), 0.975)), n_boot_ok = sum(ok)))
  Aa <- dplyr::bind_rows(Pp |> dplyr::mutate(yy = 1L), Dd |> dplyr::mutate(yy = 0L)); Aa <- as.data.frame(Aa[, c("x", "y", "yy", "slope", "log_d_igm_vias", "log_d_nucleo")])
  Xa <- rawX(Aa)
  cvl[[nm]] <- purrr::map_dfr(seq_len(n_rep), function(r) { fo <- spatial_block_folds(Aa$x, Aa$y, block_m, k_cv, seed + r)
    purrr::map_dfr(seq_len(k_cv), function(fl) { te <- fo == fl; trp <- which(!te & Aa$yy == 1L); trd <- which(!te & Aa$yy == 0L); tep <- which(te & Aa$yy == 1L); ted <- which(te & Aa$yy == 0L)
      if (length(tep) < 5L || length(ted) < 5L || length(trp) < 20L) return(NULL)
      g_tr <- st_as_sf(Aa[trp, c("x", "y")], coords = c("x", "y"), crs = CRS); near <- function(idx) { g <- st_as_sf(Aa[idx, c("x", "y")], coords = c("x", "y"), crs = CRS); as.numeric(st_distance(g, g_tr[st_nearest_feature(g, g_tr), ], by_element = TRUE)) }
      d_p <- near(tep); d_d <- near(ted); pr <- pds_idx(d_p, d_d, seed + r * 10L + fl)
      fi <- ipp_fit(Xa[trp, ], Xa[trd, ], cc); ep <- drop(Xa[tep, ] %*% fi$beta); ed <- drop(Xa[ted, ] %*% fi$beta)
      tibble::tibble(subset = nm, repeat_id = r, fold = fl, auc = auc_rank(ep, ed), auc_pds = auc_rank(ep[pr$p], ed[pr$d]), n_test_p = length(tep)) }) })
  write.csv(dplyr::bind_rows(res), file.path(OUT, if (QUICK) "rev_sensibilidad_tipo_edificio_quick.csv" else "rev_sensibilidad_tipo_edificio.csv"), row.names = FALSE)
}
eff <- dplyr::bind_rows(res); cv <- dplyr::bind_rows(cvl) |> dplyr::group_by(subset) |> dplyr::summarise(auc_mean = mean(auc), auc_pds_mean = mean(auc_pds), auc_pds_min = min(auc_pds), auc_pds_max = max(auc_pds), n_splits = dplyr::n(), .groups = "drop")
out <- eff |> dplyr::left_join(cv, by = "subset"); write.csv(out, file.path(OUT, if (QUICK) "rev_sensibilidad_tipo_edificio_quick.csv" else "rev_sensibilidad_tipo_edificio.csv"), row.names = FALSE)
print(as.data.frame(out |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))))); say("Fin.")
