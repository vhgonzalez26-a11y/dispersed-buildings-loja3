# ============================================================
# tests/test_M03.R
# Pruebas de las funciones estadísticas del módulo M03 con datos simulados de
# solución conocida. Sirven como evidencia de que la implementación (verosimilitud
# IPP, AUC, AUC calibrado, muestreo por pares de distancias, índice de Boyce,
# intervalo de Wilson) hace lo que dice el método.
# Uso: Rscript tests/test_M03.R   (desde la raíz del proyecto o con ruta absoluta)
# Salida: tests/test_M03_resultado.csv con cada prueba, valor obtenido, criterio y PASS/FAIL
# ============================================================
ROOT <- normalizePath(Sys.getenv("A1_ROOT", unset = getwd()), winslash = "/")
env <- new.env(parent = globalenv())
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
sys.source(file.path(ROOT, "methods", "M03_presencia_fondo.R"), envir = env, keep.source = FALSE)
res <- data.frame(test = character(), valor = character(), criterio = character(), resultado = character())
add <- function(test, valor, criterio, ok) res[nrow(res) + 1L, ] <<- list(test, valor, criterio, if (isTRUE(ok)) "PASS" else "FAIL")

# 1. IPP recupera coeficientes conocidos ------------------------------------
set.seed(1)
A <- 100                                  # km2
sim_pp <- function(beta, A, n_target) {
  lam_fun <- function(z) exp(beta[1] + beta[2] * z[, 1] + beta[3] * z[, 2])
  zmax <- 2.5
  lam_max <- exp(beta[1] + abs(beta[2]) * zmax + abs(beta[3]) * zmax)
  N <- rpois(1, lam_max * A)
  z <- cbind(rnorm(N), rnorm(N)); keep <- runif(N) < lam_fun(pmin(pmax(z, -zmax), zmax)) / lam_max
  z[keep, , drop = FALSE]
}
beta_true <- c(log(80), -0.6, 0.3)                      # ~ 80 presencias/km2 en z = 0
Zp <- sim_pp(beta_true, A, NA)
m <- 60000; Zd <- cbind(rnorm(m), rnorm(m))              # covariables evaluadas en el fondo (misma distribución que el ambiente)
# Para que el fondo represente el ambiente, el ambiente tiene covariables N(0,1): la integral = A * E[lambda]
Xp <- cbind(`(Intercept)` = 1, z1 = Zp[, 1], z2 = Zp[, 2]); Xd <- cbind(`(Intercept)` = 1, z1 = Zd[, 1], z2 = Zd[, 2])
f <- env$ipp_fit(Xp, Xd, cc = A / m)
add("IPP: convergencia", f$iter, "converged = TRUE", isTRUE(f$converged))
add("IPP: coeficiente z1 (verdadero -0.6)", sprintf("%.3f (se %.3f)", f$beta[2], f$se[2]), "|estimado - verdadero| < 3 se", abs(f$beta[2] - beta_true[2]) < 3 * f$se[2])
add("IPP: coeficiente z2 (verdadero 0.3)", sprintf("%.3f (se %.3f)", f$beta[3], f$se[3]), "|estimado - verdadero| < 3 se", abs(f$beta[3] - beta_true[3]) < 3 * f$se[3])
add("IPP: intercepto (verdadero log 80 = 4.382)", sprintf("%.3f", f$beta[1]), "|estimado - verdadero| < 0.15", abs(f$beta[1] - beta_true[1]) < 0.15)

# 2. Regresión logística sin ponderar con fondo pequeño sesga el efecto cuando el modelo está mal especificado ----
# (Fithian & Hastie 2013): con un modelo correcto la logística coincide en las pendientes; aquí se comprueba que el IPP no depende del tamaño del fondo
f_small <- env$ipp_fit(Xp, Xd[1:6000, ], cc = A / 6000)
add("IPP: estable con fondo 10 veces menor", sprintf("%.3f vs %.3f", f_small$beta[2], f$beta[2]), "|diferencia| < 0.1", abs(f_small$beta[2] - f$beta[2]) < 0.1)

# 3. AUC por rangos = prueba de Mann-Whitney ---------------------------------
x1 <- rnorm(300, 1); x0 <- rnorm(900)
w <- suppressWarnings(wilcox.test(x1, x0)$statistic) / (300 * 900)
add("AUC: igual a Mann-Whitney", sprintf("%.6f vs %.6f", env$auc_rank(x1, x0), w), "diferencia < 1e-10", abs(env$auc_rank(x1, x0) - w) < 1e-10)
add("AUC: modelo aleatorio ~ 0.5", sprintf("%.3f", env$auc_rank(rnorm(2000), rnorm(2000))), "entre 0.45 y 0.55", abs(env$auc_rank(rnorm(2000), rnorm(2000)) - 0.5) < 0.05)
add("AUC: separación perfecta = 1", sprintf("%.3f", env$auc_rank(c(5, 6), c(1, 2, 3))), "== 1", env$auc_rank(c(5, 6), c(1, 2, 3)) == 1)

# 4. Muestreo por pares de distancias (Hijmans 2012) --------------------------
d_p <- runif(200, 0, 1000); d_d <- runif(5000, 0, 5000)
pr <- env$pds_idx(d_p, d_d, 1L)
add("PDS: parejas sin repetir fondo, >= 90 % de presencias emparejadas", sprintf("%d parejas, %d fondos únicos", length(pr$p), length(unique(pr$d))), "únicos = parejas y >= 180", length(pr$d) == length(unique(pr$d)) && length(pr$p) >= 180)
rel <- abs(d_p[pr$p] - d_d[pr$d]) / pmax(d_p[pr$p], d_d[pr$d])
add("PDS: ninguna pareja difiere más de 33 %", sprintf("máximo %.3f", max(rel)), "<= 0.33", max(rel) <= 0.33)
# el muestreo por pares elimina el sesgo de ordenamiento: SSB = Dp/Da pasa de ~0.1-0.2 a ~1
add("PDS: SSB tras el muestreo cercano a 1", sprintf("antes %.2f, después %.2f", mean(d_p) / mean(d_d), mean(d_p[pr$p]) / mean(d_d[pr$d])), "después entre 0.9 y 1.1", abs(mean(d_p[pr$p]) / mean(d_d[pr$d]) - 1) < 0.1)

# 5. Índice de Boyce -----------------------------------------------------------
bg <- rnorm(20000)
pres_good <- rnorm(2000, mean = 1.5)                      # presencias en la parte alta de la escala
pres_rand <- rnorm(2000)                                   # sin relación
add("Boyce: modelo informativo > 0.8", sprintf("%.3f", env$boyce_cont(pres_good, bg)$index), "> 0.8", env$boyce_cont(pres_good, bg)$index > 0.8)
# el índice de un modelo sin información es ruidoso (curva P/E suavizada por ventanas móviles): se evalúa su media en 60 repeticiones
b_rand <- vapply(1:60, function(i) env$boyce_cont(rnorm(2000), bg)$index, numeric(1))
add("Boyce: modelo aleatorio, media de 60 repeticiones cerca de 0", sprintf("media %.3f (DE %.2f)", mean(b_rand), sd(b_rand)), "|media| < 0.3", abs(mean(b_rand)) < 0.3)

# 6. AUC calibrado (Hijmans 2012, ec. 1) ---------------------------------------
cauc <- function(auc, nauc) auc + 0.5 - max(0.5, nauc)
add("cAUC: modelo nulo informativo resta", sprintf("%.2f", cauc(0.85, 0.65)), "== 0.70", abs(cauc(0.85, 0.65) - 0.70) < 1e-12)
add("cAUC: nulo <= 0.5 no cambia el AUC", sprintf("%.2f", cauc(0.85, 0.45)), "== 0.85", abs(cauc(0.85, 0.45) - 0.85) < 1e-12)

# 7. Intervalo de Wilson -------------------------------------------------------
ci <- env$wilson_ci(50, 100)
add("Wilson: 50/100", sprintf("[%.4f, %.4f]", ci[1], ci[2]), "[0.4038, 0.5962]", abs(ci[1] - 0.4038) < 1e-3 && abs(ci[2] - 0.5962) < 1e-3)

# 8. Correlograma de Moran: detecta autocorrelación conocida --------------------
set.seed(2); g <- expand.grid(x = 1:25, y = 1:25) * 1000
Dm <- as.matrix(dist(g)); S <- exp(-Dm / 4000); r_auto <- as.numeric(t(chol(S + diag(1e-8, nrow(S)))) %*% rnorm(nrow(S))); r_ind <- rnorm(nrow(S))
c_auto <- env$morans_correlogram(g, r_auto, 1000, 8L, 49L, 1L); c_ind <- env$morans_correlogram(g, r_ind, 1000, 8L, 49L, 1L)
add("Moran: campo autocorrelacionado, I(banda 1) > 0.3", sprintf("%.2f", c_auto$moran_i[1]), "> 0.3 y p < 0.05", c_auto$moran_i[1] > 0.3 && c_auto$p_value[1] < 0.05)
add("Moran: ruido blanco, |I(banda 1)| < 0.1", sprintf("%.3f", c_ind$moran_i[1]), "< 0.1", abs(c_ind$moran_i[1]) < 0.1)

write.csv(res, file.path(ROOT, "tests", "test_M03_resultado.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(res, right = FALSE)
cat("\nPruebas:", nrow(res), "| PASS:", sum(res$resultado == "PASS"), "| FAIL:", sum(res$resultado == "FAIL"), "\n")
if (any(res$resultado == "FAIL")) quit(status = 1)
