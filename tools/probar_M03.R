# a1_probar_M03.R
# Objetivo : ejecutar el modulo M03_presencia_fondo FUERA de las puertas del pipeline (etapas 01-05) para probarlo con datos reales.
#            No toca el protocolo congelado ni los resultados oficiales: escribe en outputs/method_test/.
# Uso      : Rscript a1_probar_M03.R <modo: rapido|completo>
#            rapido   = fondo reducido (testing_scale 0.1), 1 repeticion, bootstrap B = 20 (comprueba que el codigo corre)
#            completo = parametros de project_card.yaml tal cual (resultado exploratorio reproducible, NO congelado)
# Entradas : project_card.yaml (bloque method_params.M03_presencia_fondo), data/raw/**
# Salidas  : outputs/method_test/<modo>/metrics_long.csv, tablas, figuras, method_facts.yaml
# Fecha    : 2026-10-20
args <- commandArgs(trailingOnly = TRUE)
modo <- if (length(args) >= 1L) args[1] else "rapido"
ROOT <- normalizePath(Sys.getenv("A1_ROOT", unset = getwd()), winslash = "/")
setwd(ROOT)
suppressMessages(source(file.path(ROOT, "00_setup.R"), encoding = "UTF-8"))
env <- new.env(parent = globalenv())
sys.source(file.path(ROOT, "methods", "M03_presencia_fondo.R"), envir = env, keep.source = FALSE)
p <- card_get(CARD, "method_params.M03_presencia_fondo")
msg <- env$validate_params(p)
if (length(msg) > 0L) stop("validate_params: ", paste(msg, collapse = "; "))
out_dir <- ensure_dir(file.path(ROOT, "outputs", "method_test", modo))
if (modo == "rapido") {
  p$testing_scale <- 0.1
  p$cv$repeats <- 1L
  p$bootstrap$B <- 20L
  p$cv$block_sizes_sensitivity <- c(4000, 8000)
}
require_pkgs(env$MODULE_INFO$packages)
ctx <- list(params = p, card = CARD, DIR = DIR, root = ROOT, seed = as.integer(card_get(CARD, "project.random_seed")), crs = card_get(CARD, "project.project_crs"))
t0 <- Sys.time()
res <- env$run_method(ctx)
cat("\nTiempo (min): ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
write_csv_safe(res$metrics, file.path(out_dir, "metrics_long.csv"))
for (id in names(res$tables)) write_csv_safe(tibble::as_tibble(res$tables[[id]]$data), file.path(out_dir, paste0(id, ".csv")))
for (fg in res$figures) save_fig(fg$plot, fg$id, fg$width_cm, fg$height_cm, dpi = 150, dir = file.path(out_dir, "figures"))
yaml::write_yaml(list(facts = res$facts, limitations = res$limitations), file.path(out_dir, "method_facts.yaml"))
cat("\nMetricas:", nrow(res$metrics), "| tablas:", length(res$tables), "| figuras:", length(res$figures), "\nSalida:", out_dir, "\n")
options(width = 200)
print(as.data.frame(res$metrics[res$metrics$subset %in% names(p$feature_sets) | grepl("^center", res$metrics$subset), c("metric", "subset", "value", "ci_low", "ci_high")]), digits = 3)
