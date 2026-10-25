# ============================================================
# 06_ejecutar_metodo.R
# Ejecuta el módulo de Materiales y Métodos definido en el protocolo
# congelado y valida que entregue el contrato común de salidas.
# Salidas : outputs/method/metrics_long.csv, method_facts.yaml,
#           outputs/tables/*.csv, outputs/figures/*.svg|png
# ============================================================

if (!exists("ROOT")) {
  ROOT <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), error = function(e) NULL)
  if (is.null(ROOT) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ROOT <- dirname(rstudioapi::getActiveDocumentContext()$path)
  }
  if (is.null(ROOT) || !nzchar(ROOT)) ROOT <- getwd()
}
source(file.path(ROOT, "00_setup.R"), encoding = "UTF-8")

cat("\n== ETAPA 06: EJECUCIÓN DEL MÓDULO DE MÉTODOS ==\n\n")
frozen <- assert_protocol_unchanged(CARD)
module_id <- frozen$payload$method_module
module_file <- file.path(DIR$methods, paste0(module_id, ".R"))
assert_file(module_file, "Módulo de métodos")

module_env <- new.env(parent = globalenv())
sys.source(module_file, envir = module_env, keep.source = FALSE)
MODULE_INFO <- module_env$MODULE_INFO
for (st in MODULE_INFO$requires_stages) require_stage(st, "06_metodo")
require_pkgs(MODULE_INFO$packages)

ctx <- list(
  params = frozen$payload$method_params,
  card = CARD, DIR = DIR, root = ROOT,
  seed = as.integer(frozen$payload$random_seed),
  crs = frozen$payload$project_crs
)

# ------------------------------------------------------------
# 1. Ejecutar
# ------------------------------------------------------------

set.seed(ctx$seed)
t0 <- Sys.time()
result <- module_env$run_method(ctx)
runtime_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
cat("Módulo ejecutado en ", runtime_min, " min.\n", sep = "")

# ------------------------------------------------------------
# 2. Validar contrato
# ------------------------------------------------------------

problems <- character()
if (!is.list(result)) problems <- c(problems, "run_method() debe devolver una lista")
if (!is.data.frame(result$metrics)) problems <- c(problems, "falta metrics (data.frame)")
if (is.data.frame(result$metrics)) {
  miss <- setdiff(c("metric", "value", "ci_low", "ci_high", "subset", "unit", "note"), names(result$metrics))
  if (length(miss) > 0L) problems <- c(problems, paste0("metrics sin columnas: ", paste(miss, collapse = ", ")))
  dup <- result$metrics |> dplyr::count(metric, subset) |> dplyr::filter(n > 1L)
  if (nrow(dup) > 0L) problems <- c(problems, paste0("métricas duplicadas: ", paste(dup$metric, dup$subset, collapse = "; ")))
}
if (!is.list(result$facts) || length(result$facts) == 0L) problems <- c(problems, "falta facts (hechos del método)")
for (fg in result$figures %||% list()) {
  if (is.null(fg$id) || !inherits(fg$plot, "ggplot") || is.null(fg$caption_draft)) problems <- c(problems, "figura sin id/plot/caption_draft")
}
if (length(problems) > 0L) stop("El módulo no cumple el contrato:\n- ", paste(problems, collapse = "\n- "), call. = FALSE)

# Hipótesis evaluables con lo producido
hyp_check <- purrr::map_dfr(frozen$payload$hypotheses, function(h) {
  subs <- c(clean_chr(h$subset), clean_chr(h$compare_subset %||% ""))
  subs <- subs[nzchar(subs)]
  found <- vapply(subs, function(s) any(result$metrics$metric == h$metric & result$metrics$subset == s), logical(1))
  tibble::tibble(id = h$id, metric = h$metric, subsets = paste(subs, collapse = " vs "), evaluable = all(found))
})
if (any(!hyp_check$evaluable)) {
  warning("Hipótesis no evaluables con las métricas producidas: ", paste(hyp_check$id[!hyp_check$evaluable], collapse = ", "),
          ". Se reportarán como NOT_EVALUABLE (no cambie el protocolo para ocultarlo).", call. = FALSE)
}

# ------------------------------------------------------------
# 3. Guardar productos
# ------------------------------------------------------------

metrics <- result$metrics |> dplyr::mutate(module = module_id, run_at = stamp())
write_csv_safe(metrics, file.path(DIR$method, "metrics_long.csv"))
register_table("tab_metrics_long", file.path(DIR$method, "metrics_long.csv"), "All performance and validation metrics produced by the analysis module", "06")

for (id in names(result$tables %||% list())) {
  tb <- result$tables[[id]]
  f <- file.path(DIR$tables, paste0(id, ".csv"))
  write_csv_safe(tibble::as_tibble(tb$data), f)
  register_table(id, f, tb$description %||% id, "06")
}

for (fg in result$figures %||% list()) {
  w <- fg$width_cm %||% 17; h <- fg$height_cm %||% 11
  paths <- save_fig(fg$plot, fg$id, w, h)
  register_figure(fg$id, paths, fg$caption_draft, "06", w, h)
}

facts <- c(result$facts, list(
  module = module_id, module_description = MODULE_INFO$description,
  protocol_version = frozen$version, protocol_hash = frozen$hash,
  runtime_minutes = runtime_min, run_at = stamp(),
  r_version = as.character(getRversion()),
  packages = paste(vapply(MODULE_INFO$packages, function(p) paste0(p, " ", utils::packageVersion(p)), character(1)), collapse = "; ")
))
yaml::write_yaml(list(facts = facts, limitations = result$limitations %||% character(), products = result$products %||% character()),
                 file.path(DIR$method, "method_facts.yaml"))
write_csv_safe(hyp_check, file.path(DIR$method, "hypothesis_evaluability.csv"))

write_md(file.path(DIR$evidence, "06_method_pack.md"), "Method as executed", list(
  "Method facts (generated from code)" = paste0("- **", names(flatten_named(facts)), "**: ", flatten_named(facts), collapse = "\n"),
  "Metrics" = metrics |> dplyr::select(subset, metric, value, ci_low, ci_high, unit, note),
  "Hypothesis evaluability" = hyp_check,
  "Design limitations stated by the module" = paste0("- ", result$limitations %||% "None stated", collapse = "\n"),
  "Registered tables" = if (file.exists(table_registry())) read_csv_safe(table_registry()) |> dplyr::select(id, description) else "None",
  "Registered figures" = if (file.exists(figure_registry())) read_csv_safe(figure_registry()) |> dplyr::select(id, caption_draft) else "None"
))

mark_stage("06_metodo", "PASSED", paste0(module_id, ", ", nrow(metrics), " métricas"))
cat("\nEtapa 06 completada.\n")
