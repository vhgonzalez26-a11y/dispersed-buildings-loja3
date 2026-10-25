# a1_correr_etapa.R
# Objetivo : ejecutar una etapa de la cadena (01 a 12) con la raiz del proyecto fijada, para poder lanzarla desde run_logged.py / ejecutar_R.py.
# Uso      : Rscript a1_correr_etapa.R <archivo de la etapa, p. ej. 01_tema_y_pregunta.R>
# Fecha    : 2026-10-22
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) >= 1L)
ROOT <- normalizePath(Sys.getenv("A1_ROOT", unset = getwd()), winslash = "/")
setwd(ROOT)
source(file.path(ROOT, args[1]), encoding = "UTF-8")
