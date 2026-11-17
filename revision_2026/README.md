# Analyses added in response to peer review (version 2.0.0)

These analyses were run in November 2026, after the original results of version 1.0.0 were known, and answer the comments of the reviewers of the article.
They use separate scripts that leave the frozen protocol and the original results unchanged; the original validation is reproduced exactly (`outputs/results/rev_reproduccion.txt`) before the new variants are run.

| Script (`scripts/`) | What it does | Main outputs (`outputs/results/`) |
|---|---|---|
| `a1_rev_00_inspeccion.R`, `a1_rev_00b_uso_edificios.R` | Inspection of the inputs and of the use labels of the buildings | `rev_sensibilidad_conteos.csv` |
| `a1_rev_01_validacion.R` | Validation variants: V0 (original), V1 (standardisation and knots from the training data only), B5 and B10 (V1 plus exclusion of training observations within 5 and 10 km of any test observation); fold maps and separation between training and test | `rev_cv_splits.csv`, `rev_cv_summary.csv`, `rev_loco.csv`, `rev_fold_blocks.csv`, `rev_fold_separation.csv`, `rev_reproduccion.txt` |
| `a1_rev_02_redes_viales.R` | Length, density, composition and coverage of the four road layers; vehicular OpenStreetMap layer | `rev_redes_composicion.csv`, `rev_redes_resumen.csv` |
| `a1_rev_03_efectos_clases_equivalentes.R` | Effects in fixed units (per 10 degrees of slope; per doubling of distance), equivalent road classes, study domains, paired differences, jackknife by block and bootstrap with B = 1000 | `rev_efectos_unidades_fijas.csv`, `rev_efectos_dominios.csv`, `rev_cv_clases_equivalentes.csv`, `rev_cv_pareada.csv`, `rev_jackknife_bloques.csv`, `rev_bootstrap_B1000.csv` |
| `a1_rev_04_sensibilidad_tipo_edificio.R` | Sensitivity to non-residential buildings and to the peri-urban fringe | `rev_sensibilidad_tipo_edificio.csv` |
| `a1_rev_05_figuras.R`, `a1_rev_06_tablas_revision.py` | Figures (validation folds; synthesis by centre) and tables built from the results above | `outputs/figures/`, `outputs/tables/` |

Reproduction: the scripts read the project outputs of version 1.0.0 (stage outputs and `data/processed/`) and need the third-party data listed in `DATA_SOURCES.md`.
Some scripts keep the absolute paths of the working machine at the top of the file; set them to your own project root before running. Tables with the prefix `tab_REV_` correspond to Tables 4 and 6 of the article and to Tables S17 to S23 of the Supplementary Materials.
