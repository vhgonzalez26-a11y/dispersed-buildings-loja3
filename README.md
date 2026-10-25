# Data processing and model code: location of dispersed buildings in the Andean periphery of Loja Province, Ecuador

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23643029.svg)](https://doi.org/10.5281/zenodo.23643029)

This repository contains the **data-processing and modelling code** and the **derived model outputs** of the study "Terrain, Roads and the Location of Dispersed Buildings in the Andean Periphery of Loja Province, Ecuador: A Spatially Cross-Validated Presence-Background Analysis".
Study design, hypotheses, interpretation and the manuscript are in the article; they are not part of this repository.
The unit of analysis is the building (census point), not the dwelling, and the model describes associations, not causal effects.

## Contents
| Folder or file | What it is |
|---|---|
| `04_inventario_datos.R`, `05_preprocesamiento_QC.R`, `06_ejecutar_metodo.R` | Data inventory (04), quality control and preprocessing (05), execution of the model (06) |
| `00_setup.R`, `R_helpers.R`, `run_stage.R` | Configuration, helper functions and a launcher (`Rscript run_stage.R 04_inventario_datos.R`) |
| `project_card.yaml` | Configuration: study area, data roles, model parameters, seed. It also holds the text that enters the hash of the frozen protocol (`outputs/protocol/`), which stage 06 checks before running |
| `methods/M03_presencia_fondo.R` | Presence-background inhomogeneous Poisson point-process model, spatial block cross-validation and block-bootstrap intervals (documentation in `methods/M03_presencia_fondo_README.md`) |
| `tests/test_M03.R` | Tests of the statistical functions on simulated data with known solutions |
| `tools/probar_M03.R` | Runs the model outside the stage gates (quick or full mode) |
| `data_preparation/` | Python scripts that built the study inputs from the original sources (see limits below) |
| `outputs/` | Inventory and quality-control tables (`04_*`, `05_*`), model tables (`tab_M03_*`), metrics, model figures, data provenance and quality-control decisions, run log |
| `DATA_SOURCES.md`, `data/` | Where to obtain each dataset; size and hash of the files used |

## Reproducing the processing and the model
1. Install R 4.4.3 and the packages in `requirements-R.txt` (versions in `outputs/logs/session_info.txt`); Python 3.12 with the pinned `requirements-python.txt` for the preparation scripts.
2. Obtain the data following `DATA_SOURCES.md` and place the files in `data/raw/`.
3. Set the project root (`A1_ROOT` environment variable, or run from this folder) and run `Rscript run_stage.R <script>` for `04_inventario_datos.R`, `05_preprocesamiento_QC.R` and `06_ejecutar_metodo.R`, in that order.
   The model stage took 43 minutes in the published run (R 4.4.3) with the fixed seed (2026). `outputs/logs/stage_status.csv` records the earlier steps as completed.
4. Compare your results with `outputs/tables/` and `outputs/method/metrics_long.csv`.

## Limits to reproducibility
- Third-party data are not redistributed (IGM, INEC and OpenStreetMap licences); they can be downloaded freely from the official sources listed in `DATA_SOURCES.md`.
- The Python scripts in `data_preparation/` were reconstructed from the working record and were not re-run as a single flow; some paths point to temporary folders. Use them as documentation of the procedure.
  The derived inputs (study area, clipped census layers, merged IGM base, OpenStreetMap snapshots) may therefore not be byte-identical to those of the published run; the hashes of the published inputs are in `data/MANIFIESTO_COPIA.csv`.
- The completeness of the census building mapping by distance to roads is not documented, so road-distance associations may be overstated; no independent coverage check was made.

## Licence and citation
Code: MIT licence (`LICENSE`). Derived results (tables, metrics and figures in `outputs/`) and documentation: Creative Commons Attribution 4.0 (`LICENSE_RESULTS.md`). The licences do not cover third-party data, which keep the terms listed in `DATA_SOURCES.md`.
Citation: see `CITATION.cff`. The archived version (v1.0.0) is on Zenodo: https://doi.org/10.5281/zenodo.23643029.

## Authorship
The software (data-processing workflow, model code and tests) was developed by Víctor González-Jaramillo (Departamento de Ingeniería Civil, Arquitectura y Geociencias, and Especialización en Gestión de Geoinformación con mención en Proyectos de Ingeniería, Universidad Técnica Particular de Loja, Loja 110101, Ecuador; vhgonzalez@utpl.edu.ec; https://orcid.org/0000-0002-2150-7690). The research reported in the article, and the contributions of all its authors, are described in the article.

## Acknowledgments
The author thanks Jhoana Paladines-Benitez (Departamento de Ciencias Empresariales, Universidad Técnica Particular de Loja) for her contributions to the methodology and to the curation of the data used in this study.

## Funding
This research was funded by Universidad Técnica Particular de Loja, grant number POA VIN-56.

## Use of AI
Claude Code (Anthropic; Claude Sonnet 5.5, accessed October 2026) was used as a programming assistant under the author's supervision; all code, parameters, model specifications and outputs were reviewed and validated by the author.
