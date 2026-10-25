# Module M03 — Presence-background as a Poisson point process

File: `methods/M03_presencia_fondo.R` · Parameters: `project_card.yaml → method_params.M03_presencia_fondo` · Tests: `tests/test_M03.R` · Run on real data outside the stage gates: `tools/probar_M03.R`

## 1. Question the module answers

What explains *where* dispersed buildings (buildings outside census blocks) are located in the Andean periphery of Loja and Catamayo, relative to the environment that is available? The module estimates the **relative intensity** of dispersed buildings as a function of slope, distance to roads and distance to the urban core (and, as a candidate, distance to rivers). It does not estimate the probability that a parcel is built on, nor the processes that produced the buildings.

## 2. Model

Dispersed buildings are treated as a point process with intensity λ(s) = exp(β₀ + Σ βₖ zₖ(s)), where zₖ are standardised covariates. The parameters are estimated by maximising the point-process likelihood:

ℓ(β) = Σᵢ xᵢ'β − c Σⱼ exp(xⱼ'β),  i = presences, j = uniform background points, c = domain area / number of background points (km²).

- The likelihood is concave; it is solved by Newton-Raphson with line search (function `ipp_fit`, base R, no dependencies).
- In the limit it is **equivalent** to logistic regression with infinitely weighted background (Fithian and Hastie 2013) and to quadrature-weighted Poisson regression (Warton and Shepherd 2010; Renner et al. 2015).
- exp(βₖ) is an **intensity ratio per standard deviation** (`intensity_ratio_<variable>`): how much the expected density of buildings per km² changes when the variable increases by one SD, holding the others fixed.
- Model-based standard errors (inverse Hessian) assume conditional independence. Because buildings cluster, the intervals of the main effects are obtained with a **spatial block bootstrap** (blocks resampled with replacement).

**Why not the 3:1 logistic regression of the exploratory analysis.** With a misspecified model, unweighted logistic regression does not coincide with the point process unless the background is very large (Fithian and Hastie 2013); the background needed is of the order of 10⁵ points (Warton and Shepherd 2010: about 86,000; Renner et al. 2015: closer to 100,000 than to the 10,000 of MaxEnt). The module keeps the 3:1 logistic regression only as a comparison (`sens_logit_3to1`, table `tab_M03_naive_vs_ipp`).

## 3. Input data (all in `data/raw`; hashes in `data/MANIFIESTO_COPIA.csv`)

| Role | File · layer · field | Notes and licence |
|---|---|---|
| Presences | `censo_cartografia_aoi.gpkg` · `edificios_2021_puntos` · `en_manzana == 0` | 2021 census buildings (INEC). A building is not a dwelling. See `DATA_SOURCES.md` |
| Background exclusion | same source · `manzanas_2021`; `areas_protegidas_snap_podocarpus_oficial.gpkg` | The background is sampled outside census blocks and outside Podocarpus National Park |
| Main domain | `aoi_estudio.gpkg` · `aoi_amplio` (sensitivity: `aoi_estricto`, `aoi_base`) | The AOI was built with slope, road-corridor and IGM-house rules: see threat A |
| Extended domain | wide AOI + 5 km (coverage of the census points) and a 10 km radius around urban cores | Reduces dependence on the AOI. Limited to the extent of the slope raster |
| Slope | `slope_srtm30_x2.tif` (SRTM 30 m; value = 2 × degrees; 0 = slope < 0.25°) | Public domain |
| Roads | IGM 1:50,000 (`via_l`; with tracks and unpaved roads: `sendero_l`, `rodera_l`), OSM 2019 (`osm_2019_01`), census axes (`ejes_2021`) | IGM: do not redistribute, cite. OSM: ODbL |
| Urban core | IGM `zona_urbana_a` (11 names in `core_names`) | Distance 0 inside the polygon |
| Rivers | IGM `rio_l` | Includes streams; candidate amenity (H7) |
| 2021 urban limit | `limite_urbano_2021.shp` | Loja, Malacatos and Vilcabamba only; restricted use under permission from the GAD of Loja |

## 4. Parameters and their justification

| Parameter (card) | Value | Justification |
|---|---|---|
| `background.ratio_main` | 20 background points per presence (341,760 in the main domain) | Order of magnitude for convergence (Warton and Shepherd 2010; Renner et al. 2015); checked with the ladder |
| `background.ladder` | 1, 3, 10, 20, 40 | Convergence: relative change < 2% in each intensity ratio and change in log-likelihood within 2 Monte Carlo standard errors between the two largest sizes (adapted from the criterion of Renner et al. 2015) |
| `background.seeds` | 2026, 7, 11 | Sensitivity to background randomness (Barbet-Massin et al. 2012; Whitford et al. 2024) |
| Background | uniform over the domain | Not weighted for sampling bias: census points are a census, not opportunistic records (Whitford et al. 2024 advise against spatial weighting when there is no sampling bias) |
| `cv.block_size_m` | `auto`: range of residual autocorrelation (Moran's I of Pearson residuals in 1-km cells), rounded to km, between 2 and 16 km | Blocks must be at least as large as the dependence range (Roberts et al. 2017; Valavi et al. 2019). Sensitivity: 2, 4, 8 and 16 km |
| `cv.folds` × `repeats` | 5 × 3 | Stability of the assignment of blocks to folds |
| Null model | distance to the training set | Hijmans (2012) |
| `bootstrap.B` | 300 | Block resampling; 2.5–97.5 percentiles |
| `nonlinear_df` | 3 | Natural splines (main model) to allow non-linear effects on log-intensity; the linear model summarises the average effect |
| Standardisation | mean and SD of the domain background; distances as log(1 + d) | Intensity ratios per SD comparable across variables; log because distances are right-skewed |
| Collinearity | Pearson r, VIF | Rule of thumb |r| > 0.7 (Dormann et al. 2013) |
| `thresholds.road_distance_m` | 250 | Descriptive value; **not an official threshold** |
| `thresholds.road_distance_sensitivity_m` | `[25, 30]` | Regulatory sensitivity thresholds: 25 m from the axis (Regulation of the Roads Law, Art. 4) and 30 m with the additional 5 m setback. They produce `share_within_25m_state_road` and `share_within_30m_state_road` only for proxies of the state road network (`state_networks`: IGM primary route; IGM primary + secondary; OSM trunk/primary; OSM + secondary), also reported in `tab_M03_shares`. The 25 m right-of-way applies to the state network; there is no official road inventory in the project and census axes have no hierarchy |
| `thresholds.slope_deg` | `null` | There is no official planning threshold; only a descriptive reference table is produced |

## 5. Procedure (sections of the code)

1. **Data**: main and sensitivity domains = layer − 2021 blocks − Podocarpus, clipped to the raster extent. Presences = points with `en_manzana == 0` inside the extended domain. Background: one uniform bag in the extended domain (bags of the smaller domains are nested by membership) and additional bags for the seeds. Exact covariates per vector layer (distance to the nearest feature; `terra::distance(pairwise = TRUE)`, identical to `sf::st_distance` in a test with 20,000 points). A cache `data/processed/M03_dataset_<key>.rds` is stored; its key depends on the data parameters.
2. **Collinearity**: VIF per feature set and correlations among all predictors.
3. **Fit** per feature set and **background convergence** (ladder).
4. **Residual autocorrelation and block size**: 1-km cells, Pearson residual (O − E)/√E, Moran correlogram with 99 permutations; range = up to the last band with significantly positive I.
5. **Spatial cross-validation** (blocks; folds assigned at random to blocks) and **random** cross-validation (reference to quantify optimism; Roberts et al. 2017; Pohjankukka et al. 2017). Per fold: AUC; AUC after pairwise distance sampling; null AUC (−distance to training); **calibrated AUC, cAUC = AUC + 0.5 − max(0.5, null AUC)** (Hijmans 2012); continuous Boyce index (Hirzel et al. 2006) on the rank scale of the test background (a variant that avoids unequal windows with long-tailed intensities).
6. **Effects**: intensity ratios with block bootstrap.
7. **Heterogeneity**: leave-one-centre-out (training without the blocks that touch the test centre), effects by centre and by distance ring from the core (0–2, 2–5, 5–10 km).
8. **Sensitivities**: background size, 3:1 logistic regression, seeds, domain (four alternatives), block size, splines.
9. **Descriptive shares**: percentage of buildings within 250 m of each road source versus percentage of area; percentage of buildings above reference slope values; percentage inside the official 2021 urban limit.
10. **Outputs** (next section).

## 6. Outputs

- **Metrics** (`metrics_long.csv`): `intensity_ratio_<slope|dist_road_igm_vias|dist_road_igm_all|dist_road_osm2019|dist_road_census2020|dist_nucleo|dist_river>`, `auc_spatial_cv`, `auc_random_cv`, `auc_pds_spatial_cv`, `auc_null_distance`, `cauc_spatial_cv`, `boyce_spatial_cv`, `log_likelihood`, `n_presences`, `n_background`, `block_size_m`, `residual_range_m`, `building_density_per_km2`, `share_within_250m_road`, `share_within_25m_state_road`, `share_within_30m_state_road`, `share_inside_urban_limit`. Subsets: name of the feature set, `center_<name>`, `center_min`, `ring_<a>_<b>km`, `sens_*`.
- **Tables** (`outputs/tables/`): `tab_M03_coefficients`, `_cv_summary`, `_cv_folds`, `_loco`, `_background_convergence`, `_loglik_by_background`, `_naive_vs_ipp`, `_collinearity_vif`, `_collinearity_pairs`, `_residual_correlogram`, `_domains`, `_block_size_sensitivity`, `_rings`, `_ring_effects`, `_centre_effects`, `_shares`, `_slope_reference`, `_urban_limit`.
- **Figures** (`outputs/figures/`): effects, validation (random vs spatial; AUC, cAUC, Boyce), background convergence, correlogram, Boyce predicted-to-expected curve, validation by centre, domain and presences.
- **Method facts** (`method_facts.yaml`): generated by the code; they document the parameters actually used.

## 7. ODMAP structure (Zurell et al. 2020) → where each element is

| ODMAP section | Elements | Where |
|---|---|---|
| Overview | Objective (explanatory inference, not prediction), model type (Poisson point process), scale (points; domain of about 756 km²), study area, hypotheses | Project card, sections 1–2 of this document |
| Data | Presence data, background, predictors, resolution, sources, biases | Section 3; `tab_M03_domains`; threats below |
| Model | Algorithm, variable selection (sets fixed a priori), collinearity, fitting, non-linearity | Sections 2, 4 and 5; `tab_M03_collinearity_*` |
| Assessment | Spatial partition, block size, metrics, null model, uncertainty (bootstrap), sensitivities | Sections 5–6; `tab_M03_cv_*`, `tab_M03_loco` |
| Prediction | This module does not produce an operational predictive map | — |

## 8. Verification of the implementation

`tests/test_M03.R` (18 tests with simulated data of known solution; result in `tests/test_M03_resultado.csv`): recovery of point-process coefficients (|error| < 3 SE), stability with a smaller background, AUC = Mann-Whitney, pairwise distance sampling, Boyce index (informative > 0.8; random ≈ 0 in 60 repetitions), cAUC formula, Wilson interval, Moran's I (autocorrelation detected; white noise not). All 18 tests pass.

## 9. Threats to validity and how the module treats them

| Threat | Treatment | Residual |
|---|---|---|
| A. Circularity of the AOI (built with slope, road corridor and IGM houses) | Four alternative domains, including AOI + 5 km and a 10 km radius around cores | A domain fully independent of the covariates would require extracting more points from the national INEC cartography |
| B. Inflated AUC in presence-background designs | cAUC, AUC after pairwise distance sampling, Boyce index, effects | — |
| C. Endogeneity of roads | IGM roads without tracks as the most exogenous source; comparison IGM/OSM/census axes; interpretation as association | An instrument (e.g. historical maps, Duranton and Turner 2012) is outside this module |
| D. Temporal mismatch (roads 2009–2013, buildings 2021) | Reported | Cannot be resolved with these data |
| E. Positional error of census points | Not assessed | — |
| F. Building ≠ dwelling | Stated | 2022 seasonal-dwelling data exist by sector, not by building |
| G. Spatial autocorrelation | Blocks, block bootstrap, correlogram | Model-based SEs are reported only as a reference |
| H. No independent validation | Leave-one-centre-out; domains; seeds | Another building inventory is needed |
| I. Sensitivity to the background | Ladder, seeds, comparison with 3:1 logistic regression | — |
| J. Linearity | Natural splines (main model) | — |

## 10. Methodological references

Warton and Shepherd 2010 (Ann. Appl. Stat.); Fithian and Hastie 2013 (Ann. Appl. Stat.); Renner et al. 2015 (Methods Ecol. Evol.); Barbet-Massin et al. 2012; Whitford et al. 2024 (Ecol. Model.); Phillips et al. 2009; Roberts et al. 2017; Valavi et al. 2019; Pohjankukka et al. 2017; Hijmans 2012; Segurado et al. 2006; Lobo et al. 2008; Hirzel et al. 2006; Dormann et al. 2013; Zurell et al. 2020; Mets et al. 2017; Sun and Robinson 2018; Chomitz and Gray 1996; Duranton and Turner 2012.

## How to run

See `README.md` at the root of the repository (stages 04 to 06, quick test with `tools/probar_M03.R`, module tests with `tests/test_M03.R`).
