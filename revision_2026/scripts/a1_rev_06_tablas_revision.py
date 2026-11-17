"""
a1_rev_06_tablas_revision.py
Objetivo : convertir los resultados de los analisis de la revision de Land (outputs/revision/rev_*.csv) en las tablas con formato de publicacion que usa el plan 08
           (outputs/tables/tab_REV_*.csv, columnas ya redactadas en ingles y valores ya redondeados). Tablas principales: tab_REV_fixed_units (Tabla 5) y tab_REV_validation_variants (Tabla 7);
           suplementarias: composicion de redes (S17), validacion completa por variante (S18), separacion y tamanos (S19), clases equivalentes (S20), diferencias pareadas (S21),
           sensibilidad por tipo de edificio (S22) y estabilidad del bootstrap (S23).
Entradas : 000_A1_research/outputs/revision/rev_*.csv; outputs/tables/tab_M03_domains.csv
Salidas  : 000_A1_research/outputs/tables/tab_REV_*.csv
Uso      : python run_logged.py --no-isolate "<proposito>" 07_infraestructura_y_registro/a1_rev_06_tablas_revision.py
Fecha    : 2026-11-07
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
RAIZ = Path(r"C:\Trabajos_y_Publicaciones\02_cambio_ciudades\000_A1_research")
REV = RAIZ / "outputs" / "revision"
TAB = RAIZ / "outputs" / "tables"
f3 = lambda x: f"{float(x):.3f}"
f2 = lambda x: f"{float(x):.2f}"
f1 = lambda x: f"{float(x):.1f}"
ci = lambda r, lo, hi: f"{f3(r)} ({f3(lo)}–{f3(hi)})"
rg = lambda m, lo, hi: f"{f3(m)} ({f3(lo)}–{f3(hi)})"
CENTRO = {"Catamayo": "Catamayo", "Loja": "Loja", "Malacatos": "Malacatos", "San_Pedro": "San Pedro de la Bendita", "Vilcabamba": "Vilcabamba"}
VARIANTE = {"V0": ("Original (whole-domain standardisation)", "None"), "V1": ("Training-only standardisation and knots", "None"),
            "B5": ("Training-only + buffer", "5"), "B10": ("Training-only + buffer", "10")}


def guardar(df, nombre):
    out = TAB / nombre
    df.to_csv(out, index=False, encoding="utf-8")
    print("guardado", out.name, df.shape)


# ---------------- Tabla 5 (principal): efectos en unidades fijas ----------------
ef = pd.read_csv(REV / "rev_efectos_unidades_fijas.csv")
dom = pd.read_csv(REV / "rev_efectos_dominios.csv")
ETIQ = {"igm_vias": "IGM roads without tracks (main set)", "igm_all": "IGM roads with tracks and paths", "osm2019": "OpenStreetMap 2019, all classes", "census2020": "Census road axes",
        "major_igm": "Primary and secondary roads: IGM", "major_osm": "Primary and secondary roads: OpenStreetMap", "vehicular_osm": "Vehicular roads: OpenStreetMap"}
filas = []
for fs, etq in ETIQ.items():
    d = ef[ef.feature_set == fs].set_index("term")
    filas.append(["A. Road layer (main domain; block-bootstrap interval)", etq, "17,088",
                  ci(d.loc["slope_10deg", "ratio"], d.loc["slope_10deg", "boot_lo"], d.loc["slope_10deg", "boot_hi"]),
                  ci(d.loc["road_doubling", "ratio"], d.loc["road_doubling", "boot_lo"], d.loc["road_doubling", "boot_hi"]),
                  ci(d.loc["core_doubling", "ratio"], d.loc["core_doubling", "boot_lo"], d.loc["core_doubling", "boot_hi"])])
DOM = {"extended": "Main domain (wide AOI + 5 km)", "aoi_amplio": "Wide AOI", "aoi_base": "Base AOI", "aoi_estricto": "Strict AOI", "core_buffer": "10 km core buffer (geometric; no slope or road rules)"}
for k, etq in DOM.items():
    d = dom[dom.domain == k].set_index("term")
    filas.append(["B. Study domain (IGM roads without tracks; model-based interval)", etq, f"{int(d['n_presences'].iloc[0]):,}",
                  ci(d.loc["slope_10deg", "ratio"], d.loc["slope_10deg", "model_lo"], d.loc["slope_10deg", "model_hi"]),
                  ci(d.loc["road_doubling", "ratio"], d.loc["road_doubling", "model_lo"], d.loc["road_doubling", "model_hi"]),
                  ci(d.loc["core_doubling", "ratio"], d.loc["core_doubling", "model_lo"], d.loc["core_doubling", "model_hi"])])
guardar(pd.DataFrame(filas, columns=["Panel", "Road layer or domain", "Buildings (n)", "Slope (per 10°)", "Distance to roads (per doubling)", "Distance to urban core (per doubling)"]), "tab_REV_fixed_units.csv")

# ---------------- Tabla 7 (principal) y S18 (completa): variantes de validacion ----------------
cv = pd.read_csv(REV / "rev_cv_splits.csv")
lo = pd.read_csv(REV / "rev_loco.csv")
filas, completa = [], []
for v, (prep, buf) in VARIANTE.items():
    c = cv[cv.variant == v]
    ig = c[c.model == "spl_igm_vias"].sort_values(["repeat_id", "fold"])
    os_ = c[c.model == "spl_osm2019"].sort_values(["repeat_id", "fold"])
    d = os_["cauc_pds"].values - ig["cauc_pds"].values
    l = lo[(lo.variant == v) & (lo.model == "spl_igm_vias")]
    imin = l.loc[l.cauc_pds.idxmin()]
    filas.append([prep, buf, f"{ig['n_train_p'].mean():,.0f}", rg(ig.cauc_pds.mean(), ig.cauc_pds.min(), ig.cauc_pds.max()), rg(os_.cauc_pds.mean(), os_.cauc_pds.min(), os_.cauc_pds.max()),
                  f"{d.mean():.3f} ({int((d > 0).sum())} of {len(d)})", f"{f3(l.cauc_pds.min())}–{f3(l.cauc_pds.max())} (lowest: {CENTRO[imin['centre']]})"])
    for m in ["igm_vias", "spl_igm_vias", "osm2019", "spl_osm2019"]:
        cm = c[c.model == m]
        lm = lo[(lo.variant == v) & (lo.model == m)]
        completa.append([prep, buf, m, "Spatial block (15 folds)", f"{cm.n_train_p.mean():,.0f}", rg(cm.auc.mean(), cm.auc.min(), cm.auc.max()), rg(cm.cauc_pds.mean(), cm.cauc_pds.min(), cm.cauc_pds.max()),
                         rg(cm.auc_null.mean(), cm.auc_null.min(), cm.auc_null.max()), rg(cm.boyce.mean(), cm.boyce.min(), cm.boyce.max())])
        completa.append([prep, buf, m, "Leave one centre out (5 centres)", f"{lm.n_train_p.mean():,.0f}", rg(lm.auc.mean(), lm.auc.min(), lm.auc.max()), rg(lm.cauc_pds.mean(), lm.cauc_pds.min(), lm.cauc_pds.max()),
                         rg(lm.auc_null.mean(), lm.auc_null.min(), lm.auc_null.max()), rg(lm.boyce.mean(), lm.boyce.min(), lm.boyce.max())])
guardar(pd.DataFrame(filas, columns=["Preprocessing", "Buffer (km)", "Training presences, spatial CV (mean)", "cAUC_pds, IGM spline model: mean (range)", "cAUC_pds, OSM spline model: mean (range)",
                                     "OSM minus IGM difference: mean (folds higher)", "Leave-one-centre-out, IGM spline model: range"]), "tab_REV_validation_variants.csv")
guardar(pd.DataFrame(completa, columns=["Preprocessing", "Buffer (km)", "Model", "Validation", "Training presences (mean)", "AUC: mean (range)", "cAUC_pds: mean (range)", "AUC_null: mean (range)", "Boyce: mean (range)"]), "tab_REV_validation_full.csv")

# ---------------- S19: separacion entre prueba y entrenamiento y tamanos ----------------
sep = pd.read_csv(REV / "rev_fold_separation.csv")
filas = []
for v in VARIANTE:
    s = sep[(sep.variant == v) & (sep.scheme == "spatial_block")]
    filas.append([VARIANTE[v][0], VARIANTE[v][1], "Spatial block (15 folds, mean and range)", f"{s.n_train_p.mean():,.0f}", f"{s.n_excl_buffer_p.mean():,.0f}", f"{s.n_test_p.mean():,.0f}", f"{s.n_pds_pairs.mean():,.0f}",
                  f"{100 * s.pds_share_retained.mean():.0f}", f"{s.dist_min_m.mean() / 1000:.2f} ({s.dist_min_m.min() / 1000:.2f}–{s.dist_min_m.max() / 1000:.2f})", f"{s.dist_p05_m.mean() / 1000:.2f} ({s.dist_p05_m.min() / 1000:.2f}–{s.dist_p05_m.max() / 1000:.2f})", f"{s.dist_median_m.mean() / 1000:.2f}"])
    for c in CENTRO:
        s = sep[(sep.variant == v) & (sep.scheme == "leave_one_centre_out") & (sep.unit == c)]
        if len(s):
            r = s.iloc[0]
            filas.append([VARIANTE[v][0], VARIANTE[v][1], f"Leave one centre out: {CENTRO[c]}", f"{int(r.n_train_p):,}", f"{int(r.n_excl_buffer_p):,}", f"{int(r.n_test_p):,}", f"{int(r.n_pds_pairs):,}", f"{100 * r.pds_share_retained:.0f}",
                          f"{r.dist_min_m / 1000:.2f}", f"{r.dist_p05_m / 1000:.2f}", f"{r.dist_median_m / 1000:.2f}"])
guardar(pd.DataFrame(filas, columns=["Preprocessing", "Buffer (km)", "Validation", "Training presences (n)", "Presences removed by the buffer (n)", "Test presences (n)", "Pairs retained by pairwise distance sampling (n)",
                                     "Test presences retained (%)", "Minimum distance from test to nearest training presence (km)", "5th percentile (km)", "Median (km)"]), "tab_REV_fold_separation.csv")

# ---------------- S17: composicion de redes ----------------
co = pd.read_csv(REV / "rev_redes_composicion.csv")
res = pd.read_csv(REV / "rev_redes_resumen.csv")
NOMBRE = {"igm_vias": "IGM roads without tracks", "igm_all": "IGM roads with tracks and paths", "osm2019": "OpenStreetMap, 1 January 2019", "census2020": "INEC road axes, 2021"}
CLASE = {"major": "Primary and secondary", "local": "Local", "track_path": "Tracks and paths", "no_hierarchy": "No hierarchy"}
filas = []
for n in NOMBRE:
    for _, r in co[co.network == n].sort_values(["clase", "detalle"]).iterrows():
        filas.append([NOMBRE[n], CLASE[r.clase], r.detalle, f"{r.length_km:,.1f}", f"{100 * r.share_of_length:.1f}", f"{r.density_km_per_km2:.3f}"])
    t = res[res.network == n].iloc[0]
    filas.append([NOMBRE[n], "Total", f"{int(t.n_cells * t.coverage_1km_cells):,} of {int(t.n_cells):,} 1 km cells with a road ({100 * t.coverage_1km_cells:.1f}%)", f"{t.length_km:,.1f}", "100.0", f"{t.density_km_per_km2:.3f}"])
guardar(pd.DataFrame(filas, columns=["Road layer", "Class group", "Class or coverage", "Length (km)", "Share of layer length (%)", "Density (km km⁻²)"]), "tab_REV_road_composition.csv")

# ---------------- S20/S21: clases equivalentes ----------------
cl = pd.read_csv(REV / "rev_cv_clases_equivalentes.csv")
par = pd.read_csv(REV / "rev_cv_pareada.csv")
ORD = [("igm_vias", "Original layers", "IGM roads without tracks"), ("igm_all", "Original layers", "IGM roads with tracks and paths"), ("osm2019", "Original layers", "OpenStreetMap 2019, all classes"),
       ("census2020", "Original layers", "Census road axes"), ("major_igm", "Primary and secondary roads", "IGM primary and secondary routes"), ("major_osm", "Primary and secondary roads", "OpenStreetMap trunk, primary and secondary"),
       ("vehicular_igm", "Vehicular roads", "IGM roads without tracks"), ("vehicular_osm", "Vehicular roads", "OpenStreetMap without tracks, paths and steps"),
       ("all_igm", "All lines", "IGM roads with tracks and paths"), ("all_osm", "All lines", "OpenStreetMap, all classes")]
filas = []
for fs, grupo, etq in ORD:
    d = ef[ef.feature_set == fs].set_index("term")
    c = cl[(cl.feature_set == fs) & (cl.spline == False)]
    cs = cl[(cl.feature_set == fs) & (cl.spline == True)]
    filas.append([grupo, etq, ci(d.loc["slope_10deg", "ratio"], d.loc["slope_10deg", "boot_lo"], d.loc["slope_10deg", "boot_hi"]), ci(d.loc["road_doubling", "ratio"], d.loc["road_doubling", "boot_lo"], d.loc["road_doubling", "boot_hi"]),
                  ci(d.loc["core_doubling", "ratio"], d.loc["core_doubling", "boot_lo"], d.loc["core_doubling", "boot_hi"]), rg(c.auc_pds.mean(), c.auc_pds.min(), c.auc_pds.max()), rg(cs.auc_pds.mean(), cs.auc_pds.min(), cs.auc_pds.max())])
guardar(pd.DataFrame(filas, columns=["Class group", "Road layer", "Slope (per 10°)", "Distance to roads (per doubling)", "Distance to urban core (per doubling)", "AUC_pds, log-linear model: mean (range)", "AUC_pds, spline model: mean (range)"]), "tab_REV_equivalent_classes.csv")
NOM_CL = {"mayor": "Primary and secondary roads", "vehicular": "Vehicular roads", "todo": "All lines"}
filas = [[NOM_CL[r["class"]], "Spline" if r.spline else "Log-linear", f3(r.mean_igm), f3(r.mean_osm), f3(r.mean_diff_osm_minus_igm), f"{f3(r.diff_min)} to {f3(r.diff_max)}", f"{int(r.n_folds_osm_higher)} of {int(r.n_folds)}"] for _, r in par.iterrows()]
guardar(pd.DataFrame(filas, columns=["Road class group", "Model", "AUC_pds, IGM (mean)", "AUC_pds, OpenStreetMap (mean)", "Difference, OpenStreetMap minus IGM (mean)", "Range of paired differences", "Folds with OpenStreetMap higher"]), "tab_REV_equivalent_paired.csv")

# ---------------- S22: sensibilidad por tipo de edificio ----------------
se = pd.read_csv(REV / "rev_sensibilidad_tipo_edificio.csv")
cn = pd.read_csv(REV / "rev_sensibilidad_conteos.csv")
NOM_S = {"S0_todos": "All dispersed buildings (reference)", "S1_sin_no_residencial": "Without buildings with a non-residential label", "S2_sin_zonas_urbanas": "Without buildings and background inside urban cores or the official urban limit",
         "S3a_mas_de_1km_del_nucleo": "Only beyond 1 km from the urban cores", "S3b_mas_de_2km_del_nucleo": "Only beyond 2 km from the urban cores"}
filas = []
for k, etq in NOM_S.items():
    d = se[se.subset == k].set_index("term")
    n = cn[cn.subset == k].iloc[0]
    filas.append([etq, f"{int(n.n_presences):,}", f"{int(n.n_background):,}", ci(d.loc["slope_10deg", "ratio"], d.loc["slope_10deg", "boot_lo"], d.loc["slope_10deg", "boot_hi"]),
                  ci(d.loc["road_doubling", "ratio"], d.loc["road_doubling", "boot_lo"], d.loc["road_doubling", "boot_hi"]), ci(d.loc["core_doubling", "ratio"], d.loc["core_doubling", "boot_lo"], d.loc["core_doubling", "boot_hi"]),
                  rg(d["auc_pds_mean"].iloc[0], d["auc_pds_min"].iloc[0], d["auc_pds_max"].iloc[0])])
guardar(pd.DataFrame(filas, columns=["Subset", "Buildings (n)", "Background points (n)", "Slope (per 10°)", "Distance to IGM roads (per doubling)", "Distance to urban core (per doubling)", "AUC_pds, spatial CV: mean (range)"]), "tab_REV_building_type.csv")

# ---------------- S23: estabilidad del bootstrap ----------------
jk = pd.read_csv(REV / "rev_jackknife_bloques.csv")
bb = pd.read_csv(REV / "rev_bootstrap_B1000.csv").set_index("term")
efm = ef[ef.feature_set == "igm_vias"].set_index("term")
filas = []
for t, etq in [("slope_10deg", "Slope (per 10°)"), ("road_doubling", "Distance to IGM roads (per doubling)"), ("core_doubling", "Distance to urban core (per doubling)")]:
    big = jk.loc[jk.n_presences_dropped.idxmax()]
    filas.append([etq, f3(efm.loc[t, "ratio"]), f"{f3(efm.loc[t, 'boot_lo'])}–{f3(efm.loc[t, 'boot_hi'])}", f"{f3(bb.loc[t, 'lo'])}–{f3(bb.loc[t, 'hi'])}", f"{f3(jk[t].min())}–{f3(jk[t].max())}", f3(big[t])])
guardar(pd.DataFrame(filas, columns=["Term", "Estimate", "Block bootstrap, B = 300: 95% interval", "Block bootstrap, B = 1,000: 95% interval", "Leaving out one block at a time: range of estimates (29 blocks)", "Leaving out the largest block (4,283 buildings)"]), "tab_REV_bootstrap_stability.csv")
print("listo")
