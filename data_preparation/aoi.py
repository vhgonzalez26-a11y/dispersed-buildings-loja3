import os, glob, pickle, csv, warnings, time
import numpy as np, geopandas as gpd, pandas as pd, rasterio
from rasterio import features
from scipy import ndimage as ndi
from shapely import STRtree, points as mkpts
from shapely.geometry import shape, Point, LineString, MultiLineString, mapping
from shapely.ops import unary_union, linemerge
from shapely import wkb
warnings.filterwarnings('ignore')

ROOT = r'<PROJECT_ROOT>'
IGM = ROOT + r'\03_cartografia_igm_base\igm_50k_2013_base.gpkg'
OSM = ROOT + r'\04_osm_vias\osm_vias_snapshots.gpkg'
INS = ROOT + r'\05_aoi\_insumos'
OUT = ROOT + r'\05_aoi'
SPK = r'<SCRATCH_DIR>'
RES = 30.0

# ---------------- malla ----------------
with rasterio.open(INS + r'\slope_srtm30_x2.tif') as d:
    TR, H, W, CRS = d.transform, d.height, d.width, d.crs
    slope = d.read(1).astype(np.float32) / 2.0
with rasterio.open(INS + r'\worldcover2021_construido_pct_30m.tif') as d:
    wc = d.read(1).astype(np.float32) / 100.0     # fraccion construida por celda
def rast(geoms):
    geoms = [g for g in geoms if g is not None and not g.is_empty]
    if not geoms: return np.zeros((H, W), bool)
    return features.rasterize([(g, 1) for g in geoms], out_shape=(H, W), transform=TR, fill=0, all_touched=True, dtype='uint8').astype(bool)
dist = lambda m: ndi.distance_transform_edt(~m, sampling=RES).astype(np.float32) if m.any() else np.full((H, W), 1e9, np.float32)

# ---------------- datos ----------------
zu = gpd.read_file(IGM, layer='zona_urbana_a')
NUC = {
    'Loja': ['LOJA', 'TURUNUMA', 'LAS PITAS', 'SAN CAYETANO BAJO', 'CRISTO REY', 'AGUAS HEDIONDAS'],
    'Catamayo': ['CATAMAYO'],
    'Malacatos': ['MALACATOS'],
    'Vilcabamba': ['VILCABAMBA', 'SAN PEDRO DE VILCABAMBA'],
    'San_Pedro': ['SAN PEDRO DE LA BENDITA'],
}
CEN = list(NUC)
nuc_geom = {c: unary_union(zu[zu['nam'].isin(NUC[c])].geometry.tolist()) for c in CEN}
for c in CEN: assert not nuc_geom[c].is_empty, c
# AOI de la tesis (solo los tres satelites), para continuidad
tes = {}
base_t = ROOT + r'\02_datos_procesados\AOI_shp'
nm = {'Malacatos_shp': 'Malacatos', 'SP_shp_1': 'San_Pedro', 'Vilcabamba_shp': 'Vilcabamba'}
for f in glob.glob(base_t + '/*/*.shp'):
    tes[nm[os.path.basename(f)[:-4]]] = gpd.read_file(f).to_crs(32717).geometry.union_all()

pod = gpd.read_file(INS + r'\areas_protegidas_podocarpus_utm17s_wgs84.gpkg', layer='podocarpus').geometry.union_all()
pod = gpd.read_file(INS + r'\areas_protegidas_podocarpus_utm17s_wgs84.gpkg', layer='podocarpus').geometry.union_all()
igm_l = {n: gpd.read_file(IGM, layer=n) for n in ('via_l', 'sendero_l', 'rodera_l')}
casas = gpd.read_file(IGM, layer='casa_p')
edif = gpd.read_file(IGM, layer='edificio_p')
osm = {k: gpd.read_file(OSM, layer=k) for k in ('osm_2019_01', 'osm_2021_01', 'osm_2023_01', 'osm_2026_01')}
for k, v in osm.items(): print(k, len(v), round(v.length.sum() / 1000), 'km')

# ---------------- red combinada (IGM 2013 + OSM al 01-01-2026) ----------------
def pieces(gdf, step=50.0):
    mids, geoms, cls = [], [], []
    for g, c in zip(gdf.geometry, gdf['clase']):
        for ln in (g.geoms if g.geom_type == 'MultiLineString' else [g]):
            n = max(int(np.ceil(ln.length / step)), 1)
            for i in range(n):
                a, b = ln.interpolate(i / n, normalized=True), ln.interpolate((i + 1) / n, normalized=True)
                if a.equals(b): continue
                seg = LineString([a, b]); geoms.append(seg); mids.append(seg.interpolate(0.5, normalized=True)); cls.append(c)
    return np.array(mids, dtype=object), geoms, cls
def near(mids, lines_gdf):
    t = STRtree(lines_gdf.geometry.values)
    idx, d = t.query_nearest(mids, return_distance=True, all_matches=False)
    out = np.full(len(mids), np.inf); out[idx[0]] = d; return out

igm_all = pd.concat([g.assign(clase=n.replace('_l', '')) for n, g in igm_l.items()], ignore_index=True)
igm_all = gpd.GeoDataFrame(igm_all[['clase', 'geometry']], crs=32717)
osm_all = osm['osm_2026_01'].assign(clase=lambda d: d['highway'].fillna('otro'))[['clase', 'geometry']]
m_i, g_i, c_i = pieces(igm_all); m_o, g_o, c_o = pieces(osm_all)
d_i_osm = near(m_i, osm['osm_2026_01']); d_o_igm = near(m_o, igm_all)
TOL = 25.0
lab_i = np.where(d_i_osm <= TOL, 'igm_2013_confirmada_en_osm', 'solo_igm_2013')
lab_o = np.full(len(m_o), 'osm_2024_2025', object)
for key, name in (('osm_2023_01', 'osm_2022_2023'), ('osm_2021_01', 'osm_2020_2021'), ('osm_2019_01', 'osm_hasta_2019')):
    lab_o[near(m_o, osm[key]) <= TOL] = name
keep_o = d_o_igm > TOL                       # OSM que ya esta en el IGM no se duplica
rows = []
for lab, geoms, cls in ((lab_i, g_i, c_i), (lab_o[keep_o], np.array(g_o, dtype=object)[keep_o], np.array(c_o, dtype=object)[keep_o])):
    df = pd.DataFrame({'fuente': lab, 'clase': cls, 'geometry': list(geoms)})
    for (f, c), grp in df.groupby(['fuente', 'clase']):
        mg = unary_union(list(grp.geometry))
        if mg.geom_type == 'MultiLineString':
            mg = linemerge(mg)
        rows.append(dict(fuente=f, clase=c, longitud_km=round(sum(x.length for x in (mg.geoms if hasattr(mg, 'geoms') else [mg])) / 1000, 1), geometry=mg))
red = gpd.GeoDataFrame(rows, geometry='geometry', crs=32717)
gp = OUT + r'\aoi_estudio.gpkg'
if os.path.exists(gp): os.remove(gp)
red.to_file(gp, layer='red_vial_combinada', driver='GPKG')
res = red.groupby('fuente')['longitud_km'].sum().round(0)
print('\n=== COMPARACION OSM vs IGM (km) ===\n' + res.to_string())
tot_i = igm_all.length.sum() / 1000; tot_o = osm['osm_2026_01'].length.sum() / 1000
print(f'IGM 2013 total {tot_i:.0f} km | confirmado en OSM-2026: {(lab_i == "igm_2013_confirmada_en_osm").mean() * 100:.0f}% de los tramos')
print(f'OSM 2026 total {tot_o:.0f} km | ya presente en IGM: {(~keep_o).mean() * 100:.0f}% de los tramos | nuevo respecto al IGM: {keep_o.mean() * 100:.0f}%')

# ---------------- mascaras ----------------
net_m = rast(list(igm_all.geometry) + list(osm['osm_2026_01'].geometry))
house_pts = list(casas.geometry) + list(edif.geometry)
house_m = rast(house_pts)
d_net, d_house = dist(net_m), dist(house_m)
nuc_m = {c: rast([nuc_geom[c]]) for c in CEN}
allnuc = np.zeros((H, W), bool)
lab = np.zeros((H, W), np.int8)
for i, c in enumerate(CEN, 1):
    allnuc |= nuc_m[c]; lab[nuc_m[c]] = i
d_nuc, (ir, ic) = ndi.distance_transform_edt(~allnuc, sampling=RES, return_indices=True)
owner = lab[ir, ic]                        # centro del nucleo mas cercano (reparto)
tes_m = {c: rast([tes[c]]) & (owner == (CEN.index(c) + 1)) for c in tes}
pod_m = rast([pod])      # Parque Nacional Podocarpus: no se puede construir, se excluye del AOI
for i, c in enumerate(CEN, 1):
    assert not (nuc_m[c] & pod_m).any(), 'nucleo dentro de area protegida: ' + c
pod_m = rast([pod])      # Parque Nacional Podocarpus: no se puede construir, se excluye del AOI
for i, c in enumerate(CEN, 1):
    assert not (nuc_m[c] & pod_m).any(), 'nucleo dentro de area protegida: ' + c

SC = {
    'estricto': dict(pend=30, corr=250, casa=250, R=3000),
    'base':     dict(pend=35, corr=1000, casa=500, R=5000),
    'amplio':   dict(pend=40, corr=1500, casa=750, R=7000),
}
cell_ha = RES * RES / 1e4
hm = np.zeros((H, W), bool)
# casas dentro de la malla para medir captura
hp = np.array([(g.x, g.y) for g in house_pts])
rr, cc = rasterio.transform.rowcol(TR, hp[:, 0], hp[:, 1])
rr, cc = np.array(rr), np.array(cc); ok = (rr >= 0) & (rr < H) & (cc >= 0) & (cc < W)
rr, cc = rr[ok], cc[ok]
srows = []
for sc, p in SC.items():
    reach = (d_nuc <= p['R'])
    for c in tes_m: reach |= tes_m[c]
    access = (d_net <= p['corr']) | (d_house <= p['casa']) | allnuc
    A = reach & (slope <= p['pend']) & access
    A &= (owner > 0)
    A_sin_excl = A.copy()
    A &= ~pod_m
    recs = []
    for i, c in enumerate(CEN, 1):
        excl_ha = float((A_sin_excl & pod_m & (owner == i)).sum() * cell_ha)
        n_h_pod = int(pod_m[rr, cc][owner[rr, cc] == i].sum())
        Ac = A & (owner == i)
        # nucleo siempre incluido
        Ac |= nuc_m[c]
        shp = [(shape(g), v) for g, v in features.shapes(Ac.astype('uint8'), mask=Ac, transform=TR)]
        poly = unary_union([s for s, v in shp])
        poly = poly.buffer(15).buffer(-15).simplify(15)
        # areas
        ha = Ac.sum() * cell_ha
        cons = (wc * Ac).sum() * cell_ha; cons_tot = (wc * (owner == i) * (d_nuc <= 12000)).sum() * cell_ha
        capt_h = Ac[rr, cc].mean() * 100 if False else float((Ac[rr, cc]).sum())
        n_h_dom = float(((owner[rr, cc] == i) & (d_nuc[rr, cc] <= p['R'])).sum())
        recs.append(dict(centro=c, escenario=sc, area_ha=round(ha, 0), construido_wc2021_ha=round(cons, 1),
                         casas_igm_dentro=int(capt_h), casas_igm_en_alcance=int(n_h_dom),
                         excluido_podocarpus_ha=round(excl_ha), casas_igm_en_podocarpus=n_h_pod, geometry=poly))
        srows.append(dict(escenario=sc, centro=c, area_ha=round(ha), construido_wc2021_ha=round(cons, 1),
                          casas_dentro=int(capt_h), casas_alcance=int(n_h_dom),
                          excluido_podocarpus_ha=round(excl_ha), casas_igm_en_podocarpus=n_h_pod))
    gdf = gpd.GeoDataFrame(recs, geometry='geometry', crs=32717)
    gdf.to_file(gp, layer='aoi_' + sc, driver='GPKG')
    # construido de WC2021 a <= R+2km del nucleo del centro y que queda fuera
    tot_c = (wc * (d_nuc <= p['R'] + 2000) * (owner > 0)).sum() * cell_ha
    print(f'\n=== escenario {sc}: pend<={p["pend"]}, corredor {p["corr"]} m, casas {p["casa"]} m, alcance {p["R"] / 1000:.0f} km ===')
    print(gdf.drop(columns='geometry').to_string(index=False))
    print(f'AOI total {A.sum() * cell_ha:,.0f} ha | construido WC2021 dentro de AOI: {(wc * A).sum() * cell_ha:.0f} ha de {(wc * (d_nuc <= p["R"]) * (owner > 0)).sum() * cell_ha:.0f} ha en el alcance ({(wc * A).sum() / max((wc * (d_nuc <= p["R"]) * (owner > 0)).sum(), 1e-9) * 100:.1f}%)')
pd.DataFrame(srows).to_csv(OUT + r'\resumen_escenarios.csv', index=False, encoding='utf-8-sig')
pickle.dump({s: {c: None for c in CEN} for s in SC}, open(SPK + r'\aoi_done.pkl', 'wb'))
print('\nGPKG:', gp, round(os.path.getsize(gp) / 1e6, 1), 'MB')
