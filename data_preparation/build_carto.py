import glob, os, warnings, numpy as np, pandas as pd, geopandas as gpd, pyogrio, shapely
from shapely.ops import unary_union
warnings.filterwarnings('ignore')
S = r'<SCRATCH_DIR>'
ROOT = r'<PROJECT_ROOT>'
G14 = glob.glob(S + r'\gdb\g2014\**\GEODATABASE2014.gdb', recursive=True)[0]
G21 = glob.glob(S + r'\gdb\g2021\**\GEODATABASE_NACIONAL_2021.gdb', recursive=True)[0]
OUT = ROOT + r'\06_censos\02_cartografia\censo_cartografia_aoi.gpkg'
aoi = gpd.read_file(ROOT + r'\05_aoi\aoi_estudio.gpkg', layer='aoi_amplio').set_index('centro')
AU = aoi.geometry.union_all()
REG = AU.buffer(5000)                      # area de estudio + margen de 5 km
bb = REG.bounds
if os.path.exists(OUT): os.remove(OUT)

def centro_de(geoms):
    pts = geoms if geoms.geom_type.iloc[0] == 'Point' else geoms.representative_point()
    out = pd.Series('', index=geoms.index)
    for c, r in aoi.iterrows(): out[pts.within(r.geometry)] = c
    return out

def clip_write(df, name, ren=None):
    if ren: df = df.rename(columns=ren)
    pts = df.geometry if df.geom_type.iloc[0] in ('Point', 'MultiPoint') else df.geometry
    sel = df.geometry.intersects(REG)
    df = df[sel].copy()
    df['centro_aoi'] = centro_de(df.geometry)
    df.to_file(OUT, layer=name, driver='GPKG')
    n_in = int((df.centro_aoi != '').sum())
    print(f'  {name:34} {len(df):>8} elementos ({n_in} dentro del AOI amplio)', flush=True)
    return df

def rd(g, layer, cols=None):
    return pyogrio.read_dataframe(g, layer=layer, bbox=bb, columns=cols, use_arrow=False)

print('--- cartografia censo 2010 (ciclo 2010, ajustada 2014 con ortofoto) ---')
c10 = rd(G14, 'Ca04'); c10['TOTAL_VIV'] = pd.to_numeric(c10.TOTAL_VIV.astype(str).str.strip(), errors='coerce')
clip_write(c10, 'edificios_2010', {'DPA_MAN': 'manzana', 'TOTAL_VIV': 'total_viv', 'TIPO_EDIF': 'tipo_edif', 'NOM_EDIF': 'nom_edif', 'FUENTE': 'fuente'})
clip_write(rd(G14, 'GEO_MAN2014'), 'manzanas_2010')
clip_write(rd(G14, 'GEO_SEC2014'), 'sectores_2010')
clip_write(rd(G14, 'GEO_EJES2014'), 'ejes_2010')
clip_write(rd(G14, 'GEO_LOC2014'), 'localidades_2010')
print('--- cartografia precenso 2021 (capa etiquetada 2020 / CPV2020) ---')
v = rd(G21, 'viv_p'); v['geometry'] = v.geometry.apply(lambda g: g.geoms[0] if g.geom_type == 'MultiPoint' else g)
v = gpd.GeoDataFrame(v, geometry='geometry', crs=32717)
mz = rd(G21, 'man_a'); mu = mz.geometry.union_all()
v['en_manzana'] = v.geometry.within(mu).astype(int)
clip_write(v, 'edificios_2021_puntos')
clip_write(rd(G21, 'ca04_a'), 'edificios_2021_poligonos')
clip_write(mz, 'manzanas_2021')
clip_write(rd(G21, 'ejes_l'), 'ejes_2021')
clip_write(rd(G21, 'loc_p'), 'localidades_2021')
print('--- sectores censo 2022 con indicadores ---')
sec = pyogrio.read_dataframe(S + r'\censos\CapaSectores\sectores_anonimizados.gpkg', bbox=bb, use_arrow=False)
sec = sec.rename(columns={'pob_t': 'poblacion_2022', 'v_pres': 'viv_ocupadas_2022', 'p_hog': 'tam_hogar_prom', 'serv_b': 'viv_serv_basicos', 'pob_serv_b': 'pob_serv_basicos'})
clip_write(sec, 'sectores_2022_indicadores')
print('--- limites de parroquias INEC 2001 (provincia Loja) ---')
par = gpd.read_file(S + r'\censos\SHP\SHP\nxparroquias.shp'); par = par[par.DPA_PROVIN == '11']
par.to_file(OUT, layer='limites_parroquias_2001', driver='GPKG'); print('  limites_parroquias_2001', len(par))

# ---------- cobertura de los edificios censales por las AOI ----------
print('\n=== Cobertura de los edificios censales 2021 por el AOI amplio (area de estudio +5 km) ===')
pts = gpd.read_file(OUT, layer='edificios_2021_puntos')
pod = gpd.read_file(ROOT + r'\05_aoi\_insumos\areas_protegidas_snap_podocarpus_oficial.gpkg').geometry.union_all()
pts['en_podocarpus'] = pts.within(pod)
pts['dist_aoi'] = shapely.distance(pts.geometry.values, AU)
den = pts[~pts.en_podocarpus]
tot = len(den)
print(f'edificios en la region (sin Podocarpus): {tot}')
for lab, m in (('dentro del AOI amplio', den.dist_aoi == 0), ('a 0-1 km fuera del AOI', (den.dist_aoi > 0) & (den.dist_aoi <= 1000)),
               ('a 1-3 km fuera', (den.dist_aoi > 1000) & (den.dist_aoi <= 3000)), ('a 3-5 km fuera', (den.dist_aoi > 3000))):
    print(f'  {lab:28} {int(m.sum()):>7} ({m.mean() * 100:5.1f}%)')
print(f'  dentro de Podocarpus (excluidos): {int(pts.en_podocarpus.sum())} edificios')
disp = den[den.en_manzana == 0]
print(f'edificios DISPERSOS (fuera de manzana) sin Podocarpus: {len(disp)} | dentro del AOI: {(disp.dist_aoi == 0).mean() * 100:.1f}% | a <=1 km fuera: {((disp.dist_aoi > 0) & (disp.dist_aoi <= 1000)).mean() * 100:.1f}% | a >1 km: {(disp.dist_aoi > 1000).mean() * 100:.1f}%')
print('\nEdificios fuera del AOI pero cerca (<=1 km): donde estan (localidad mas frecuente)')
near = den[(den.dist_aoi > 0) & (den.dist_aoi <= 1000)]
print(near.n_loc.str.strip().replace('', 'sin localidad').value_counts().head(6).to_dict())
