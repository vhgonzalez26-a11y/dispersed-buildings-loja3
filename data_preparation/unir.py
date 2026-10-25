import os, glob, shutil, sys, warnings, csv
from collections import defaultdict
import pandas as pd, geopandas as gpd
warnings.filterwarnings('ignore')
from shapely.geometry import MultiPolygon, MultiLineString, MultiPoint, Polygon, LineString, Point

C = r'<PROJECT_ROOT>\Cartas'
OUT = r'<PROJECT_ROOT>\03_cartografia_igm_base'
DRY = '--run' not in sys.argv

sheets = {}
for d in sorted(x for x in os.listdir(C) if os.path.isdir(os.path.join(C, x))):
    cand = [r for r, _, f in os.walk(os.path.join(C, d)) if any(i.endswith('.shp') for i in f)]
    if cand:
        sheets[d.replace(' (1)', '')] = cand[0]

layers = defaultdict(dict)   # capa -> {carta: ruta}
for carta, p in sheets.items():
    for s in glob.glob(p + '/*.shp'):
        layers[os.path.basename(s)[:-4].lower()][carta] = s

def multi(g):
    if g is None or g.is_empty: return g
    t = g.geom_type
    if t == 'Polygon': return MultiPolygon([g])
    if t == 'LineString': return MultiLineString([g])
    return g

if not DRY:
    os.makedirs(OUT, exist_ok=True)
    gp = os.path.join(OUT, 'igm_50k_2013_base.gpkg')
    if os.path.exists(gp): os.remove(gp)

rows = []
for name in sorted(layers):
    parts = []
    for carta, s in layers[name].items():
        g = gpd.read_file(s)
        n0 = len(g)
        g = g[~g.geometry.isna() & ~g.geometry.is_empty].copy()
        if g.crs is None or g.crs.to_epsg() != 32717:
            g = g.set_crs(32717, allow_override=True) if g.crs is None else g.to_crs(32717)
        g['carta'] = carta.replace('_8417S', '')
        parts.append(g)
        rows.append([name, carta, n0, len(g)])
    allg = pd.concat(parts, ignore_index=True)
    allg = gpd.GeoDataFrame(allg, geometry='geometry', crs=32717)
    n_in = len(allg)
    if n_in == 0:
        continue
    allg['geometry'] = allg.geometry.apply(multi)
    # duplicados exactos entre cartas (franja de solape): misma geometria y mismos atributos
    attr = [c for c in allg.columns if c not in ('geometry', 'carta')]
    wk = allg.geometry.apply(lambda g: g.wkb_hex).tolist()
    at = allg[attr].values.tolist()
    key = pd.Series([w + '|' + '|'.join(map(str, a)) for w, a in zip(wk, at)], index=allg.index)
    dup = key.duplicated(keep='first')
    allg = allg[~dup].copy()
    rows.append([name, 'TOTAL', n_in, len(allg)])
    print(f'{name:34} cartas={len(parts)} entrada={n_in:6d} duplicados_borde={int(dup.sum()):4d} final={len(allg):6d}')
    if not DRY:
        allg.to_file(gp, layer=name, driver='GPKG')

if not DRY:
    # metadatos originales por carta
    for carta, p in sheets.items():
        dst = os.path.join(OUT, 'metadatos', carta.replace('_8417S', ''))
        os.makedirs(dst, exist_ok=True)
        for x in glob.glob(p + '/*.xml'):
            shutil.copy2(x, dst)
    with open(os.path.join(OUT, 'inventario_capas.csv'), 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.writer(fh); w.writerow(['capa', 'carta', 'elementos_origen', 'elementos_validos_o_finales']); w.writerows(rows)
    print('GPKG:', gp, round(os.path.getsize(gp) / 1e6, 1), 'MB')
