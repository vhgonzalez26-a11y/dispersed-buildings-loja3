import os, time, json, requests, warnings, itertools
import geopandas as gpd
from shapely.geometry import LineString
warnings.filterwarnings('ignore')
H = {'User-Agent': 'cambio-ciudades-loja-research/1.0 (academic study of road network history)'}
EP = "https://overpass-api.de/api/interpreter"
OUT = r'<PROJECT_ROOT>\04_osm_vias'
CACHE = os.path.join(OUT, '_cache'); os.makedirs(CACHE, exist_ok=True)
S, W, N, E = -4.30, -79.51, -3.88, -79.00
NX, NY = 3, 2
xs = [W + (E - W) * i / NX for i in range(NX + 1)]
ys = [S + (N - S) * j / NY for j in range(NY + 1)]
tiles = [(ys[j], xs[i], ys[j + 1], xs[i + 1]) for j in range(NY) for i in range(NX)]
# franjas extra: AOI amplio de Vilcabamba (sur) y San Pedro (oeste y norte) superaban el bbox original
tiles += [(-4.35, -79.525, -4.30, -79.00),    # franja sur
          (-3.88, -79.525, -3.85, -79.00),    # franja norte
          (-4.30, -79.525, -3.88, -79.51)]    # franja oeste
# franjas extra: AOI amplio de Vilcabamba (sur) y San Pedro (oeste y norte) superaban el bbox original
tiles += [(-4.35, -79.525, -4.30, -79.00),    # franja sur
          (-3.88, -79.525, -3.85, -79.00),    # franja norte
          (-4.30, -79.525, -3.88, -79.51)]    # franja oeste
DATES = [('2019-01-01T00:00:00Z', 'osm_2019_01'), ('2021-01-01T00:00:00Z', 'osm_2021_01'),
         ('2023-01-01T00:00:00Z', 'osm_2023_01'), ('2026-01-01T00:00:00Z', 'osm_2026_01')]


def fetch(body, tag):
    fn = os.path.join(CACHE, tag + '.json')
    if os.path.exists(fn):
        return json.load(open(fn, encoding='utf-8'))
    for k in range(12):
        try:
            r = requests.post(EP, data={'data': body}, headers=H, timeout=(20, 400))
            if r.status_code == 200 and r.text.lstrip().startswith('{'):
                js = r.json()
                json.dump(js, open(fn, 'w', encoding='utf-8'))
                return js
            msg = f'{r.status_code}'
        except Exception as e:
            msg = type(e).__name__
        wait = min(30 * (k + 1), 240)
        print(f'   {tag}: intento {k + 1} fallo ({msg}); espero {wait}s', flush=True)
        time.sleep(wait)
    raise RuntimeError('sin respuesta ' + tag)


gp = os.path.join(OUT, 'osm_vias_snapshots.gpkg')
if os.path.exists(gp): os.remove(gp)
for date, lay in DATES:
    t0 = time.time(); rows = {}
    for n, (s, w, nn, e) in enumerate(tiles):
        d = f'[date:"{date}"]' if date else ''
        body = f'[out:json][timeout:300]{d};way["highway"]({s:.4f},{w:.4f},{nn:.4f},{e:.4f});out geom tags;'
        js = fetch(body, f'{lay}_{n}')
        for el in js['elements']:
            if el.get('type') == 'way' and len(el.get('geometry', [])) >= 2:
                tg = el.get('tags', {})
                rows[el['id']] = dict(osm_id=el['id'], highway=tg.get('highway'), name=tg.get('name'), surface=tg.get('surface'),
                                      geometry=LineString([(p['lon'], p['lat']) for p in el['geometry']]))
        print(f'{lay} mosaico {n + 1}/{len(tiles)} ok ({len(rows)} vias acumuladas)', flush=True)
    g = gpd.GeoDataFrame(list(rows.values()), geometry='geometry', crs=4326).to_crs(32717)
    g.to_file(gp, layer=lay, driver='GPKG')
    print(f'== {lay}: {len(g)} vias, {g.length.sum() / 1000:.0f} km, {time.time() - t0:.0f}s', flush=True)
print('LISTO', gp, flush=True)
