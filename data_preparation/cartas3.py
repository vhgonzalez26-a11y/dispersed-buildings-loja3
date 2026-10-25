import os, glob, re, pickle, warnings
import numpy as np, geopandas as gpd
warnings.filterwarnings('ignore')
from shapely import wkb
from shapely.geometry import box, Point
from shapely.ops import unary_union

C = r'<PROJECT_ROOT>\Cartas'
SP = r'<SCRATCH_DIR>'
sheets = {}
for d in sorted(x for x in os.listdir(C) if os.path.isdir(os.path.join(C, x))):
    cand = [r for r, _, f in os.walk(os.path.join(C, d)) if any(i.endswith('.shp') for i in f)]
    if not cand:
        continue
    p = cand[0]
    lay, bnds, crs = {}, [], None
    for s in sorted(glob.glob(p + '/*.shp')):
        n = os.path.basename(s)[:-4]
        try:
            g = gpd.read_file(s)
        except Exception:
            lay[n] = ('ERR', 0); continue
        crs = crs or g.crs
        lay[n] = (len(g), g.geom_type.iloc[0] if len(g) else '-')
        if len(g):
            bnds.append(g.total_bounds)
    b = np.array(bnds)
    ext = (b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max())
    dates = []
    for x in glob.glob(p + '/*.xml'):
        t = open(x, encoding='utf-8', errors='ignore').read()
        dates += re.findall(r'<(?:pubdate|caldate|begdate|enddate)[^>]*>(\d{8})', t)
    sheets[d] = dict(ext=ext, crs=crs, p=p, lay=lay, dates=sorted(set(dates)))
    print(f'\n===== {d} | {crs.to_string() if crs else None} | {len(lay)} capas')
    print('extension UTM17S: E %d-%d  N %d-%d  (%.1f x %.1f km)' % (ext[0], ext[2], ext[1], ext[3], (ext[2]-ext[0])/1000, (ext[3]-ext[1])/1000))
    key = [k for k in lay if re.search(r'via|sendero|casa|edificio|poblad|curva', k)]
    print({k: lay[k][0] for k in key})
    print('fechas en metadatos (min/max):', (sheets[d]['dates'][0], sheets[d]['dates'][-1]) if dates else None)
    v = gpd.read_file(os.path.join(p, 'via_l.shp'))
    print('via_l campos:', [c for c in v.columns if c != 'geometry'])
    for c in v.columns:
        if v[c].dtype == object and v[c].nunique() < 12 and c != 'geometry':
            print('   ', c, v[c].value_counts().to_dict())

# ---- cobertura de las zonas ----
Z = {k: wkb.loads(bytes.fromhex(v)) for k, v in pickle.load(open(SP + r'\zonas.pkl', 'rb')).items()}
pt = lambda lon, lat: gpd.GeoSeries([Point(lon, lat)], crs=4326).to_crs(32717).iloc[0]
targets = {}
for p_ in ['Malacatos', 'San_Pedro', 'Vilcabamba']:
    targets[f'{p_} nucleo'] = Z[(p_, 'nucleo')]
    targets[f'{p_} nucleo+3km'] = unary_union([Z[(p_, 'nucleo')], Z[(p_, 'anillo_0_2km')], Z[(p_, 'anillo_2_3km')]])
targets['Loja ciudad (r=6 km)'] = pt(-79.214, -3.986).buffer(6000)
targets['Catamayo (r=4 km)'] = pt(-79.373, -3.985).buffer(4000)
targets['Control Loja (canton-satelites)'] = Z[('Loja', 'control_canton_sin_satelites')]
boxes = {d: box(*s['ext']) for d, s in sheets.items()}
allsh = unary_union(list(boxes.values()))
print('\n===== COBERTURA de las zonas de estudio por las 4 cartas =====')
for n, g in targets.items():
    cov = g.intersection(allsh).area / g.area * 100
    per = {d.split('_8417S')[0][:14]: round(g.intersection(b).area / g.area * 100) for d, b in boxes.items() if g.intersects(b)}
    print(f'{n:34} cubierto {cov:5.1f}%  | por carta: {per}')
