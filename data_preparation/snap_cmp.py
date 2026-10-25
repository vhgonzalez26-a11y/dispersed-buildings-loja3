import glob, warnings, numpy as np, geopandas as gpd, pyogrio
from shapely import force_2d
from shapely.geometry import box
from shapely.affinity import translate
from shapely.ops import unary_union
warnings.filterwarnings('ignore')
SP = r'<SCRATCH_DIR>'
ROOT = r'<PROJECT_ROOT>'
f = glob.glob(SP + r'\snap\**\*.shp', recursive=True)[0]
print(open(f[:-4] + '.prj').read()[:200])
n = gpd.read_file(f, encoding='latin1')
n['geometry'] = n.geometry.apply(force_2d)
print('CRS nuevo:', n.crs.to_string(), '| n:', len(n))
BB = box(661980, 9516000, 726000, 9580020)
aoi = gpd.read_file(ROOT + r'\05_aoi\aoi_estudio.gpkg', layer='aoi_amplio')
hit = n[n.intersects(BB.buffer(10000))]
print('\nAreas SNAP oficiales a <=10 km de la malla de estudio:')
for _, r in hit.iterrows():
    print(f'  {r.nam} ({r.map}) csnap={r.csnap} | area oficial {r.are:.0f} ha | dentro de la malla: {r.geometry.intersection(BB).area / 1e4:.0f} ha')

nm = n[n.nam.str.upper().str.contains('PODOCARPUS')]
print('\nPodocarpus nuevo: n poligonos =', len(nm), '| area ha =', round(nm.geometry.area.sum() / 1e4), '| campos:', nm[['nam', 'map', 'ror', 'rom', 'edel', 'esc']].to_dict('records'))
new = nm.geometry.union_all()

# antiguo: reproyectado (PSAD56->WGS84) y sin transformar (coordenadas PSAD56 leidas como WGS84)
old_raw = gpd.read_file(ROOT + r'\Areas_protegidas\areas_p.shp')
old_raw = old_raw[old_raw.NOMBRE.fillna('').str.contains('Podocarpus')]
old_tr = old_raw.to_crs(32717).geometry.union_all()
old_nt = old_raw.set_crs(32717, allow_override=True).geometry.union_all()

def m(a, b, lab):
    inter = a.intersection(b).area; uni = a.union(b).area
    print(f'  {lab:34} area {a.area / 1e4:8.0f} ha | IoU con nuevo {inter / uni * 100:5.1f}% | fuera del nuevo {a.difference(b).area / 1e4:6.0f} ha | nuevo no cubierto {b.difference(a).area / 1e4:6.0f} ha | Hausdorff {a.hausdorff_distance(b):6.0f} m')
print('\nComparacion de geometrias:')
m(old_tr, new, 'antiguo CON transformacion de datum')
m(old_nt, new, 'antiguo SIN transformacion (UTM igual)')
for nm_, g in (('antiguo transformado', old_tr), ('antiguo sin transformar', old_nt)):
    cx, cy = g.centroid.x - new.centroid.x, g.centroid.y - new.centroid.y
    print(f'  desfase de centroides {nm_:24}: dE={cx:7.0f} m dN={cy:7.0f} m')

# desfase por traslacion optima (maximiza IoU) del antiguo respecto al nuevo
best = (0, 0, 0)
for dx in range(-600, 601, 50):
    for dy in range(-600, 601, 50):
        t = translate(old_nt, dx, dy)
        iou = t.intersection(new).area / t.union(new).area
        if iou > best[0]: best = (iou, dx, dy)
print(f'  traslacion que mejor alinea el antiguo SIN transformar con el nuevo: dE={best[1]} m dN={best[2]} m (IoU {best[0] * 100:.1f}%)')
best = (0, 0, 0)
for dx in range(-600, 601, 50):
    for dy in range(-600, 601, 50):
        t = translate(old_tr, dx, dy)
        iou = t.intersection(new).area / t.union(new).area
        if iou > best[0]: best = (iou, dx, dy)
print(f'  traslacion residual del antiguo TRANSFORMADO respecto al nuevo: dE={best[1]} m dN={best[2]} m (IoU {best[0] * 100:.1f}%)')

# dentro de la zona de estudio
for lab, g in (('antiguo transformado', old_tr), ('NUEVO oficial', new)):
    for _, c in aoi.iterrows():
        a = g.intersection(c.geometry).area / 1e4
        if a > 0: print(f'  {lab:22} {c.centro:11} dentro del AOI amplio actual: {a:7.0f} ha')
unary_union([new]) if False else None
gpd.GeoDataFrame({'nam': ['Podocarpus'], 'fuente': ['SNAP v_fa210_snap_d'], 'geometry': [new]}, crs=32717).to_file(SP + r'\podocarpus_nuevo.gpkg', driver='GPKG')
