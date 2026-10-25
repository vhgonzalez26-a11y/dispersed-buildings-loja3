import glob, os, warnings, numpy as np, pandas as pd, geopandas as gpd, pyogrio
from shapely import STRtree
warnings.filterwarnings('ignore')
B = r'<SCRATCH_DIR>\gdb'
ROOT = r'<PROJECT_ROOT>'
G14 = glob.glob(B + r'\g2014\**\GEODATABASE2014.gdb', recursive=True)[0]
G21 = glob.glob(B + r'\g2021\**\GEODATABASE_NACIONAL_2021.gdb', recursive=True)[0]
aoi = gpd.read_file(ROOT + r'\05_aoi\aoi_estudio.gpkg', layer='aoi_amplio').set_index('centro')
bb = tuple(aoi.total_bounds)
pd.set_option('display.width', 230); pd.set_option('display.max_columns', 30)

def load(gdb, layer, cols=None):
    g = pyogrio.read_dataframe(gdb, layer=layer, bbox=bb, columns=cols, use_arrow=False)
    return g

def assign(g):
    """asigna cada elemento al centro cuyo AOI contiene su punto representativo"""
    pts = g.geometry.representative_point() if not g.geom_type.iloc[0].endswith('Point') else g.geometry
    out = pd.Series('', index=g.index)
    for c, r in aoi.iterrows():
        m = pts.within(r.geometry); out[m] = c
    return out

res = {}
print('=== Lectura (bbox de las AOI) ===')
L = {
 '2010 Ca04 (edificios, polig.)': (G14, 'Ca04', None),
 '2010 GEO_EDIF_P (edificios importantes, pts)': (G14, 'GEO_EDIF_P', None),
 '2010 GEO_MAN2014 (manzanas)': (G14, 'GEO_MAN2014', ['DPA_MANZAN', 'DPA_ANIO', 'FUENTE']),
 '2010 GEO_EJES2014 (ejes viales)': (G14, 'GEO_EJES2014', ['DPA_PARROQ', 'FUENTE']),
 '2021 viv_p (viviendas/edificios, pts)': (G21, 'viv_p', None),
 '2021 ca04_a (edificios, polig.)': (G21, 'ca04_a', None),
 '2021 man_a (manzanas con tot_edif/tot_viv)': (G21, 'man_a', None),
 '2021 ejes_l (ejes viales)': (G21, 'ejes_l', ['parroquia', 'tipo_eje', 'anio', 'fuente']),
}
data = {}
for k, (g, l, c) in L.items():
    d = load(g, l, c); data[k] = d
    print(f'{k:48} {len(d):>8} elementos en el bbox')
rows = []
for k, d in data.items():
    if len(d) == 0: continue
    d = d.copy(); d['centro'] = assign(d)
    d = d[d.centro != '']
    if 'EJES' in k or 'ejes' in k:
        s = d.groupby('centro').apply(lambda x: x.length.sum() / 1000).round(0)
        rows.append(pd.DataFrame({k: s}))
    else:
        rows.append(pd.DataFrame({k: d.groupby('centro').size()}))
    data[k] = d
T = pd.concat(rows, axis=1).fillna(0).astype(int)
print('\n=== Elementos dentro de cada AOI amplio (ejes en km) ===')
print(T.T.to_string())

# viviendas / edificios totales
c04 = data['2010 Ca04 (edificios, polig.)']; v21 = data['2021 man_a (manzanas con tot_edif/tot_viv)']
print('\n=== Ca04 2010: atributos ===')
print(c04.drop(columns='geometry').head(4).to_string())
print('TIPO_EDIF:', c04.TIPO_EDIF.value_counts().head(8).to_dict(), '| FUENTE:', c04.FUENTE.value_counts().head(4).to_dict())
print('TOTAL_VIV por centro (2010, suma sobre Ca04):'); print(c04.groupby('centro').TOTAL_VIV.agg(['size', 'sum']).to_string())
print('\n=== man_a 2021: atributos y totales por centro ===')
print(v21.drop(columns='geometry').head(3).to_string()); print('anio:', v21.anio.value_counts().to_dict(), '| fuente:', v21.fuente.value_counts().head(4).to_dict())
print(v21.groupby('centro').agg(manzanas=('tot_viv', 'size'), tot_edif=('tot_edif', 'sum'), tot_viv=('tot_viv', 'sum')).to_string())
p21 = data['2021 viv_p (viviendas/edificios, pts)']
print('\n=== viv_p 2021: atributos ===')
print(p21.drop(columns='geometry').head(4).to_string()); print('edif_censo:', p21.edif_censo.value_counts().head(5).to_dict() if 'edif_censo' in p21 else '')
print('cod_otros:', p21.cod_otros.value_counts().head(6).to_dict())
a21 = data['2021 ca04_a (edificios, polig.)']
print('\n=== ca04_a 2021: atributos ===')
print(a21.drop(columns='geometry').head(3).to_string()); print('edif_censo:', a21.edif_censo.value_counts().head(5).to_dict())
pd.to_pickle({k: v.drop(columns='geometry', errors='ignore') for k, v in data.items()}, r'<SCRATCH_DIR>\gdb_attrs.pkl')
