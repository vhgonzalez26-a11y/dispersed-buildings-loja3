import glob, warnings, numpy as np, pandas as pd, geopandas as gpd, pyogrio
from shapely import STRtree, get_num_geometries, get_parts
warnings.filterwarnings('ignore')
B = r'<SCRATCH_DIR>\gdb'
ROOT = r'<PROJECT_ROOT>'
G21 = glob.glob(B + r'\g2021\**\GEODATABASE_NACIONAL_2021.gdb', recursive=True)[0]
G14 = glob.glob(B + r'\g2014\**\GEODATABASE2014.gdb', recursive=True)[0]
aoi = gpd.read_file(ROOT + r'\05_aoi\aoi_estudio.gpkg', layer='aoi_amplio').set_index('centro')
bb = tuple(aoi.total_bounds)
IGM = ROOT + r'\03_cartografia_igm_base\igm_50k_2013_base.gpkg'
pd.set_option('display.width', 220)

v = pyogrio.read_dataframe(G21, layer='viv_p', bbox=bb, use_arrow=False)
man = pyogrio.read_dataframe(G21, layer='man_a', bbox=bb, use_arrow=False)
v['n_pts'] = get_num_geometries(v.geometry.values)
print('viv_p: features', len(v), '| puntos totales', int(v.n_pts.sum()), '| features con >1 punto:', int((v.n_pts > 1).sum()), '| max pts/feature', int(v.n_pts.max()))
# explotar a puntos
pts = []
ids = []
for i, g in zip(v.index, v.geometry.values):
    for p in get_parts(g):
        pts.append(p); ids.append(i)
P = gpd.GeoDataFrame({'fid': ids}, geometry=pts, crs=v.crs)
P['centro'] = ''
for c, r in aoi.iterrows():
    P.loc[P.within(r.geometry), 'centro'] = c
P = P[P.centro != '']
mu = man.geometry.union_all()
P['en_manzana'] = P.within(mu)
t = P.groupby('centro').agg(puntos=('fid', 'size'), en_manzana=('en_manzana', 'sum'))
t['dispersos'] = t.puntos - t.en_manzana
t['%_dispersos'] = (t.dispersos / t.puntos * 100).round(0)
man['centro'] = ''
for c, r in aoi.iterrows(): man.loc[man.geometry.representative_point().within(r.geometry), 'centro'] = c
mm = man[man.centro != ''].groupby('centro').agg(tot_edif=('tot_edif', 'sum'), tot_viv=('tot_viv', 'sum'))
t = t.join(mm)
# referencia IGM 2013
igm = pd.concat([gpd.read_file(IGM, layer='casa_p'), gpd.read_file(IGM, layer='edificio_p')])
igm['centro'] = ''
for c, r in aoi.iterrows(): igm.loc[igm.within(r.geometry), 'centro'] = c
t['casas_IGM_2013'] = igm[igm.centro != ''].groupby('centro').size()
t['pob_2022_AOI'] = pd.Series({'Loja': 213018, 'Catamayo': 28531, 'Malacatos': 6961, 'Vilcabamba': 7434, 'San_Pedro': 2071})
t['viv_ocup_2022_AOI'] = pd.Series({'Loja': 62575, 'Catamayo': 8236, 'Malacatos': 2210, 'Vilcabamba': 2375, 'San_Pedro': 634})
t['puntos/viv_ocup'] = (t.puntos / t.viv_ocup_2022_AOI).round(2)
print(t.T.to_string())
# donde estan los dispersos: distancia al nucleo IGM2013 y pendiente? solo conteo por pendiente rapida
print('\natributos extra de viv_p:', [c for c in v.columns if c != 'geometry'])
print('n_loc (localidades) principales:', v.n_loc.value_counts().head(6).to_dict())
P.drop(columns='geometry').assign(x=P.geometry.x, y=P.geometry.y).to_pickle(r'<SCRATCH_DIR>\viv2021_pts.pkl')
