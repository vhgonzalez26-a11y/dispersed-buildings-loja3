import glob, warnings, numpy as np, pandas as pd, geopandas as gpd, pyogrio
from shapely import STRtree
warnings.filterwarnings('ignore')
B = r'<SCRATCH_DIR>'
ROOT = r'<PROJECT_ROOT>'
G14 = glob.glob(B + r'\gdb\g2014\**\GEODATABASE2014.gdb', recursive=True)[0]
IGM = ROOT + r'\03_cartografia_igm_base\igm_50k_2013_base.gpkg'
aoi = gpd.read_file(ROOT + r'\05_aoi\aoi_estudio.gpkg', layer='aoi_amplio').set_index('centro'); bb = tuple(aoi.total_bounds)
pd.set_option('display.width', 220)

def inc(gs):
    out = pd.Series('', index=gs.index)
    for c, r in aoi.iterrows(): out[gs.within(r.geometry)] = c
    return out

igm = pd.concat([gpd.read_file(IGM, layer='casa_p'), gpd.read_file(IGM, layer='edificio_p')], ignore_index=True)
igm['centro'] = inc(igm.geometry); igm = igm[igm.centro != ''].reset_index(drop=True)
c10 = pyogrio.read_dataframe(G14, layer='Ca04', bbox=bb, use_arrow=False)
c10['geometry'] = c10.geometry.representative_point(); c10['centro'] = inc(c10.geometry); c10 = c10[c10.centro != ''].reset_index(drop=True)
p21 = pd.read_pickle(B + r'\viv2021_pts.pkl')
c21 = gpd.GeoDataFrame(p21, geometry=gpd.points_from_xy(p21.x, p21.y), crs=32717)

def nn(a, b):
    """distancia de cada punto de a al punto mas cercano de b"""
    t = STRtree(b.geometry.values); idx, d = t.query_nearest(a.geometry.values, return_distance=True, all_matches=False)
    out = np.full(len(a), np.inf); out[idx[0]] = d; return out

rows = []
for c in aoi.index:
    ig = igm[igm.centro == c]; a10 = c10[c10.centro == c]; a21 = c21[c21.centro == c]
    disp = a21[~a21.en_manzana]; urb = a21[a21.en_manzana]
    r = dict(centro=c, igm_2013=len(ig), censo_2010_manzanas=len(a10), censo_2021_total=len(a21), censo_2021_en_manzana=len(urb), censo_2021_dispersos=len(disp))
    for tol in (30, 60):
        r[f'recall_IGM_vs_censo2010@{tol}m'] = round((nn(a10, ig) <= tol).mean() * 100) if len(a10) else np.nan
        r[f'recall_IGM_vs_2021_manz@{tol}m'] = round((nn(urb, ig) <= tol).mean() * 100) if len(urb) else np.nan
        r[f'recall_IGM_vs_2021_disp@{tol}m'] = round((nn(disp, ig) <= tol).mean() * 100) if len(disp) else np.nan
    r['precision_IGM(existe en censo2021)@60m'] = round((nn(ig, a21) <= 60).mean() * 100)
    r['precision_IGM(existe en censo2010)@60m'] = round((nn(ig, a10) <= 60).mean() * 100) if len(a10) else np.nan
    rows.append(r)
R = pd.DataFrame(rows).set_index('centro')
print(R.T.to_string())

# cota de omision en dispersos bajo escenarios de crecimiento 2013 -> 2021
print('\n=== Cota de captura del IGM en edificios DISPERSOS: IGM 2013 / (dispersos 2021 / (1+g)) ===')
rr = []
for c in ('Malacatos', 'Vilcabamba', 'San_Pedro'):
    disp = int(R.loc[c, 'censo_2021_dispersos']); igm_n = int(R.loc[c, 'igm_2013'])
    rr.append(dict(centro=c, disp_2021=disp, igm_2013=igm_n, **{f'g={int(g * 100)}%': f'{igm_n / (disp / (1 + g)) * 100:.0f}%' for g in (0, 0.25, 0.5, 1.0)}))
print(pd.DataFrame(rr).to_string(index=False))
