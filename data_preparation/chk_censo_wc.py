import warnings, numpy as np, pandas as pd, geopandas as gpd, rasterio
from scipy import ndimage as ndi
warnings.filterwarnings('ignore')
R = r'<PROJECT_ROOT>'
b = gpd.read_file(R + r'\06_censos\02_cartografia\censo_cartografia_aoi.gpkg', layer='edificios_2021_puntos', columns=['centro_aoi', 'en_manzana'])
b = b[b.centro_aoi != ''].reset_index(drop=True)
with rasterio.open(R + r'\05_aoi\_insumos\worldcover2021_construido_pct_30m.tif') as d:
    wc = d.read(1).astype(float) / 100; tr = d.transform; H, W = wc.shape
rr, cc = rasterio.transform.rowcol(tr, b.geometry.x.values, b.geometry.y.values)
rr, cc = np.array(rr), np.array(cc); ok = (rr >= 0) & (rr < H) & (cc >= 0) & (cc < W)
b = b[ok].reset_index(drop=True); rr, cc = rr[ok], cc[ok]
cnt = np.zeros((H, W)); np.add.at(cnt, (rr, cc), 1)
cnt3 = ndi.uniform_filter(cnt, 3) * 9   # edificios en la celda y sus 8 vecinas
b['wc'] = wc[rr, cc]; b['wc3'] = ndi.maximum_filter(wc, 3)[rr, cc]
print('=== Edificios censales 2021 frente a WorldCover 2021 (celdas de 30 m) ===')
for lab, m in (('en manzana', b.en_manzana == 1), ('dispersos', b.en_manzana == 0)):
    s = b[m]
    print(f'{lab:11} n={len(s):6d} | celda con algo de construido WC: {(s.wc > 0).mean() * 100:5.1f}% | con >=30%: {(s.wc >= .3).mean() * 100:5.1f}% | construido WC en la celda o sus vecinas: {(s.wc3 > 0).mean() * 100:5.1f}%')
print('\nPor centro, dispersos: % con construido WC en la celda o vecinas')
print(b[b.en_manzana == 0].groupby('centro_aoi').wc3.apply(lambda x: round((x > 0).mean() * 100, 1)).to_string())
# recall inverso: celdas WC construidas (>=30 %) con al menos un edificio censal a <=1 celda
cells = wc >= 0.3
near = ndi.maximum_filter(cnt, 3) > 0
aoi = gpd.read_file(R + r'\05_aoi\aoi_estudio.gpkg', layer='aoi_amplio')
from rasterio import features
am = features.rasterize([(g, 1) for g in aoi.geometry], out_shape=(H, W), transform=tr, fill=0, dtype='uint8').astype(bool)
c = cells & am
print(f'\nCeldas WC construidas (>=30 %) dentro del AOI amplio: {c.sum():,} | con un edificio censal en la celda o sus vecinas: {(c & near).sum() / c.sum() * 100:.1f}%')
print('Ejemplo: construido WC en celdas SIN edificio censal cercano (ha):', round(((wc * (c & ~near)).sum()) * 0.09, 1), 'de', round((wc * c).sum() * 0.09, 1))
