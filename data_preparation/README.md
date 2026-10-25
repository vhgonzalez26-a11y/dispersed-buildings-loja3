# Legacy preparation scripts (Python): building the study inputs

Status: RECONSTRUCTED on 20 Oct 2026 from the record of the working session (the writes and edits of each file were re-applied). They were not re-run as a single flow. Several paths point to a temporary session folder (`<SCRATCH_DIR>`) that no longer exists and must be adjusted before running. The scripts that had edits that could not be re-applied (`aoi.py`, `unir.py`) may differ from the version that produced the data; use them as documentation of the procedure, not as a guarantee of exact reproduction.

The results of these steps are the input files listed in `DATA_SOURCES.md` (hashes in `data/MANIFIESTO_COPIA.csv`).

| Script | What it does |
|---|---|
| `cartas3.py` | Checks the coverage of the IGM sheets over the study area |
| `unir.py` | Merges the IGM 1:50,000 sheets into a single base (`03_cartografia_igm_base`) |
| `osm_dl2.py` | Downloads the OSM road network (Overpass, four snapshots) in tiles with retries |
| `snap_cmp.py` | Compares the old protected-areas layer (PSAD56) with the official SNAP layer |
| `aoi.py` | Builds the AOI with 30 m grid rules (three scenarios) |
| `gdb_aoi.py`, `gdb_aoi2.py` | Extract buildings, blocks and axes from the INEC geodatabases inside the AOI |
| `build_carto.py` | Assembles `censo_cartografia_aoi.gpkg` with the 2010, 2021 and 2022 layers |
| `chk_censo_wc.py` | Compares census buildings with ESA WorldCover |
| `valid_igm.py` | Compares 2013 IGM houses with census buildings (IGM capture) |

Note on inputs: `aoi.py` also reads polygons of the three satellite centres (Malacatos, San Pedro de la Bendita and Vilcabamba) delineated in earlier work by the author, and a 30 m built-up fraction raster derived from ESA WorldCover 2021 that is used only to measure how much built-up area each scenario captures; neither input is redistributed.

Pending: port to R or consolidate into a reproducible flow what is needed for the processing stages (04 and 05).
