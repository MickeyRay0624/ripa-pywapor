# RIPA pilot

## Wad Helal

`wad_helal.py` downloads the public FAO WaPOR v3 Level 3 AETI and NPP rasters and the IHE Delft boundaries. It writes a three-band GeoTIFF on the native 20 m grid, a map preview, one CSV record per wheat parcel, an edge coverage check and a summary with source URLs. The pilot follows the course dates from 1 October 2022 through 30 April 2023. [IHE Delft seasonal analysis notebook](https://github.com/wateraccounting/WaPORMOOC/blob/main/3_Python_for_WaPOR/N04_Seasonal_Aggregation_and_WP.ipynb)

```sh
python3.12 -m venv work/wad-helal/venv
work/wad-helal/venv/bin/python -m pip install -r pilot/requirements.txt
work/wad-helal/venv/bin/python pilot/wad_helal.py \
  --cache work/wad-helal/inputs --output work/wad-helal/results
```

The raster shows the command area wherever both datasets are present throughout the season. Field statistics are reported for wheat parcels with complete coverage. AETI is the estimated water consumed through evapotranspiration and canopy interception. Biomass is total plant dry matter estimated from NPP using a carbon fraction of 0.45. [FAO NPP methodology](https://github.com/un-fao/wapor-et-look/wiki/Net-Primary-Production)

## KGM

The AOI GeoJSON holds the command polygon and its enclosing rectangle. The top-level `bbox` is the input to `pywapor.Project` in west, south, east, north order. The standard `WaPOR3_level_3` configuration requires access to Copernicus Data Space, Copernicus Climate Data Store and NASA Earthdata for the VIIRS downloader. Configure these through pyWaPOR in the execution environment. [FAO pyWaPOR usage](https://github.com/un-fao/fao-wetl-pywapor)

```sh
work/wad-helal/venv/bin/python pilot/prepare_kgm.py \
  KGM_Irrigation_Command.shp work/kgm/staged/kgm_pywapor_extent.geojson
conda create -p work/kgm/model -c conda-forge python=3.12 pywapor=3.7.3
```

Pass inclusive reporting dates to `run_kgm.py` in the model environment. The dates below are for checking the input and model configuration. Use an interactive terminal and remove `--preflight` when starting a run for an agreed period with data access configured.

```sh
conda run --no-capture-output -p work/kgm/model python pilot/run_kgm.py \
  --aoi work/kgm/staged/kgm_pywapor_extent.geojson --start 2023-02-01 --end 2023-02-10 \
  --folder work/kgm/preflight --preflight
```

Each run caches its native NetCDF in a subfolder keyed by the model version, bounding box, configuration and model dates. The NetCDF retains evapotranspiration with interception, evaporation, transpiration, interception, reference evapotranspiration, root-zone moisture stress and NPP. The GeoTIFF export uses the AETI and NPP daily rates, weights each dekad by the number of days inside the reporting dates, and estimates dry biomass with a carbon fraction of 0.45. The December dekad uses its calendar midpoint for day of year.

The export uses a 20 × 20 m EPSG:32642 grid and the command polygon. Nearest-neighbour reprojection retains the model's missing-value pattern. The seasonal totals use pixels where AETI and NPP are valid in every included dekad.

- `kgm_dekads_20m.tif` has AETI, estimated dry biomass and biomass water productivity bands for each included dekad.
- `kgm_season_20m.tif` has the same three layers aggregated across the reporting dates.
- `kgm_edge_coverage.csv` counts valid and missing pixels within 20 m inside the boundary, by dekad and for the full reporting period. It uses the model grid footprint to identify any edge pixels outside the raster.
- `kgm_summary.json` records the reporting period, covered area, coverage and scheme-level results.

Use `--netcdf` with an existing native **dekadal** pyWaPOR ETLook output to rerun the export and coverage check. Choose `--output` to write the resulting files to a specific folder.

### Public-input KGM pilot for 1–10 February 2023

This run uses pyWaPOR 3.7.3 SE_ROOT v3 and ETLook v3 with staged Sentinel-2 reflectance, USGS Landsat Collection 2 surface temperature, ERA5 hourly weather, AgERA5 solar radiation, CHIRPS v2 rainfall and Copernicus DEM. It uses the FAO static model coefficients. Landsat temperature is screened at 3 K uncertainty. No thermal sharpening is applied. The input fetch scripts are specific to the KGM pilot and its dates. [pyWaPOR file inputs](https://www.fao.org/aquastat/py-wapor/notebooks/4_sideloading.html)

```sh
python3.12 -m venv work/kgm/era5-env
work/kgm/era5-env/bin/python -m pip install icechunk==2.3.1 zarr==3.4.2 xarray==2026.9.0 pcodec==1.1.0
work/kgm/era5-env/bin/python pilot/fetch_kgm_era5.py work/kgm/public-inputs/era5.npz

work/kgm/model/bin/python pilot/prepare_kgm_public_inputs.py \
  --aoi work/kgm/staged/kgm_pywapor_extent.geojson --folder work/kgm/public-inputs
work/kgm/model/bin/python pilot/run_kgm.py \
  --aoi work/kgm/staged/kgm_pywapor_extent.geojson --start 2023-02-01 --end 2023-02-10 \
  --folder work/kgm/public-runs --configuration-json work/kgm/public-inputs/kgm_public_inputs.json \
  --output work/kgm/public-results
```

### Run the staged pilot on an analysis server

On a Linux analysis VM, create the same model environment with conda or micromamba. Copy the `pilot` scripts, the AOI GeoJSON and the staged input folder with its source manifest. Run the two model-environment commands above, adding `--reuse-staged` to the first command to regenerate the configuration with paths on that server. The server runs SE_ROOT, ETLook and the command-area analysis. The raster, NetCDF and coverage table can be kept on the VM for Python or Jupyter analysis and archived in object storage for download. Use a fresh model run folder for the first server benchmark.

The 4-vCPU, 16-GiB Linux VM with a 100-GiB balanced disk discussed for this pilot costs about US$15.36 per month for 40 running hours or US$107.84 at 730 hours in Google Cloud Iowa. The disk is about US$10 per month in both examples. Regional Standard object storage is about US$0.02 per GiB per month. Network transfer, requests and other network resources are separate charges. A fresh local run took about three minutes from staged inputs through export and reached about 1.2 GiB in resident memory. A cloud benchmark is needed for server-side runtime and any estimate for multiple years. [VM pricing](https://cloud.google.com/products/compute/pricing/general-purpose) · [Disk pricing](https://cloud.google.com/products/block-storage) · [Object storage pricing](https://cloud.google.com/storage/pricing)
