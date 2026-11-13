"""Prepare public observations for the February 2023 KGM pyWaPOR pilot."""

import argparse
import copy
import json
import math
from concurrent.futures import ThreadPoolExecutor
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import requests
import rasterio
import rioxarray
import xarray as xr
from osgeo import gdal
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.warp import Resampling
from rasterio.windows import Window, from_bounds

import pywapor
from pywapor.collect.product.SENTINEL2 import apply_qa, calc_normalized_difference, calc_r0
from pywapor.enhancers.wind import adjust_wind_height


START, END = date(2023, 2, 1), date(2023, 2, 10)
SOURCE_START, SOURCE_END = date(2023, 1, 26), date(2023, 2, 17)
OPTICAL_START = date(2023, 1, 23)
CRS = "EPSG:4326"
FAO_STATICS = "https://storage.googleapis.com/fao-cog-data/WaPOR3.tif"
DEM_TILE = "Copernicus_DSM_COG_10_N25_00_E069_00_DEM"
DEM = f"https://copernicus-dem-30m.s3.amazonaws.com/{DEM_TILE}/{DEM_TILE}.tif"
RADIATION = "https://storage.googleapis.com/fao-gismgr-c3s-data/DATA/C3S/MAPSET/AGERA5-SRF"
CHIRPS = "https://data.chc.ucsb.edu/products/CHIRPS/v2.0/global_daily/cogs/p05"
EARTH_SEARCH = "https://earth-search.aws.element84.com/v1/search"
PLANETARY_COMPUTER = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS_SERVICE = "https://planetarycomputer.microsoft.com/api/sas/v1/token/landsateuwest/landsat-c2"
GDAL_HTTP = {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR", "GDAL_HTTP_TIMEOUT": "30",
             "GDAL_HTTP_MAX_RETRY": "2", "CPL_VSIL_CURL_USE_HEAD": "NO"}
INPUT_GROUPS = {
    "optical": ("ndvi", "r0"), "thermal": ("lst",),
    "meteo": ("t_air", "t_air_min", "t_air_max", "t_dew", "u2m", "v2m", "p_air", "p_air_0", "wv"),
    "terrain": ("z", "slope", "aspect"),
    "statics": ("lw_offset", "lw_slope", "rn_offset", "rn_slope", "t_amp", "t_opt", "vpd_slope"),
    "p": ("p",), "ra_flat": ("ra_flat",),
}
MAX_LST_UNCERTAINTY_K = 3


def buffered(bbox, degrees):
    west, south, east, north = bbox
    return west - degrees, south - degrees, east + degrees, north + degrees


def coordinates(transform, height, width):
    return {"x": transform.c + (np.arange(width) + 0.5) * transform.a,
            "y": transform.f + (np.arange(height) + 0.5) * transform.e}


def model_grid(bbox):
    # ETLook uses its x/y coordinates as longitude/latitude for solar geometry.
    dx, dy = 0.0002, 0.00018
    west, south, east, north = buffered(bbox, 0.0004)
    left, bottom = math.floor(west / dx), math.floor(south / dy)
    right, top = math.ceil(east / dx), math.ceil(north / dy)
    return from_origin(left * dx, top * dy, dx, dy), top - bottom, right - left


def search(endpoint, collection, bbox, start=SOURCE_START):
    params = {"collections": collection, "bbox": ",".join(map(str, bbox)), "limit": 40,
              "datetime": f"{start}T00:00:00Z/{SOURCE_END}T23:59:59Z"}
    response = requests.get(endpoint, params=params, timeout=30)
    response.raise_for_status()
    return sorted(response.json()["features"], key=lambda item: item["properties"]["datetime"])


def raster_array(source, **kwargs):
    return source.read(out_dtype="float32", masked=True, **kwargs).filled(np.nan)


def read_aligned(url, grid, method):
    transform, height, width = grid
    try:
        with rasterio.Env(**GDAL_HTTP), rasterio.open("/vsicurl/" + url) as source:
            with WarpedVRT(source, crs=CRS, transform=transform, height=height, width=width,
                           dtype="float32", nodata=np.nan, resampling=method) as raster:
                return raster_array(raster, indexes=1)
    except rasterio.errors.RasterioError:
        raise RuntimeError(f"Could not read raster {url.split('?')[0]}") from None


def read_native(url, bbox):
    with rasterio.Env(**GDAL_HTTP), rasterio.open("/vsicurl/" + url) as source:
        if source.crs != rasterio.crs.CRS.from_epsg(4326):
            raise ValueError(f"Unexpected CRS in {url}")
        window = from_bounds(*bbox, transform=source.transform)
        left, top = math.floor(window.col_off), math.floor(window.row_off)
        right = math.ceil(window.col_off + window.width)
        bottom = math.ceil(window.row_off + window.height)
        window = Window(left, top, right - left, bottom - top)
        data = raster_array(source, window=window)
        data = data * np.asarray(source.scales)[:, None, None] + np.asarray(source.offsets)[:, None, None]
        return data, source.window_transform(window)


def observations(arrays, dates, grid):
    transform, height, width = grid
    coords = coordinates(transform, height, width) | {"time": np.asarray(dates, dtype="datetime64[ns]")}
    data = xr.Dataset({key: (("time", "y", "x"), np.stack(values)) for key, values in arrays.items()}, coords=coords)
    return data.rio.write_crs(CRS).rio.write_transform(transform)


def scene_time(item):
    return item["properties"]["datetime"].removesuffix("Z")


def sentinel2(item, grid):
    bands = {"blue": "blue", "green": "green", "red": "red", "nir": "nir08"}
    images = {"qa": read_aligned(item["assets"]["scl"]["href"], grid, Resampling.nearest)}
    for variable, key in bands.items():
        asset = item["assets"][key]
        calibration = asset["raster:bands"][0]
        applied = item["properties"]["earthsearch:boa_offset_applied"]
        offset = 0 if applied else calibration.get("offset", 0)
        raw = read_aligned(asset["href"], grid, Resampling.average)
        images[variable] = raw * calibration["scale"] + offset
    ds = observations({key: [value] for key, value in images.items()}, [scene_time(item)], grid)
    for variable in bands:
        ds = apply_qa(ds, variable)
        ds[variable] = ds[variable].where((ds[variable] >= 0) & (ds[variable] <= 1))
    ds = calc_normalized_difference(ds, "ndvi")
    ds = calc_r0(ds, "r0")
    return ds[["ndvi", "r0"]]


def landsat(item, grid, token):
    temperature, flags, error = landsat_bands(item, grid, token)
    band = item["assets"]["lwir11"]["raster:bands"][0]
    valid = usable_temperature(temperature, flags, error)
    temperature = np.where(valid, temperature * band["scale"] + band["offset"], np.nan)
    return observations({"lst": [temperature]}, [scene_time(item)], grid)


def landsat_bands(item, grid, token):
    assets = item["assets"]
    temperature, qa_pixel, qa = (
        read_aligned(f"{assets[key]['href']}?{token}", grid, Resampling.nearest)
        for key in ("lwir11", "qa_pixel", "qa")
    )
    flags = np.nan_to_num(qa_pixel, nan=1).astype(np.uint16)
    uncertainty = qa * assets["qa"]["raster:bands"][0]["scale"]
    return temperature, flags, uncertainty


def usable_temperature(temperature, flags, uncertainty):
    clear = (flags & 0b111111) == 0  # Fill, dilated cloud, cirrus, cloud, cloud shadow and snow.
    return clear & (temperature > 0) & (uncertainty >= 0) & (uncertainty <= MAX_LST_UNCERTAINTY_K)


def daily_rasters(folder, bbox, executor):
    dates = [SOURCE_START + timedelta(days=i) for i in range((SOURCE_END - SOURCE_START).days + 1)]
    recipes = {
        "ra_flat": [(f"{RADIATION}/C3S.AGERA5-SRF.{day}.tif", 1 / 86400) for day in dates],
        "p": [(f"{CHIRPS}/{day:%Y}/chirps-v2.0.{day:%Y.%m.%d}.cog", 1) for day in dates],
    }
    for key, entries in recipes.items():
        if (folder / f"{key}.nc").is_file():
            continue
        outputs = list(executor.map(lambda item: read_native(item[0], buffered(bbox, 0.18)), entries))
        first, transform = outputs[0]
        if any(affine != transform or array.shape != first.shape for array, affine in outputs):
            raise ValueError(f"Daily source grid changed for {key}")
        arrays = [array[0] * scale for (array, _), (_, scale) in zip(outputs, entries)]
        ds = observations({key: arrays}, [f"{day}T12:00:00" for day in dates], (transform, *first.shape[-2:]))
        write_dataset(ds, folder / f"{key}.nc", (key,))
    return {key: [url for url, _ in entries] for key, entries in recipes.items()}


def era5(input_path):
    with np.load(input_path, allow_pickle=False) as raw:
        # pyWaPOR uses right-closed time bins. Put each hourly record just after
        # its start time so midnight falls in the correct reporting period.
        time = raw["valid_time"] + np.timedelta64(1, "us")
        coords = {"x": raw["longitude"], "y": raw["latitude"], "time": time}
        dims = ("time", "y", "x")
        variables = {"t_air": raw["t2m"] - 273.15, "t_dew": raw["d2m"] - 273.15,
                     "p_air": raw["sp"] / 1000, "p_air_0": raw["msl"] / 1000, "wv": raw["tcwv"],
                     "u10m": raw["u10"], "v10m": raw["v10"]}
    ds = xr.Dataset({key: (dims, value) for key, value in variables.items()}, coords=coords)
    for var in ("u2m", "v2m"):
        ds = adjust_wind_height(ds, var)
    ds["t_air_min"] = ds["t_air"]
    ds["t_air_max"] = ds["t_air"]
    ds = ds.drop_vars(["u10m", "v10m"]).sortby("y", ascending=False).rio.write_crs(CRS)
    return ds.rio.write_transform(ds.rio.transform(recalc=True))


def write_dataset(dataset, target, required):
    # pyWaPOR reopens side-loaded files using xarray's CF coordinate decoder.
    # An explicit per-variable grid mapping survives the custom NetCDF encoding.
    dataset = dataset.rio.write_crs(CRS).rio.write_transform(dataset.rio.transform(recalc=True))
    dataset["x"].attrs.update({"standard_name": "longitude", "units": "degrees_east", "axis": "X"})
    dataset["y"].attrs.update({"standard_name": "latitude", "units": "degrees_north", "axis": "Y"})
    for var in required:
        if not np.isfinite(dataset[var]).any():
            raise ValueError(f"No usable {var} observations for {target.name}")
        dataset[var].attrs["grid_mapping"] = "spatial_ref"
    staging = target.with_suffix(".tmp.nc")
    encoding = {var: {"zlib": True, "complevel": 3, "dtype": "float32"} for var in dataset.data_vars}
    dataset.to_netcdf(staging, encoding=encoding)
    staging.replace(target)
    print(f"Prepared {target.name}", flush=True)


def slope_and_aspect(dataset, transform):
    elevation = dataset["z"].values
    height, width = elevation.shape
    source = gdal.GetDriverByName("MEM").Create("", width, height, 1, gdal.GDT_Float32)
    source.SetGeoTransform(transform.to_gdal())
    source.SetProjection(rasterio.crs.CRS.from_epsg(4326).to_wkt())
    source.GetRasterBand(1).SetNoDataValue(-9999)
    source.GetRasterBand(1).WriteArray(np.where(np.isfinite(elevation), elevation, -9999))
    for variable, options in (("slope", {"scale": 111120}), ("aspect", {"zeroForFlat": True})):
        raster = gdal.DEMProcessing("", source, variable, format="MEM", computeEdges=True, **options)
        values = raster.ReadAsArray()
        dataset[variable] = (("y", "x"), np.where(values == raster.GetRasterBand(1).GetNoDataValue(), np.nan, values))
    return dataset


def configuration(folder):
    config = pywapor.Configuration.from_name("WaPOR3_level_3")
    full = config.full
    full["lst"] = copy.deepcopy(full["bt"])
    full["lst"]["variable_enhancers"] = []
    full["lst"]["temporal_interp"] = None
    del full["bt"]

    for filename, variables in INPUT_GROUPS.items():
        source = f"FILE:{folder / (filename + '.nc')}"
        for variable in variables:
            full[variable]["products"] = [{"source": source, "product_name": "none", "enhancers": [],
                                            "is_example": False}]
            if filename == "optical":
                full[variable]["temporal_interp"] = "linear"
    for variable in list(full):
        if variable not in {"se_root", *(v for variables in INPUT_GROUPS.values() for v in variables)}:
            del full[variable]
    config.update_se_root_config()
    config.update_et_look_config()
    config.summarize()
    config.to_json(str(folder / "kgm_public_inputs.json"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aoi", type=Path, required=True)
    parser.add_argument("--folder", type=Path, required=True)
    parser.add_argument("--reuse-staged", action="store_true", help="Configure a transferred set of staged KGM inputs")
    args = parser.parse_args()
    folder = args.folder.resolve()
    folder.mkdir(parents=True, exist_ok=True)
    bbox = json.loads(args.aoi.read_text(encoding="utf-8"))["bbox"]
    grid = model_grid(bbox)
    manifest_path = folder / "kgm_input_sources.json"
    if manifest_path.is_file():
        previous_bbox = json.loads(manifest_path.read_text(encoding="utf-8"))["bounds_wgs84"]
        if previous_bbox != bbox:
            parser.error("Use a new input folder for a different bounding box")
    if args.reuse_staged:
        if not manifest_path.is_file():
            parser.error("A transferred input folder needs kgm_input_sources.json")
        for name, required in INPUT_GROUPS.items():
            with xr.open_dataset(folder / f"{name}.nc", decode_coords="all") as ds:
                if ds.rio.crs is None or ds.rio.crs.to_epsg() != 4326 or not set(required) <= ds.data_vars.keys():
                    raise ValueError(f"Invalid staged input {name}.nc")
        configuration(folder)
        return

    with ThreadPoolExecutor(max_workers=5) as pool:
        optical = pool.submit(search, EARTH_SEARCH, "sentinel-2-l2a", bbox, OPTICAL_START)
        thermal = pool.submit(search, PLANETARY_COMPUTER, "landsat-c2-l2", bbox)
        dem = pool.submit(read_native, DEM, buffered(bbox, 0.005))
        statics = pool.submit(read_native, FAO_STATICS, buffered(bbox, 0.34))
        s2 = [item for item in optical.result() if "_42RWP_" in item["id"]]
        ls = [item for item in thermal.result() if item["id"].startswith(("LC08_L2SP_151043", "LC09_L2SP_151043"))
              and "lwir11" in item["assets"]]
        if len(s2) < 2 or not any(START <= date.fromisoformat(scene_time(item)[:10]) <= END for item in ls):
            raise ValueError("The pilot requires Sentinel-2 observations and at least one Landsat overpass inside the dekad")
        if not (folder / "optical.nc").exists():
            ds = xr.concat(list(pool.map(lambda item: sentinel2(item, grid), s2)), dim="time")
            write_dataset(ds, folder / "optical.nc", ("ndvi", "r0"))
        if not (folder / "thermal.nc").exists():
            token_response = requests.get(SAS_SERVICE, timeout=30)
            token_response.raise_for_status()
            token = token_response.json()["token"]
            ds = xr.concat(list(pool.map(lambda item: landsat(item, grid, token), ls)), dim="time")
            write_dataset(ds, folder / "thermal.nc", ("lst",))
        daily_sources = daily_rasters(folder, bbox, pool)

    for name, (data, transform), variables in (
        ("terrain", dem.result(), ("z",)),
        ("statics", statics.result(), INPUT_GROUPS["statics"]),
    ):
        coords = coordinates(transform, *data.shape[-2:])
        ds = xr.Dataset({var: (("y", "x"), band) for var, band in zip(variables, data)}, coords=coords)
        ds = ds.rio.write_crs(CRS).rio.write_transform(transform)
        target = folder / f"{name}.nc"
        if name == "terrain":
            ds = slope_and_aspect(ds, transform)
            variables += ("slope", "aspect")
        write_dataset(ds, target, variables)
    weather = era5(folder / "era5.npz")
    write_dataset(weather, folder / "meteo.nc", tuple(weather.data_vars))
    configuration(folder)
    manifest = {
        "pilot": f"KGM {START} to {END}", "model": f"pyWaPOR {pywapor.__version__} SE_ROOT and ETLook v3",
        "bounds_wgs84": bbox,
        "sentinel2": {"catalogue": EARTH_SEARCH, "items": [item["id"] for item in s2]},
        "landsat": {"catalogue": PLANETARY_COMPUTER, "items": [item["id"] for item in ls],
                    "method": "Collection 2 Level 2 land surface temperature with QA masks, no thermal sharpening",
                    "thermal_sensor_resolution_m": 100, "max_temperature_uncertainty_k": MAX_LST_UNCERTAINTY_K},
        "era5": "C3S/ECMWF ERA5, Earthmover Icechunk open mirror, single/temporal group",
        "daily_sources": daily_sources, "dem": DEM, "model_coefficients": FAO_STATICS,
        "optical_time_interpolation": "Linear between Sentinel-2 observations",
        "biomass_carbon_fraction_for_reporting": 0.45,
        "reporting_grid": "EPSG:32642, 20 m; Landsat thermal input has a coarser native resolution",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
