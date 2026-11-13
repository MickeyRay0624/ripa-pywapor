"""Trace missing KGM cells to the 5 February 2023 Landsat quality bands."""

import argparse
import json
from pathlib import Path

import numpy as np
import rasterio
import requests

from prepare_kgm_public_inputs import (
    CRS, MAX_LST_UNCERTAINTY_K, PLANETARY_COMPUTER, SAS_SERVICE,
    landsat_bands, model_grid, search, usable_temperature,
)
from pywapor_results import command_grid, regrid


SCENE_ID = "LC08_L2SP_151043_20230205_02_T1"
QA_FLAGS = ("fill_flag", "dilated_cloud", "cirrus", "cloud", "cloud_shadow", "snow")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aoi", required=True, type=Path)
    parser.add_argument("--raster", required=True, type=Path, help="Pilot GeoTIFF with AETI and biomass")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    aoi = json.loads(args.aoi.read_text(encoding="utf-8"))
    grid = command_grid(aoi)
    with rasterio.open(args.raster) as ds:
        if ds.crs != grid.crs or ds.transform != grid.transform or ds.shape != grid.inside.shape:
            raise ValueError("The GeoTIFF does not use the command-area grid")
        missing = grid.inside & np.ma.getmaskarray(ds.read((1, 2), masked=True)).any(axis=0)

    items = search(PLANETARY_COMPUTER, "landsat-c2-l2", aoi["bbox"])
    item = next((item for item in items if item["id"] == SCENE_ID), None)
    if item is None:
        raise ValueError(f"Landsat scene {SCENE_ID} is missing from the catalogue")
    token = requests.get(SAS_SERVICE, timeout=30)
    token.raise_for_status()
    input_grid = model_grid(aoi["bbox"])
    bands = landsat_bands(item, input_grid, token.json()["token"])
    temperature, raw_flags, uncertainty = (regrid(values, input_grid[0], CRS, grid) for values in bands)
    flags = np.nan_to_num(raw_flags, nan=1).astype(np.uint16)
    screened = grid.inside & ~usable_temperature(temperature, flags, uncertainty)
    if not np.array_equal(missing, screened):
        raise ValueError("The output missing cells do not match the Landsat quality mask")

    flags_by_reason = {key: (flags & (1 << bit)) > 0 for bit, key in enumerate(QA_FLAGS)}
    reasons = {
        "missing_temperature_dn": ~np.isfinite(temperature) | (temperature <= 0),
        **flags_by_reason,
        "uncertainty_over_3_k": uncertainty > MAX_LST_UNCERTAINTY_K,
        "missing_uncertainty": ~np.isfinite(uncertainty),
        "invalid_negative_uncertainty": uncertainty < 0,
    }
    audit = {"landsat_item": SCENE_ID, "surface_temperature_asset": item["assets"]["lwir11"]["href"].rsplit("/", 1)[-1]}
    for name, area in (("boundary_edge", grid.edge), ("command_area", grid.inside)):
        affected = missing & area
        audit[name] = {
            "missing_model_pixels": int(affected.sum()),
            "matching_flags_may_overlap": {reason: int((affected & mask).sum()) for reason, mask in reasons.items()},
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
