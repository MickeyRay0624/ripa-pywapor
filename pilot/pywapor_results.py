"""Convert native pyWaPOR dekads to command-area rasters and coverage statistics."""

import csv
import json
import math
import tempfile
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import rasterio
import rioxarray  # Registers the rio accessor on xarray datasets.
import xarray as xr
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_geom
from shapely import segmentize
from shapely.geometry import mapping, shape

from pywapor.general.compositer import time_bins


DAY = timedelta(days=1)
PIXEL_SIZE = 20
CARBON_FRACTION = 0.45
NODATA = -9999.0
LAYERS = (
    ("aeti", "Actual evapotranspiration and interception", "mm"),
    ("dry_biomass", "Estimated dry biomass", "t/ha"),
    ("biomass_water_productivity", "Biomass water productivity", "kg/m3"),
)


@dataclass(frozen=True)
class Period:
    model_start: date
    model_end: date  # Exclusive bound of the model's dekad.
    start: date
    end: date  # Inclusive bound of the reporting period.

    @property
    def days(self):
        return (self.end - self.start).days + 1


@dataclass(frozen=True)
class Grid:
    crs: rasterio.crs.CRS
    transform: rasterio.Affine
    inside: np.ndarray
    edge: np.ndarray


def dekadal_periods(start, end):
    """Use pyWaPOR's own bin boundaries and weight by days inside the request."""
    edges = time_bins([start.isoformat(), (end + DAY).isoformat()], "DEKAD")
    dates = edges.astype("datetime64[D]").astype(object)
    return [Period(left, right, max(start, left), min(end, right - DAY))
            for left, right in zip(dates[:-1], dates[1:])]


def dekad_indices(dataset, periods):
    if "time_bins" not in dataset.coords:
        raise ValueError("Expected a pyWaPOR NetCDF with a time_bins coordinate")
    timestamps = dataset["time_bins"].values
    if timestamps.ndim != 1 or not np.issubdtype(timestamps.dtype, np.datetime64):
        raise ValueError("Expected one-dimensional dates for time_bins")
    days = timestamps.astype("datetime64[D]")
    if np.isnat(days).any() or np.any(timestamps != days):
        raise ValueError("Expected midnight dates for pyWaPOR dekadal bins")
    labels = days.astype(object)
    indices = {day: index for index, day in enumerate(labels)}
    if len(indices) != len(labels) or any(day.day not in (1, 11, 21) for day in labels):
        raise ValueError("Expected unique pyWaPOR dekad start dates")
    missing = [period.model_start.isoformat() for period in periods if period.model_start not in indices]
    if missing:
        raise ValueError(f"NetCDF lacks modelled dekads starting {', '.join(missing)}")
    return [indices[period.model_start] for period in periods]


def fix_year_boundary_doy(dataset, periods):
    """Use the calendar midpoint for dekads that cross New Year."""
    crossing = [period for period in periods if period.model_start.year != period.model_end.year]
    if not crossing:
        return
    indices = dekad_indices(dataset, crossing)
    corrected = dataset["doy"].values.copy()
    for period, index in zip(crossing, indices):
        midpoint = period.model_start + (period.model_end - period.model_start) / 2
        corrected[index] = midpoint.timetuple().tm_yday
    dataset["doy"] = dataset["doy"].copy(data=corrected)


def command_grid(aoi):
    [feature] = [feature for feature in aoi["features"]
                 if feature["properties"].get("role") == "analysis_boundary"]
    boundary = shape(feature["geometry"])
    if boundary.is_empty or not boundary.is_valid or boundary.geom_type not in ("Polygon", "MultiPolygon"):
        raise ValueError("Expected a valid command-area polygon in the AOI GeoJSON")
    crs = rasterio.crs.CRS.from_user_input(aoi["pywapor"]["export_crs"])
    if not crs.is_projected or crs.linear_units != "metre":
        raise ValueError("The 20 m output grid requires a projected CRS measured in metres")
    boundary = shape(transform_geom("EPSG:4326", crs, mapping(segmentize(boundary, 0.0001))))
    left, bottom, right, top = boundary.bounds
    left, bottom = (math.floor(value / PIXEL_SIZE) * PIXEL_SIZE for value in (left, bottom))
    right, top = (math.ceil(value / PIXEL_SIZE) * PIXEL_SIZE for value in (right, top))
    dimensions = ((top - bottom) // PIXEL_SIZE, (right - left) // PIXEL_SIZE)
    transform = from_origin(left, top, PIXEL_SIZE, PIXEL_SIZE)
    inside = geometry_mask([mapping(boundary)], dimensions, transform, invert=True)
    if not inside.any():
        raise ValueError("The command area contains no 20 m pixel centres")
    edge_ring = boundary.difference(boundary.buffer(-PIXEL_SIZE))
    edge = geometry_mask([mapping(edge_ring)], dimensions, transform, invert=True) & inside
    return Grid(crs, transform, inside, edge)


def regrid(values, source_transform, source_crs, grid):
    source = np.asarray(values, dtype=np.float32)
    source = np.where(np.isfinite(source), source, np.nan)
    target = np.full(grid.inside.shape, np.nan, dtype=np.float32)
    reproject(source, target, src_transform=source_transform, src_crs=source_crs, src_nodata=np.nan,
              dst_transform=grid.transform, dst_crs=grid.crs, dst_nodata=np.nan,
              resampling=Resampling.nearest)
    target[~grid.inside] = np.nan
    return target.astype(np.float64)


def productivity(aeti_mm, biomass_t_ha):
    result = np.full(aeti_mm.shape, np.nan)
    np.divide(100 * biomass_t_ha, aeti_mm, out=result,
              where=np.isfinite(biomass_t_ha) & np.isfinite(aeti_mm) & (aeti_mm > 0))
    return result


def coverage_row(label, start, end, grid, arrays, model_footprint):
    aeti, biomass, water_productivity = [np.isfinite(array) for array in arrays]
    both = aeti & biomass
    edge_total = int(grid.edge.sum())
    edge_complete = int((grid.edge & both).sum())
    return {
        "period": label, "start": start.isoformat(), "end": end.isoformat(),
        "inside_pixels": int(grid.inside.sum()), "inside_complete_pixels": int(both.sum()),
        "edge_pixels": edge_total,
        "edge_valid_aeti_pixels": int((grid.edge & aeti).sum()),
        "edge_valid_biomass_pixels": int((grid.edge & biomass).sum()),
        "edge_complete_pixels": edge_complete,
        "edge_missing_pixels": edge_total - edge_complete,
        "edge_outside_model_footprint_pixels": int((grid.edge & ~model_footprint).sum()),
        "edge_missing_inside_model_footprint_pixels": int((grid.edge & model_footprint & ~both).sum()),
        "edge_water_productivity_pixels": int((grid.edge & water_productivity).sum()),
        "edge_complete_percent": round(100 * edge_complete / edge_total, 3) if edge_total else None,
    }


def tiff_profile(grid, count):
    height, width = grid.inside.shape
    return {"driver": "GTiff", "height": height, "width": width, "count": count,
            "crs": grid.crs, "transform": grid.transform, "dtype": "float32", "nodata": NODATA,
            "tiled": True, "compress": "deflate", "predictor": 3}


def write_layers(output, arrays, start, end, offset=0):
    for band, (array, (key, label, unit)) in enumerate(zip(arrays, LAYERS), offset + 1):
        output.write(np.where(np.isfinite(array), array, NODATA).astype(np.float32), band)
        output.set_band_description(band, f"{label} ({start} to {end})")
        output.set_band_unit(band, unit)
        output.update_tags(band, variable=key, start_date=start.isoformat(), end_date=end.isoformat(),
                           days=(end - start).days + 1)


def model_rates(dataset):
    accepted_units = {
        "aeti_24_mm": {"mm day-1", "mm d-1", "mm/day"},
        "npp": {"gc m-2 day-1", "gc m-2 d-1", "gc/m2/day"},
    }
    rates = []
    for variable, allowed in accepted_units.items():
        if variable not in dataset:
            raise ValueError(f"pyWaPOR did not produce {variable}. Check the upstream inputs.")
        rate = dataset[variable]
        unit = rate.attrs.get("units", "").lower().replace("²", "2").replace("^", "").replace("−", "-")
        unit = " ".join(unit.split())
        # The ETLook biomass routine integrates radiation over one day, even when its label omits day-1.
        if variable == "npp" and unit == "gc m-2" and rate.attrs.get("et_look_module") == "pywapor.et_look_v2_v3.biomass":
            unit = "gc m-2 day-1"
        if unit not in allowed:
            raise ValueError(f"Unexpected units on pyWaPOR {variable}: {unit!r}")
        rates.append(rate.transpose("time_bins", "y", "x"))
    return rates


def export_results(netcdf, aoi, periods, output_dir):
    grid = command_grid(aoi)
    output_dir.mkdir(parents=True, exist_ok=True)
    totals = np.zeros((2, *grid.inside.shape), dtype=np.float64)
    season_complete = grid.inside.copy()
    rows = []
    start, end = periods[0].start, periods[-1].end

    with tempfile.TemporaryDirectory(dir=output_dir, prefix=".kgm-") as temporary:
        staging = Path(temporary)
        with xr.open_dataset(netcdf, decode_coords="all", mask_and_scale=True) as dataset:
            aeti_rate, npp_rate = model_rates(dataset)
            indices = dekad_indices(dataset, periods)
            source_crs = dataset.rio.crs
            if source_crs is None:
                raise ValueError("The pyWaPOR NetCDF has no coordinate reference system")
            source_transform = dataset.rio.transform(recalc=True)
            # Test geometry coverage independently of the model's missing data.
            footprint = np.ones((dataset.sizes["y"], dataset.sizes["x"]), dtype=np.float32)
            model_footprint = np.isfinite(regrid(footprint, source_transform, source_crs, grid))
            with rasterio.open(staging / "kgm_dekads_20m.tif", "w", **tiff_profile(grid, 3 * len(periods))) as output:
                output.update_tags(source="pyWaPOR ETLook v3", resampling="nearest")
                for index, (position, period) in enumerate(zip(indices, periods)):
                    aeti = regrid(aeti_rate.isel(time_bins=position).values, source_transform, source_crs, grid) * period.days
                    biomass = regrid(npp_rate.isel(time_bins=position).values, source_transform, source_crs, grid)
                    biomass *= period.days * 0.01 / CARBON_FRACTION  # gC/m²/day to t dry matter/ha.
                    arrays = (aeti, biomass, productivity(aeti, biomass))
                    write_layers(output, arrays, period.start, period.end, offset=3 * index)
                    rows.append(coverage_row("dekad", period.start, period.end, grid, arrays, model_footprint))
                    complete = np.isfinite(aeti) & np.isfinite(biomass)
                    season_complete &= complete
                    totals[0] += np.where(complete, aeti, 0)
                    totals[1] += np.where(complete, biomass, 0)

        aeti, biomass = np.where(season_complete, totals, np.nan)
        arrays = (aeti, biomass, productivity(aeti, biomass))
        rows.append(coverage_row("season", start, end, grid, arrays, model_footprint))
        with rasterio.open(staging / "kgm_season_20m.tif", "w", **tiff_profile(grid, 3)) as output:
            output.update_tags(source="pyWaPOR ETLook v3", resampling="nearest",
                               coverage="Both source variables are valid in every modelled dekad")
            write_layers(output, arrays, start, end)
        with (staging / "kgm_edge_coverage.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

        cell_area = PIXEL_SIZE ** 2
        volume_m3 = float(np.nansum(aeti) * cell_area / 1000)
        dry_mass_kg = float(np.nansum(biomass) * cell_area / 10)
        summary = {
            "source": "pyWaPOR ETLook v3", "netcdf_file": Path(netcdf).name,
            "start": start.isoformat(), "end": end.isoformat(), "days": sum(period.days for period in periods),
            "dekads": len(periods), "crs": str(grid.crs), "pixel_size_m": PIXEL_SIZE,
            "edge_width_inside_boundary_m": PIXEL_SIZE, "resampling": "nearest", "season_coverage": rows[-1],
            "valid_area_ha": int(season_complete.sum()) * cell_area / 10_000,
            "mean_aeti_mm": float(np.nanmean(aeti)) if season_complete.any() else None,
            "mean_dry_biomass_t_ha": float(np.nanmean(biomass)) if season_complete.any() else None,
            "aeti_volume_m3": volume_m3 if season_complete.any() else None,
            "dry_biomass_kg": dry_mass_kg if season_complete.any() else None,
            "biomass_water_productivity_kg_m3": dry_mass_kg / volume_m3 if volume_m3 > 0 else None,
        }
        (staging / "kgm_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        for name in ("kgm_dekads_20m.tif", "kgm_season_20m.tif", "kgm_edge_coverage.csv", "kgm_summary.json"):
            (staging / name).replace(output_dir / name)
    return summary
