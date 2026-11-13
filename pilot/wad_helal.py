#!/usr/bin/env python3
"""Build a seasonal 20 m WaPOR analysis for Wad Helal's wheat fields."""

import argparse
import calendar
import csv
import json
import math
import shutil
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import rasterio
from rasterio.features import geometry_mask
from rasterio.warp import transform_geom
from rasterio.windows import Window, bounds, from_bounds
from shapely.geometry import mapping, shape


START = date(2022, 10, 1)
END = date(2023, 4, 30)
COURSE = "https://raw.githubusercontent.com/wateraccounting/WaPORMOOC/main/data"
BUCKET = "https://storage.googleapis.com/fao-gismgr-wapor-3-data/DATA/WAPOR-3/MOSAICSET"
PRODUCTS = {"L3-AETI-D": 0.1, "L3-NPP-D": 0.001}
CARBON_FRACTION = 0.45
NODATA = -9999.0
BANDS = (
    ("Seasonal water use (AETI)", "mm", "YlGnBu"),
    ("Estimated dry biomass", "t/ha", "YlGn"),
    ("Biomass water productivity", "kg/m³", "magma"),
)


def dekads(start, end):
    month = start.replace(day=1)
    while month <= end:
        last_day = calendar.monthrange(month.year, month.month)[1]
        for part, (first, last) in enumerate(((1, 10), (11, 20), (21, last_day)), 1):
            overlap = (min(end, month.replace(day=last)) - max(start, month.replace(day=first))).days + 1
            if overlap > 0:
                yield f"{month:%Y-%m}-D{part}", overlap
        month = date(month.year + (month.month == 12), month.month % 12 + 1, 1)


def raster_name(product, period):
    return f"WAPOR-3.{product}.GEZ.{period}.tif"


def fetch(url, cache):
    target = cache / url.rsplit("/", 1)[-1]
    if target.is_file():
        return
    for attempt in range(3):
        with tempfile.NamedTemporaryFile(dir=cache, delete=False) as stream:
            temporary = Path(stream.name)
            try:
                request_url = url if attempt == 0 else f"{url}?download={attempt}"
                request = Request(request_url, headers={"User-Agent": "RIPA WaPOR pilot"})
                with urlopen(request, timeout=30) as response:
                    if response.status != 200:
                        raise OSError(f"Expected a complete file, received HTTP {response.status}")
                    shutil.copyfileobj(response, stream)
                    expected = response.headers.get("Content-Length")
                    if expected and stream.tell() != int(expected):
                        raise OSError("Incomplete download")
                stream.close()
                temporary.replace(target)
                return
            except OSError:
                if attempt == 2:
                    raise
                time.sleep(attempt + 1)
            finally:
                temporary.unlink(missing_ok=True)


def read_features(path, source_crs, target_crs):
    collection = json.loads(path.read_text(encoding="utf-8"))
    reproject = rasterio.crs.CRS.from_user_input(source_crs) != target_crs
    features = []
    for feature in collection["features"]:
        geometry = shape(feature["geometry"])
        if reproject and not geometry.is_empty:
            geometry = shape(transform_geom(source_crs, target_crs, mapping(geometry)))
        if geometry.geom_type not in ("Polygon", "MultiPolygon") or not geometry.is_valid:
            raise ValueError(f"Invalid polygon in {path}: {feature['properties']}")
        features.append((geometry, feature["properties"]))
    return features


def enclosing_window(geometry, transform):
    fractional = from_bounds(*geometry.bounds, transform=transform)
    col = math.floor(fractional.col_off)
    row = math.floor(fractional.row_off)
    return Window(
        col, row,
        math.ceil(fractional.col_off + fractional.width) - col,
        math.ceil(fractional.row_off + fractional.height) - row,
    )


def seasonal_totals(cache, periods, window, reference):
    dimensions = (int(window.height), int(window.width))
    complete = np.ones(dimensions, dtype=bool)
    totals = []
    for product, expected_scale in PRODUCTS.items():
        total = np.zeros(dimensions, dtype=np.float64)
        for period, days in periods:
            path = cache / raster_name(product, period)
            with rasterio.open(path) as raster:
                if (raster.crs, raster.transform, raster.shape) != reference:
                    raise ValueError(f"Misaligned source raster: {path.name}")
                if not math.isclose(raster.scales[0], expected_scale) or raster.offsets[0] != 0:
                    raise ValueError(f"Unexpected scale or offset in {path.name}")
                rates = raster.read(1, window=window, boundless=True, masked=True)
            valid = ~np.ma.getmaskarray(rates) & np.isfinite(rates.data)
            complete &= valid
            total += np.where(valid, rates.data, 0) * expected_scale * days
        totals.append(total)
    aeti, npp = totals
    return aeti, npp, complete


def audit_edge(command, complete, transform, window, raster_shape, cache, periods):
    ring = command.difference(command.buffer(-20))
    edge = geometry_mask([mapping(ring)], complete.shape, transform, invert=True)
    rows = np.arange(window.height) + window.row_off
    cols = np.arange(window.width) + window.col_off
    height, width = raster_shape
    footprint = ((rows[:, None] >= 0) & (rows[:, None] < height)
                 & (cols[None, :] >= 0) & (cols[None, :] < width))
    missing = edge & footprint & ~complete
    missing_everywhere = missing.copy()
    if missing.any():
        for product in PRODUCTS:
            for period, _ in periods:
                with rasterio.open(cache / raster_name(product, period)) as raster:
                    values = raster.read(1, window=window, boundless=True, masked=True)
                missing_everywhere &= np.ma.getmaskarray(values) | ~np.isfinite(values.data)

    total = int(edge.sum())
    valid = int((edge & complete).sum())
    return {
        "start": START.isoformat(), "end": END.isoformat(), "dekads": len(periods),
        "edge_definition": "Pixel centres within 20 m inside the command boundary on the native WaPOR grid",
        "edge_pixels": total,
        "edge_complete_pixels": valid,
        "edge_complete_percent": round(100 * valid / total, 2) if total else None,
        "edge_outside_source_footprint_pixels": int((edge & ~footprint).sum()),
        "edge_no_data_in_both_source_products_every_dekad_pixels": int(missing_everywhere.sum()),
    }


def field_statistics(fields, command_mask, complete, transform, aeti, biomass, cell_area):
    rows = []
    for geometry, attributes in fields:
        mask = np.zeros(aeti.shape, dtype=bool)
        if not geometry.is_empty:
            mask = geometry_mask([mapping(geometry)], aeti.shape, transform, invert=True) & command_mask
        valid = mask & complete
        cells = int(mask.sum())
        observed = int(valid.sum())
        if geometry.is_empty:
            status = "missing_geometry"
        elif cells and cells == observed:
            status = "complete"
        else:
            status = "partial" if observed else "missing"
        row = {
            "field_id": attributes["id"],
            "tertiary_unit": attributes["layer"],
            "coverage_status": status,
            "field_pixels": cells,
            "complete_season_pixels": observed,
            "complete_season_pct": round(100 * observed / cells, 2) if cells else None,
            "mapped_area_ha": round(cells * cell_area / 10_000, 4),
            "covered_area_ha": round(observed * cell_area / 10_000, 4),
            "mean_aeti_mm": None,
            "mean_dry_biomass_t_ha": None,
            "consumed_water_m3": None,
            "dry_biomass_kg": None,
            "biomass_water_productivity_kg_m3": None,
        }
        if status == "complete":
            volume = float(aeti[valid].sum() * cell_area / 1000)
            mass = float(biomass[valid].sum() * cell_area)
            row.update(
                mean_aeti_mm=round(float(aeti[valid].mean()), 3),
                mean_dry_biomass_t_ha=round(float(biomass[valid].mean() * 10), 4),
                consumed_water_m3=round(volume, 2),
                dry_biomass_kg=round(mass, 2),
                biomass_water_productivity_kg_m3=round(mass / volume, 5) if volume > 0 else None,
            )
        rows.append(row)
    return rows


def write_tiff(path, arrays, crs, transform):
    height, width = arrays[0].shape
    with rasterio.open(
        path, "w", driver="GTiff", width=width, height=height, count=len(arrays), dtype="float32",
        crs=crs, transform=transform, nodata=NODATA, tiled=True, compress="deflate", predictor=3,
    ) as output:
        for band, (array, (label, unit, _)) in enumerate(zip(arrays, BANDS), 1):
            output.write(np.where(np.isfinite(array), array, NODATA).astype(np.float32), band)
            output.set_band_description(band, label)
            output.set_band_unit(band, unit.replace("³", "3"))
        output.update_tags(
            location="Wad Helal, Gezira, Sudan", start_date=START.isoformat(), end_date=END.isoformat(),
            source="FAO WaPOR v3 Level 3, GEZ, L3-AETI-D and L3-NPP-D",
            coverage="Irrigation block pixels with both inputs for every dekad",
        )


def plot_maps(path, arrays, extent, command, parcels, incomplete, field_counts):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    def draw_outline(axis, geometry, **style):
        polygons = geometry.geoms if geometry.geom_type == "MultiPolygon" else (geometry,)
        for polygon in polygons:
            for ring in (polygon.exterior, *polygon.interiors):
                axis.plot(*ring.xy, **style)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.4), dpi=160)
    fig.patch.set_facecolor("#fbfcfc")
    fig.subplots_adjust(left=0.045, right=0.975, top=0.79, bottom=0.24, wspace=0.14)
    fig.text(0.045, 0.94, "Wad Helal  |  2022–23 wheat season", fontsize=23, weight="bold", color="#142c37")
    fig.text(0.046, 0.878, "Gezira, Sudan  •  1 Oct 2022 – 30 Apr 2023  •  FAO WaPOR v3 Level 3  •  20 m pixels",
             fontsize=11, color="#415966")
    left, right, bottom, top = extent
    for axis, array, (title, unit, cmap) in zip(axes, arrays, BANDS):
        values = array[np.isfinite(array)]
        low, high = np.percentile(values, [2, 98])
        if low == high:
            low, high = float(values.min()), float(values.max())
        axis.set_facecolor("#fbfcfc")
        axis.imshow(np.ma.masked_where(~incomplete, incomplete), extent=extent,
                    cmap=ListedColormap(["#cbd2d5"]), interpolation="nearest")
        plotted = axis.imshow(np.ma.masked_invalid(array), extent=extent, cmap=cmap,
                             vmin=low, vmax=high, interpolation="nearest")
        for parcel, _ in parcels:
            draw_outline(axis, parcel, color="#273c47", linewidth=0.22, alpha=0.75)
        draw_outline(axis, command, color="#253b46", linewidth=0.85)
        axis.set(xlim=(left, right), ylim=(bottom, top), aspect="equal")
        axis.set_axis_off()
        axis.set_title(title, fontsize=12, weight="bold", color="#142c37", pad=8)
        colorbar = fig.colorbar(plotted, ax=axis, orientation="horizontal", pad=0.06,
                                shrink=0.83, fraction=0.078, extend="both")
        colorbar.set_label(unit, fontsize=9, color="#415966")
        colorbar.ax.tick_params(labelsize=8, length=3, colors="#415966")
        colorbar.outline.set_visible(False)
    x0, y0 = left + 0.035 * (right - left), bottom + 0.06 * (top - bottom)
    axes[0].plot([x0, x0 + 2000], [y0, y0], color="#253b46", linewidth=1.4)
    axes[0].text(x0 + 1000, y0 + 200, "2 km", fontsize=8, ha="center", color="#253b46")
    fig.legend(handles=[Patch(facecolor="#cbd2d5", label="Incomplete raster coverage"),
                        Line2D([0], [0], color="#273c47", linewidth=0.75, label="Wheat parcel outlines")],
               loc="lower left", bbox_to_anchor=(0.042, 0.097), frameon=False, fontsize=9, ncols=2)
    fig.text(0.046, 0.065,
             f"Field statistics cover {field_counts['complete']} of {sum(field_counts.values())} wheat parcel records. "
             "Colour ranges use the 2nd to 98th percentiles. The dark line marks the irrigation block.",
             fontsize=9, color="#415966")
    fig.savefig(path, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    periods = list(dekads(START, END))
    urls = [f"{COURSE}/{name}.geojson" for name in ("Wad_Helal", "WH_Fields")]
    urls.extend(f"{BUCKET}/{product}/{raster_name(product, period)}"
                for product in PRODUCTS for period, _ in periods)
    download_start = time.perf_counter()
    cache_hits = sum((args.cache / url.rsplit("/", 1)[-1]).is_file() for url in urls)
    with ThreadPoolExecutor(max_workers=6) as pool:
        list(pool.map(lambda url: fetch(url, args.cache), urls))
    download_seconds = time.perf_counter() - download_start

    analysis_start = time.perf_counter()
    with rasterio.open(args.cache / raster_name("L3-AETI-D", periods[0][0])) as source:
        crs, transform = source.crs, source.transform
        if crs.to_epsg() != 32636 or source.res != (20.0, 20.0) or not transform.is_rectilinear:
            raise ValueError("Expected the native Gezira 20 m UTM grid")
        reference = (crs, transform, source.shape)
    [(command, _)] = read_features(args.cache / "Wad_Helal.geojson", "EPSG:4326", crs)
    fields = read_features(args.cache / "WH_Fields.geojson", "EPSG:32636", crs)
    mapped_fields = [(geometry, attrs) for geometry, attrs in fields if not geometry.is_empty]
    window = enclosing_window(command, transform)
    transform = rasterio.windows.transform(window, transform)
    aeti, npp, complete = seasonal_totals(args.cache, periods, window, reference)
    edge_report = audit_edge(command, complete, transform, window, reference[2], args.cache, periods)
    cell_area = abs(transform.a * transform.e)
    command_mask = geometry_mask([mapping(command)], aeti.shape, transform, invert=True)
    wheat_mask = geometry_mask([mapping(p) for p, _ in mapped_fields], aeti.shape, transform, invert=True) & command_mask
    map_covered = command_mask & complete
    covered = wheat_mask & complete
    if not covered.any():
        raise ValueError("No wheat pixels have a complete AETI and NPP series")

    biomass = npp / CARBON_FRACTION / 1000  # kg of dry matter per square metre
    productivity = np.full(aeti.shape, np.nan)
    np.divide(biomass, aeti / 1000, out=productivity, where=map_covered & (aeti > 0))
    arrays = (np.where(map_covered, aeti, np.nan), np.where(map_covered, 10 * biomass, np.nan), productivity)
    rows = field_statistics(fields, command_mask, complete, transform, aeti, biomass, cell_area)
    field_counts = Counter(row["coverage_status"] for row in rows)
    volume = float(aeti[covered].sum() * cell_area / 1000)
    mass = float(biomass[covered].sum() * cell_area)
    analysis_seconds = time.perf_counter() - analysis_start

    write_tiff(args.output / "wad_helal_seasonal_20m.tif", arrays, crs, transform)
    with (args.output / "wad_helal_fields.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output / "wad_helal_edge_coverage.json").write_text(json.dumps(edge_report, indent=2) + "\n", encoding="utf-8")
    west, south, east, north = bounds(Window(0, 0, window.width, window.height), transform)
    plot_maps(args.output / "wad_helal_maps.png", arrays, (west, east, south, north), command, mapped_fields,
              command_mask & ~complete, field_counts)
    manifest = {
        "site": "Wad Helal, Gezira irrigation scheme, Sudan",
        "season_start": START.isoformat(),
        "season_end": END.isoformat(),
        "dekads": dict(periods),
        "days": sum(days for _, days in periods),
        "crs": str(crs),
        "pixel_size_m": 20,
        "command_area_ha": round(command.area / 10_000, 2),
        "complete_command_area_ha": round(int(map_covered.sum()) * cell_area / 10_000, 2),
        "mapped_wheat_area_ha": round(int(wheat_mask.sum()) * cell_area / 10_000, 2),
        "complete_wheat_area_ha": round(int(covered.sum()) * cell_area / 10_000, 2),
        "field_counts": dict(field_counts),
        "mean_aeti_mm": round(float(aeti[covered].mean()), 3),
        "mean_dry_biomass_t_ha": round(float(biomass[covered].mean() * 10), 4),
        "aeti_volume_m3": round(volume, 2),
        "dry_biomass_kg": round(mass, 2),
        "biomass_water_productivity_kg_m3": round(mass / volume, 5) if volume > 0 else None,
        "summary_statistics_scope": "Complete pixels in the IHE Delft wheat parcels within the irrigation block.",
        "method": {
            "AETI": "Actual evapotranspiration plus interception. Multiply each mean daily rate by the dekad's days.",
            "biomass": "Multiply each daily NPP rate by the dekad's days. Divide cumulative carbon by 0.45 to estimate dry matter.",
            "water_productivity": "Dry biomass divided by AETI water volume. Field and overall values use the ratio of totals.",
            "coverage": "Map uses the irrigation block boundary. Field statistics use the pixel-centre intersection of the wheat parcels and the block. Both products must be valid in every dekad. Incomplete fields have blank seasonal statistics.",
        },
        "timing": {
            "source_files_cached_at_start": cache_hits,
            "source_check_and_download_seconds": round(download_seconds, 3),
            "local_analysis_seconds": round(analysis_seconds, 3),
            "local_analysis_and_export_seconds": round(time.perf_counter() - analysis_start, 3),
        },
        "source_urls": urls,
    }
    (args.output / "wad_helal_summary.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: manifest[key] for key in (
        "field_counts", "complete_wheat_area_ha", "mean_aeti_mm", "mean_dry_biomass_t_ha",
        "biomass_water_productivity_kg_m3", "timing",
    )}, indent=2))


if __name__ == "__main__":
    main()
