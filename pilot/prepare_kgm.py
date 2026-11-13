#!/usr/bin/env python3
"""Export the KGM command boundary and its pyWaPOR bounding box."""

import argparse
import json
from pathlib import Path

import shapefile
from rasterio.crs import CRS
from rasterio.warp import transform_geom
from shapely.geometry import box, mapping, shape
from shapely.ops import unary_union


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", type=Path, help="KGM command-area shapefile")
    parser.add_argument("output", type=Path, help="Output GeoJSON with the boundary and rectangle")
    args = parser.parse_args()

    source_crs = CRS.from_wkt(args.command.with_suffix(".prj").read_text(encoding="utf-8"))
    with shapefile.Reader(str(args.command)) as source:
        polygons = [shape(record.__geo_interface__) for record in source.iterShapes()]
    if not polygons or any(polygon.geom_type not in ("Polygon", "MultiPolygon") or not polygon.is_valid
                           for polygon in polygons):
        raise ValueError("The command shapefile must contain valid polygons")
    native = unary_union(polygons)
    command = shape(transform_geom(source_crs, "EPSG:4326", mapping(native)))
    rectangle = box(*command.bounds)
    target_crs = "EPSG:32642"
    command_utm = shape(transform_geom("EPSG:4326", target_crs, mapping(command)))
    result = {
        "type": "FeatureCollection",
        "bbox": list(command.bounds),
        "pywapor": {
            "configuration": "WaPOR3_level_3",
            "bbox_order": ["west", "south", "east", "north"],
            "bbox_crs": "EPSG:4326",
            "export_crs": target_crs,
            "export_pixel_size_m": 20,
        },
        "features": [
            {
                "type": "Feature",
                "geometry": mapping(command),
                "properties": {"name": "KGM irrigation command", "role": "analysis_boundary",
                               "area_km2": round(command_utm.area / 1_000_000, 3)},
            },
            {
                "type": "Feature",
                "geometry": mapping(rectangle),
                "properties": {"name": "KGM model rectangle", "role": "processing_extent"},
            },
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"KGM box [west, south, east, north]: {result['bbox']}")
    print(f"KGM command area: {result['features'][0]['properties']['area_km2']} km²")


if __name__ == "__main__":
    main()
