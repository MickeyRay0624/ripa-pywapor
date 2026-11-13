#!/usr/bin/env python3
"""Run pyWaPOR and export 20 m command-area results for KGM."""

import argparse
import json
from datetime import date
from hashlib import sha256
from pathlib import Path


MODEL_VARIABLES = ("int_mm", "t_24_mm", "e_24_mm", "aeti_24_mm", "et_ref_24_mm", "se_root", "npp")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aoi", type=Path, required=True)
    parser.add_argument("--start", type=date.fromisoformat, required=True)
    parser.add_argument("--end", type=date.fromisoformat, required=True, help="Last day to include")
    parser.add_argument("--folder", type=Path, required=True, help="Directory for model caches and results")
    parser.add_argument("--output", type=Path, help="Override the results directory")
    parser.add_argument("--configuration-json", type=Path, help="Use a pyWaPOR configuration with staged input files")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--preflight", action="store_true", help="Load the AOI and configuration only")
    mode.add_argument("--netcdf", type=Path, help="Export an existing pyWaPOR dekadal ETLook result")
    args = parser.parse_args()
    if args.start > args.end:
        parser.error("The end date must be on or after the start date")
    if args.netcdf and args.configuration_json:
        parser.error("Use either an existing NetCDF or a model configuration")

    import pywapor
    from pywapor.collect import accounts
    from pywapor_results import dekadal_periods, export_results, fix_year_boundary_doy

    aoi = json.loads(args.aoi.read_text(encoding="utf-8"))
    west, south, east, north = aoi["bbox"]
    if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
        parser.error("The AOI must have a valid longitude and latitude bounding box")
    periods = dekadal_periods(args.start, args.end)
    model_dates = [periods[0].model_start.isoformat(), periods[-1].model_end.isoformat()]
    result_name = f"results_{args.start}_to_{args.end}"

    if args.netcdf:
        netcdf = args.netcdf.resolve()
        output = args.output or args.folder / result_name
    else:
        config = aoi["pywapor"]["configuration"]
        config_fingerprint = sha256(args.configuration_json.read_bytes()).hexdigest() if args.configuration_json else config
        run_key = sha256(json.dumps([pywapor.__version__, config_fingerprint, aoi["bbox"], model_dates]).encode()).hexdigest()[:12]
        model_folder = args.folder if args.preflight else args.folder / f"model_{model_dates[0]}_{run_key}"
        project = pywapor.Project(str(model_folder.resolve()), [west, south, east, north], model_dates)
        if args.configuration_json:
            project.load_configuration(json=str(args.configuration_json.resolve()))
        else:
            project.load_configuration(name=config)
        if args.preflight:
            print(f"Loaded KGM area and {args.configuration_json or config} in pyWaPOR {pywapor.__version__}.")
            return

        result = project.et_look_out
        if result is None or not set(MODEL_VARIABLES).issubset(result.data_vars):
            if result is not None:
                result.close()
            if project.et_look_in is None:
                try:
                    project.set_passwords()
                    if any(product["source"] == "VIIRSL1" for variable in project.configuration.full.values()
                           for product in variable["products"]):
                        # pyWaPOR 3.7.3's default VIIRS downloader uses Earthdata OPeNDAP.
                        accounts.get("NASA")
                except EOFError:
                    parser.exit(2, "\nConfigure the Copernicus Data Space, Climate Data Store and NASA VIIRS "
                                  "credentials in the pyWaPOR runtime before retrying.\n")
                project.download_data()
                project.run_pre_se_root()
                project.run_se_root()
            prepared = project.run_pre_et_look(bin_length="DEKAD")
            fix_year_boundary_doy(prepared, periods)
            result = project.run_et_look(et_look_version="v3", export_vars=list(MODEL_VARIABLES),
                                         chunks={"time_bins": 1, "x": 500, "y": 500})
        if result is None:
            raise RuntimeError("pyWaPOR produced no ETLook results. Check the model log.")
        netcdf = Path(result.encoding["source"])
        result.close()
        output = args.output or model_folder / result_name

    summary = export_results(netcdf, aoi, periods, output)
    edge = summary["season_coverage"]
    print(f"KGM 20 m results: {output.resolve()}")
    print(f"Complete season data in {edge['edge_complete_pixels']}/{edge['edge_pixels']} edge pixels.")


if __name__ == "__main__":
    main()
