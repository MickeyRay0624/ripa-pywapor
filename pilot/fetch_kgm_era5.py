"""Stage the public ERA5 hourly input for the February 2023 KGM pilot."""

import argparse
from pathlib import Path

import icechunk
import numpy as np
import xarray as xr


BUCKET = "earthmover-icechunk-era5"
PREFIX = "icechunkV2"
GROUP = "single/temporal"
VARIABLES = ("t2m", "d2m", "sp", "msl", "u10", "v10", "tcwv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    storage = icechunk.s3_storage(bucket=BUCKET, prefix=PREFIX, region="us-east-1", anonymous=True)
    session = icechunk.Repository.open(storage).readonly_session("main")
    with xr.open_zarr(session.store, group=GROUP, consolidated=False, chunks=None) as source:
        # Select contiguous index ranges to avoid fetching a global map for each hour.
        lat = np.flatnonzero((source.latitude.values >= 25) & (source.latitude.values <= 25.75))
        lon = np.flatnonzero((source.longitude.values >= 69) & (source.longitude.values <= 69.75))
        data = source[list(VARIABLES)].isel(latitude=slice(lat[0], lat[-1] + 1),
                                           longitude=slice(lon[0], lon[-1] + 1))
        data = data.sel(valid_time=slice("2023-01-26T00:00:00", "2023-02-17T23:00:00")).load()
        if data.sizes["valid_time"] != 23 * 24 or any(not np.isfinite(data[key]).all() for key in VARIABLES):
            raise ValueError("ERA5 is incomplete for the KGM pilot")
        values = {key: data[key].values for key in (*VARIABLES, "valid_time", "latitude", "longitude")}

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **values)
    temporary.replace(args.output)
    print(f"Saved 552 hours of ERA5 over KGM to {args.output}")


if __name__ == "__main__":
    main()
