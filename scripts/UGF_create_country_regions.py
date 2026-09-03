#!/usr/bin/env python
"""Generate regions.tif / regions_tbl.csv for the UGF study area, with
regions defined by country boundaries.

Modeled directly on scripts/create_toy_regions.py, so the output format
matches what train_policy.py / run_inference.py expect:
    - env_layers/regions.tif       (float32 region-ID raster)
    - env_layers/regions_tbl.csv   (REGION_ID -> NAME lookup table)

"""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize

# --- CONFIG ---------------------------------------------------------------
# Reference raster: defines the exact grid (CRS, transform, shape) regions.tif
# must match. Using one of your AOH rasters, per your species-only pipeline.
REFERENCE_TIF = Path("/Users/kalena/captain3preview/ugf3_data/present_habitat_suitability/Apalis_sharpii.tif").expanduser()

# Existing area mask (built via cn.data_loader.create_mask_from_map()).
MASK_FILE = Path("/Users/kalena/captain3preview/ugf3_data/environmental_layers/area_mask.npy").expanduser()

# Country boundary vector file.
BOUNDARY_FILE = Path(
    "/Users/kalena/captain3preview/ugf3_data/environmental_layers/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp"
).expanduser()


# Attribute column holding the country name.
COUNTRY_NAME_FIELD = "ADMIN"

OUTPUT_TIF = Path("/Users/kalena/captain3preview/ugf3_data/environmental_layers/regions.tif").expanduser()
OUTPUT_CSV = Path("/Users/kalena/captain3preview/ugf3_data/environmental_layers/regions_tbl.csv").expanduser()
NODATA_VALUE = -1.0
# ---------------------------------------------------------------------------


def main():
    with rasterio.open(REFERENCE_TIF) as src:
        profile = src.profile.copy()
        height, width = src.height, src.width
        ref_crs, ref_transform = src.crs, src.transform

    mask = np.load(MASK_FILE)

    boundaries = gpd.read_file(BOUNDARY_FILE)
    if boundaries.crs != ref_crs:
        boundaries = boundaries.to_crs(ref_crs)

    # Clip to the reference raster's extent so we don't rasterize the whole world.
    minx, miny, maxx, maxy = rasterio.transform.array_bounds(
        height, width, ref_transform
    )
    boundaries = boundaries.cx[minx:maxx, miny:maxy].reset_index(drop=True)
    if boundaries.empty:
        raise ValueError(
            "No country boundaries intersect the reference raster's extent. "
            "Check BOUNDARY_FILE CRS / geometry against REFERENCE_TIF."
        )

    # Assign a stable integer ID per country (1..N), preserving name lookup.
    names = boundaries[COUNTRY_NAME_FIELD].tolist()
    ids = list(range(1, len(names) + 1))
    shapes = list(zip(boundaries.geometry, ids))

    region_grid = rasterize(
        shapes=shapes,
        out_shape=(height, width),
        transform=ref_transform,
        fill=NODATA_VALUE,
        dtype="float32",
    )

    # Respect the area mask: cells outside it are not part of any region,
    # regardless of which country polygon they fall in.
    region_grid[mask == 0] = NODATA_VALUE

    profile.update(dtype="float32", count=1, nodata=NODATA_VALUE)
    with rasterio.open(OUTPUT_TIF, "w", **profile) as dst:
        dst.write(region_grid, 1)

    pd.DataFrame({"REGION_ID": ids, "NAME": names}).to_csv(OUTPUT_CSV, index=False)

    # Self-check: confirm the written raster matches the reference georeferencing.
    with rasterio.open(OUTPUT_TIF) as check_src:
        assert check_src.shape == (height, width), "shape mismatch"
        assert check_src.crs == ref_crs, "CRS mismatch"
        assert check_src.transform == ref_transform, "transform mismatch"

    print(f"Wrote {OUTPUT_TIF} and {OUTPUT_CSV} ({len(ids)} regions)")
    for region_id, name in zip(ids, names):
        count = int(np.sum(region_grid == region_id))
        print(f"  {name} (id={region_id}): {count} cells")


if __name__ == "__main__":
    main()
