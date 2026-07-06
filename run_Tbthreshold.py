#!/usr/bin/env python
# coding: utf-8
"""
================================================================================
HPC Tb-Threshold Workflow — F16 SSMIS 37V
================================================================================

PURPOSE
-------
Computes per-pixel brightness temperature (Tb) thresholds from CETB passive
microwave SIR cubes. Thresholds are used downstream for snowmelt onset detection
(MOD) via the SSD-DAV framework.

For each valid (non-Ocean, non-Maritime) land pixel in each cube, the script:
  1. Extracts the full multi-year Tb time series
  2. Computes a smoothed histogram of Tb values
  3. Identifies the valley threshold between the frozen/melting bimodal peaks
  4. Identifies the start and end of the transition range

OUTPUT
------
One CSV per cube:
    TbThresholds_{year_start}-{year_end}_{sensor}_{SiteLabel}_{region}V2Data.csv

Columns:
    Site          — pixel identifier "row,col"
    Snow Class    — Sturm/Liston snow classification
    Latitude      — pixel centre latitude (°N)
    Longitude     — pixel centre longitude (°E)
    x             — cube column index
    y             — cube row index
    Tb Threshold  — valley threshold between frozen/melt Tb peaks (K)
    StartofRange  — lower bound of threshold search range (K)
    EndofRange    — upper bound of threshold search range (K)

USAGE
-----
Single job (all cubes sequentially):
    python run_Tbthresholds.py

SLURM array (one cube per task, recommended):
    sbatch --array=0-N run_thresholds.slurm

DEPENDENCIES
------------
Custom:
    CETB_IO.py                    — read_Tb_whole, years_for, find_cube_offset
    CETB_algorithmsMB_Final.py    — extract_relevant_data, compute_smoothed_histogram,
                                    analyze_histogram, analyze_and_segment_histogram
Data:
    CETB SIR cubes                — /pl/active/PMESDR/CETB_cubesv2/
    SnowClass_Global_1km.nc       — NSIDC-0768 global snow classification

NOTES
-----
- The script auto-discovers all N25_* cube directories for the configured sensor.
- Cubes whose output CSV already exists are skipped automatically — safe to rerun.
- Snow class is loaded once (cached) from the 4GB global NetCDF to save memory.
- Pixels classified as Ocean or Maritime get no threshold (blanked in CSV).

Author: Mahboubeh (Nava) Boueshagh
================================================================================
"""

import os
import csv
import math
import logging
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol
import glob
import xarray as xr
from pyproj import Transformer


# ============================================================
# LOGGING
# Log messages go to both the console and SLURM .err file.
# INFO  = normal progress updates
# WARNING = something unexpected but non-fatal (e.g. no valid sites)
# ERROR = a cube failed — processing continues with the next one
# ============================================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)


# ============================================================
# USER / PATHS
# All output goes to ThresholdDir on scratch storage.
# The script changes to scriptDir so custom imports are found.
# ============================================================
USER = "Mahboubeh"

if USER == "Mahboubeh":
    dataDir      = "/pl/active/PMESDR/CETB_cubesv2"      # CETB cube storage (read-only)
    scriptDir    = "/projects/mabo7890/Melt_Codes/scripts" # where CETB_IO.py etc. live
    outDir       = "/scratch/alpine/mabo7890/Melt_Codes/ipynb_melt_onset_plots"
    ThresholdDir = "/scratch/alpine/mabo7890/Tb_Thresholds_NewCubes"  # output CSVs
    os.makedirs(ThresholdDir, exist_ok=True)
else:
    raise ValueError(f"Unknown user = {USER}")

os.chdir(scriptDir)
log.info("Working directory set to: %s", scriptDir)


# ============================================================
# CUSTOM IMPORTS
# These must be in scriptDir (set by os.chdir above).
# ============================================================
from CETB_IO import read_Tb_whole, coords, years_for, get_sir_info, get_site_boundaries, find_cube_offset
from CETB_algorithmsMB_Final import (
    extract_relevant_data,          # filters Tb time series to relevant obs
    compute_smoothed_histogram,     # builds + smooths Tb histogram
    analyze_histogram,              # finds valley threshold
    analyze_and_segment_histogram,  # finds start/end of threshold range
)


# ============================================================
# SENSOR / SATELLITE CONFIG
# To switch satellites:
#   F13/SSMI  (1995–2009): sat_GRD="F13", sat_SIR="F13", sensor_GRD="SSMI",  sensor_SIR="SSMI"
#   F16/SSMIS (2005–2019): sat_GRD="F16", sat_SIR="F16", sensor_GRD="SSMIS", sensor_SIR="SSMIS"
# ============================================================
sat_GRD    = "F16"    # satellite identifier for GRD files
sat_SIR    = "F16"    # satellite identifier for SIR files (same as GRD here)
sensor_GRD = "SSMIS"  # sensor name for GRD files
sensor_SIR = "SSMIS"  # sensor name for SIR files

channel_GRD = "37V"   # passive microwave channel (37 GHz vertical polarization)
channel_SIR = "37V"
alg_GRD     = "GRD"   # gridded resolution product
alg_SIR     = "SIR"   # scatterometer image reconstruction (higher resolution)
hemName     = "N"      # N = Northern Hemisphere

# Snow class raster — global NSIDC-0768 NetCDF (4GB, loaded once and cached)
snowclass_raster_path = "/projects/mabo7890/Melt_Codes/SnowClass_Global_1km.nc"

# Histogram segmentation bin range (K) — values outside this range are excluded
# from the segmentation step. Extend upper bound if TB values exceed 300K.
SEG_BIN_RANGE = (150, 300)

# Provider / version strings — auto-set from sensor, do not edit
if sensor_GRD in ["SSMI", "SSMIS"]:
    provider = "L1C"
    version  = "v2.*"
elif sensor_GRD == "AMSRE":
    provider = "RSS"
    version  = "v1.3"
elif sensor_GRD == "AMSR2":
    provider = "PPS_XCAL"
    version  = "v1.*"
else:
    raise ValueError(f"Unknown sensor_GRD: {sensor_GRD}")

# Years to process — years_for() returns all valid years for this satellite.
# subYears[0:15] = first 15 years (useful for testing).
# Replace with Years for the full satellite record.
Years     = years_for(sat_GRD)
subYears  = Years[0:15]   # ← change to Years for full run
year_list = list(subYears)
log.info("Years to process: %s → %s (%d years)", year_list[0], year_list[-1], len(year_list))

# Snow class category mapping — from Sturm & Liston (2021) NSIDC-0768
categories = {
    1: "Tundra",
    2: "Boreal Forest",
    3: "Maritime",
    4: "Ephemeral",
    5: "Prairie",
    6: "Montane Forest",
    7: "Ice",
    8: "Ocean",
}


# ============================================================
# SLURM ARRAY HELPERS
# When submitted as a SLURM array job (--array=0-N), each task
# processes exactly one cube. The task index is normalized to
# 0-based using SLURM_ARRAY_TASK_MIN.
# In single-job mode (no SLURM env vars), all cubes are processed
# sequentially by one job.
# ============================================================
def get_slurm_task_info():
    """
    Returns (task_id_0based, n_tasks).
    Works for both SLURM array jobs and single interactive jobs.
    """
    tid    = os.environ.get("SLURM_ARRAY_TASK_ID")
    tmin   = os.environ.get("SLURM_ARRAY_TASK_MIN", "0")
    tmax   = os.environ.get("SLURM_ARRAY_TASK_MAX")
    tcount = os.environ.get("SLURM_ARRAY_TASK_COUNT")

    if tid is None:
        return 0, 1  # not a SLURM array job

    try:
        tid_i  = int(tid)
        tmin_i = int(tmin)
        task0  = tid_i - tmin_i  # normalize to 0-based index

        if tcount is not None:
            n_tasks = int(tcount)
        elif tmax is not None:
            n_tasks = int(tmax) - tmin_i + 1
        else:
            n_tasks = 1

        return max(0, task0), max(1, n_tasks)
    except Exception:
        return 0, 1


def split_work_indices(n_items: int, n_tasks: int, task_id_0based: int) -> np.ndarray:
    """
    Block-partition n_items across n_tasks.
    Returns the array indices belonging to this task.
    Used when splitting pixels within a cube across sub-tasks.
    """
    if n_tasks <= 1:
        return np.arange(n_items, dtype=int)
    block = int(math.ceil(n_items / n_tasks))
    start = task_id_0based * block
    end   = min(n_items, start + block)
    if start >= n_items:
        return np.array([], dtype=int)
    return np.arange(start, end, dtype=int)


# ============================================================
# CUBE DISCOVERY
# Scans the cube storage directory for all N25_* subdirectories
# matching the configured sensor. Skips cubes whose output CSV
# already exists in ThresholdDir — safe to rerun after failures.
# Skips cubes with no matching .nc files (e.g. missing 37V-SIR).
# ============================================================
def discover_cubes():
    """
    Auto-discover all N25_* cube directories for the configured sensor.
    Returns list of (region, SiteLabel) tuples sorted alphabetically.

    SiteLabel naming convention:
        N25_44  → Cube44   (numeric suffix)
        N25_d34 → Cube34d  (letter prefix moved to end)
        N25_a45 → Cube45a  (letter prefix moved to end)
    """
    base_path = f"{dataDir}/{sat_SIR}_{sensor_SIR}/{hemName}/nc_cubes/"
    cube_dirs = sorted(glob.glob(os.path.join(base_path, "cubes_N25_*")))

    log.info("Scanning %d N25_* directories in %s", len(cube_dirs), base_path)

    discovered = []
    skipped    = []

    for cube_dir in cube_dirs:
        region = os.path.basename(cube_dir).replace("cubes_", "")  # e.g. N25_44

        # Build SiteLabel from region suffix
        suffix = region.replace("N25_", "")
        if suffix[0].isalpha():
            site_label = f"Cube{suffix[1:]}{suffix[0]}"  # d34 → Cube34d
        else:
            site_label = f"Cube{suffix}"                  # 44  → Cube44

        # Skip if output CSV already exists (allows safe reruns)
        csv_filename = (
            f"TbThresholds_{year_list[0]}-{year_list[-1]}_"
            f"{sensor_SIR}_{site_label}_{region}V2Data.csv"
        )
        csv_path = os.path.join(ThresholdDir, csv_filename)
        if os.path.exists(csv_path):
            log.info("SKIP %s — CSV already exists: %s", region, csv_filename)
            skipped.append(region)
            continue

        # Verify .nc files exist for this cube
        prefix   = (f"CETB.cubefile.{region}.{sat_SIR}_{sensor_SIR}"
                    f"-{channel_SIR}-{alg_SIR}-{provider}-{version}")
        nc_files = glob.glob(os.path.join(cube_dir, f"{prefix}*.nc"))
        if not nc_files:
            log.warning("No .nc files matching prefix in %s — skipping", cube_dir)
            continue

        discovered.append((region, site_label))

    log.info("Cubes to process : %d", len(discovered))
    log.info("Cubes skipped    : %d %s", len(skipped), skipped)
    for r, s in discovered:
        log.info("  Queued: %s -> %s", r, s)
    return discovered


# ============================================================
# SPATIAL HELPERS
# ============================================================

def latlon_boundary_to_meters(lat_start, lat_end, lon_start, lon_end):
    """
    Convert a lat/lon bounding box to EASE-Grid 2.0 North metre bounds.
    Used when computing cube row/col envelopes from geographic bounds.
    """
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:6931", always_xy=True)
    corners_lon = [lon_start, lon_start, lon_end, lon_end]
    corners_lat = [lat_start, lat_end,   lat_start, lat_end]
    xs, ys = transformer.transform(corners_lon, corners_lat)
    return min(xs), max(xs), min(ys), max(ys)


def meters_to_env_rows_cols(ds_sir, x_min, x_max, y_min, y_max):
    """
    Find (row0, row1, col0, col1) slice indices in a SIR cube
    from metre bounding box coordinates.
    """
    x_coords = ds_sir["x"].values
    y_coords = ds_sir["y"].values  # descending (north-first)
    col_indices = np.where((x_coords >= x_min) & (x_coords <= x_max))[0]
    row_indices = np.where((y_coords >= y_min) & (y_coords <= y_max))[0]
    if len(col_indices) == 0 or len(row_indices) == 0:
        raise ValueError(
            f"No pixels in bounds x=[{x_min:.0f},{x_max:.0f}] y=[{y_min:.0f},{y_max:.0f}]"
        )
    return (
        int(row_indices[0]), int(row_indices[-1] + 1),
        int(col_indices[0]), int(col_indices[-1] + 1),
    )


def get_full_cube_extent(datadir_SIR, prefix_SIR):
    """
    Return (env_rows_cols, lat2d, lon2d) for the ENTIRE cube.
    Uses the first matching .nc file as a spatial reference.
    env_rows_cols = (0, ny, 0, nx) — full extent.
    lat2d, lon2d  — 2D arrays of pixel centre coordinates (ny, nx).
    """
    sir_files = sorted(glob.glob(os.path.join(datadir_SIR, f"{prefix_SIR}*.nc")))
    if not sir_files:
        raise FileNotFoundError(
            f"No .nc files for prefix '{prefix_SIR}' in {datadir_SIR}"
        )
    log.info("Reference file: %s", os.path.basename(sir_files[0]))
    ds = xr.open_dataset(sir_files[0])
    ny, nx = ds.sizes["y"], ds.sizes["x"]
    env_rows_cols = (0, ny, 0, nx)
    lat2d = ds["latitude"].values
    lon2d = ds["longitude"].values
    ds.close()
    log.info("Full cube: ny=%d nx=%d -> %d pixels total", ny, nx, ny * nx)
    return env_rows_cols, lat2d, lon2d


def build_sites_dataframe_for_box(lat2d, lon2d, env_rows_cols):
    """
    Build a flat DataFrame of all pixels in the cube envelope.
    Columns: lat_start, lon_start, y (row), x (col), Site ("row,col").
    This is the pixel table that gets filtered by snow class and
    iterated over in the main processing loop.
    """
    r0, r1, c0, c1 = env_rows_cols
    y_abs = np.arange(r0, r1, dtype=int)
    x_abs = np.arange(c0, c1, dtype=int)
    Y, X  = np.meshgrid(y_abs, x_abs, indexing="ij")
    sites_df = pd.DataFrame({
        "lat_start": lat2d.ravel(),
        "lon_start": lon2d.ravel(),
        "y":         Y.ravel(),
        "x":         X.ravel(),
    })
    sites_df["Site"] = sites_df["y"].astype(str) + "," + sites_df["x"].astype(str)
    return sites_df


# ============================================================
# SNOW CLASS — cached global load
# The NSIDC-0768 NetCDF is 4GB. Loading it once per script run
# (not once per cube) avoids repeated I/O.
# _snow_array is a module-level cache; it persists across
# process_cube() calls within the same job.
# ============================================================
_snow_ds    = None
_snow_array = None

def attach_snow_class(sites_df: pd.DataFrame, raster_path: str) -> pd.DataFrame:
    """
    Attach a SnowClasses column to sites_df using the global snow classification.
    Supports both NetCDF (.nc) and GeoTIFF (.tif) rasters.

    NetCDF path: vectorized xarray nearest-neighbour lookup (fast).
    GeoTIFF path: rasterio pixel sampling (fallback for regional TIFFs).

    Pixels outside the raster coverage return NaN and are dropped
    downstream (e.g. Eurasian cubes when using the NA-only TIF).
    Use SnowClass_Global_1km.nc for full hemispheric coverage.
    """
    global _snow_ds, _snow_array

    if raster_path.endswith(".nc"):
        if _snow_array is None:
            log.info("Loading snow class NetCDF (one-time, ~4GB)...")
            _snow_ds    = xr.open_dataset(raster_path)
            _snow_array = _snow_ds["SnowClass"]
            log.info("Snow class NetCDF loaded and cached.")

        lats = xr.DataArray(sites_df["lat_start"].values, dims="points")
        lons = xr.DataArray(sites_df["lon_start"].values, dims="points")
        vals = _snow_array.sel(lat=lats, lon=lons, method="nearest").values

    else:
        # Fallback: rasterio sampling for GeoTIFF
        coords_xy = list(zip(sites_df["lon_start"].values, sites_df["lat_start"].values))
        with rasterio.open(raster_path) as src:
            vals = np.fromiter(
                (v[0] for v in src.sample(coords_xy)),
                dtype=np.float32, count=len(coords_xy),
            )

    codes = vals.astype(np.int32)
    sites_df = sites_df.copy()
    sites_df["Snow_Class_Code"] = codes
    sites_df["SnowClasses"]     = sites_df["Snow_Class_Code"].map(categories)
    return sites_df


# ============================================================
# PROCESS ONE CUBE
# This is the main processing unit. Called once per cube per task.
# Steps:
#   1. Build pixel table (lat/lon/row/col for all pixels)
#   2. Attach snow class — drop Ocean/Maritime/unclassified
#   3. Read TB cube once for all years (big read, done once)
#   4. Loop over pixels — extract time series, run histogram analysis
#   5. Write results to CSV
# ============================================================
def process_cube(region, SiteLabel):
    """
    Run the full Tb threshold pipeline for one cube.
    Returns the path to the saved CSV, or None if no valid sites.
    """
    cubeType_SIR = f"{channel_SIR}-{alg_SIR}"
    datadir_SIR  = f"{dataDir}/{sat_SIR}_{sensor_SIR}/{hemName}/nc_cubes/cubes_{region}/"
    prefix_SIR   = (f"CETB.cubefile.{region}.{sat_SIR}_{sensor_SIR}"
                    f"-{channel_SIR}-{alg_SIR}-{provider}-{version}")

    # find_cube_offset validates the cube and reads the row/col offset
    # of this cube relative to the global EASE-Grid (needed by read_Tb_whole)
    find_cube_offset(region, cubeDir=datadir_SIR, cubeType=cubeType_SIR, verbose=False)

    # Step 1: Get full cube spatial extent
    env_rows_cols, lat2d, lon2d = get_full_cube_extent(datadir_SIR, prefix_SIR)
    log.info("Cube rows/cols = %s", env_rows_cols)

    # Step 2: Build pixel table and attach snow class
    sites_df = build_sites_dataframe_for_box(lat2d, lon2d, env_rows_cols)
    sites_df = attach_snow_class(sites_df, snowclass_raster_path)

    # Drop pixels with no snow class (outside coverage, or ocean/invalid)
    sites_df = sites_df.dropna(subset=["SnowClasses"]).reset_index(drop=True)
    log.info("Valid sites after snow class filter: %d", len(sites_df))

    if len(sites_df) == 0:
        log.warning("No valid sites for %s — all pixels are Ocean/unclassified. Skipping.", region)
        return None

    # Step 3: Read TB cube ONCE for all years
    # TB shape: (time, ny, nx) — this is the main memory cost (~14GB for 640x640x15yr)
    # Reading once and indexing by pixel avoids repeated disk I/O.
    log.info("Reading TB data for years %s–%s...", year_list[0], year_list[-1])
    data_SIR = read_Tb_whole(
        datadir_SIR, prefix_SIR, year_list,
        env_rows_cols[0], env_rows_cols[1],
        env_rows_cols[2], env_rows_cols[3],
    )
    TB        = np.asarray(data_SIR["TB"])   # (time, ny, nx)
    cal_date  = data_SIR["cal_date"]         # array of dates (used as DataFrame index)
    cal_year  = data_SIR["cal_year"]         # year for each time step
    cal_month = data_SIR["cal_month"]        # month for each time step

    ny_tb, nx_tb = TB.shape[1], TB.shape[2]
    y0, x0       = env_rows_cols[0], env_rows_cols[2]  # origin of the subset
    log.info("TB loaded: %d time steps, %d rows, %d cols", TB.shape[0], ny_tb, nx_tb)

    # Step 4: Per-pixel threshold analysis
    # For each pixel:
    #   - Extract 1D TB time series via direct NumPy indexing (fast)
    #   - Wrap in DataFrame for downstream API compatibility
    #   - Run histogram → smoothing → valley detection → range segmentation
    results_rows = []
    n_sites      = len(sites_df)

    for k, row in enumerate(sites_df.itertuples(index=False), start=1):
        Site_local = row.Site        # "row,col" string
        snow_class = row.SnowClasses # e.g. "Tundra"
        lat_site   = float(row.lat_start)
        lon_site   = float(row.lon_start)
        y_abs      = int(row.y)      # absolute row in full cube
        x_abs      = int(row.x)      # absolute col in full cube

        # Convert to local indices within the TB subset
        yi = y_abs - y0
        xi = x_abs - x0
        if yi < 0 or yi >= ny_tb or xi < 0 or xi >= nx_tb:
            # Pixel is outside the loaded subset — should not happen but guard anyway
            continue

        # Extract time series for this pixel (direct NumPy indexing — O(1))
        tb_series  = TB[:, yi, xi]
        Tb_nearest = pd.DataFrame({"TB": tb_series}, index=cal_date)

        # extract_relevant_data: filters to melt-season months, removes invalid obs
        data = extract_relevant_data(
            Tb_nearest, year_list, cal_year, cal_month, snow_class, Site_local
        )

        # compute_smoothed_histogram: bins Tb values and applies KDE smoothing
        hist, bin_edges, hist_smooth = compute_smoothed_histogram(
            data, snow_class, Site_local
        )

        # analyze_histogram: finds the valley between frozen/melting peaks
        # Returns dict with key "threshold" (K)
        out = analyze_histogram(
            hist, bin_edges, hist_smooth, data,
            Site_local, year_list, sensor_SIR, channel_SIR, snow_class, ThresholdDir
        )

        # analyze_and_segment_histogram: finds start/end of the transition range
        # Returns dict with "effective_end_peak1" (StartofRange) and
        # "effective_start_peak2" (EndofRange)
        seg = analyze_and_segment_histogram(
            data=data, Site=Site_local, year=year_list,
            snow_class=snow_class, ThresholdDir=ThresholdDir,
            bin_range=SEG_BIN_RANGE, plot=False
        )

        StartofRange = EndofRange = None
        if seg and seg.get("peak_distance") is not None:
            StartofRange = seg.get("effective_end_peak1")
            EndofRange   = seg.get("effective_start_peak2")

        results_rows.append([
            Site_local, snow_class, lat_site, lon_site,
            x_abs, y_abs,
            out.get("threshold"), StartofRange, EndofRange,
        ])

        # Progress log every 500 pixels (avoids flooding the log)
        if k % 500 == 0:
            log.info("  %d/%d sites done in %s", k, n_sites, region)

    # Step 5: Write CSV
    csv_filename = (
        f"TbThresholds_{year_list[0]}-{year_list[-1]}_"
        f"{sensor_SIR}_{SiteLabel}_{region}V2Data.csv"
    )
    csv_path = os.path.join(ThresholdDir, csv_filename)
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "Site", "Snow Class", "Latitude", "Longitude",
            "x", "y", "Tb Threshold", "StartofRange", "EndofRange"
        ])
        w.writerows(results_rows)
    log.info("Saved: %s (%d rows)", csv_path, len(results_rows))
    return csv_path


# ============================================================
# MAIN
# Discovers cubes, splits them across SLURM tasks (or runs all
# sequentially in single-job mode), and calls process_cube()
# for each assigned cube. Errors on one cube are caught and
# logged — processing continues with the next cube.
# ============================================================
def main():
    log.info("=" * 60)
    log.info("Tb Threshold Workflow — %s %s %s", sat_SIR, sensor_SIR, channel_SIR)
    log.info("Years: %s → %s (%d years)", year_list[0], year_list[-1], len(year_list))
    log.info("Output: %s", ThresholdDir)
    log.info("=" * 60)

    # Discover all cubes — skips already-processed ones
    cubes = discover_cubes()
    if not cubes:
        log.warning("No cubes to process. Either all done or none found.")
        return

    # SLURM array: each task processes exactly ONE cube
    # Single job: all cubes processed sequentially
    task0, n_tasks = get_slurm_task_info()

    if n_tasks > 1:
        # Array job mode — one cube per task
        if task0 >= len(cubes):
            log.warning(
                "Task %d has no cube assigned (only %d cubes available). Exiting.",
                task0, len(cubes)
            )
            return
        cubes_this_task = [cubes[task0]]
        log.info("SLURM array task %d/%d → cube: %s", task0 + 1, n_tasks, cubes[task0][0])
    else:
        # Single job mode — process all cubes in sequence
        cubes_this_task = cubes
        log.info("Single job mode — processing all %d cubes sequentially", len(cubes))

    # Process each assigned cube
    for region, SiteLabel in cubes_this_task:
        log.info("=" * 60)
        log.info("START  cube: region=%s  SiteLabel=%s", region, SiteLabel)
        try:
            csv_path = process_cube(region, SiteLabel)
            if csv_path:
                log.info("DONE   cube: %s → %s", region, csv_path)
            else:
                log.warning("EMPTY  cube: %s — no valid sites, no CSV written.", region)
        except Exception as e:
            log.error("FAILED cube: %s — %s", region, e)
            import traceback
            traceback.print_exc()
            continue  # move on to the next cube even if this one failed

    log.info("=" * 60)
    log.info("All assigned cubes processed.")


if __name__ == "__main__":
    main()