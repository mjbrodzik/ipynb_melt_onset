#!/usr/bin/env python
# coding: utf-8
"""
================================================================================
HPC DAV-Threshold Workflow — F16 SSMIS 37V
================================================================================

PURPOSE
-------
Computes per-pixel Daily Amplitude Variation (DAV) thresholds from CETB passive
microwave SIR cubes. DAV thresholds are used downstream for snowmelt onset
detection (MOD) via the SSD-DAV framework, alongside the Tb thresholds produced
by run_Tbthresholds.py.

DAV is defined as the absolute value of the day-to-day difference in Tb:
    DAV[t] = |TB[t] - TB[t-1]|

High DAV values indicate rapidly changing surface conditions — characteristic
of the melt transition. The p90 and p95 thresholds of the annual DAV distribution
define the detection threshold for melt onset.

For each valid (non-Ocean, non-Maritime) land pixel in each cube, the script:
  1. Reads the full multi-year Tb time series
  2. Computes DAV for the entire time axis (once, vectorized)
  3. Loops over each year:
       - Extracts the annual DAV series
       - Computes p90, p95, and average thresholds
       - Writes results to a per-year CSV

OUTPUT
------
One CSV per year per cube:
    DAVThresholds_{year}_{sensor}_{SiteLabel}_{region}V2Data.csv

Columns:
    Site                          — pixel identifier "row,col"
    Snow Class                    — Sturm/Liston snow classification
    Latitude                      — pixel centre latitude (°N)
    Longitude                     — pixel centre longitude (°E)
    x                             — cube column index
    y                             — cube row index
    DAV_Threshold_90_Percentile   — p90 of annual DAV distribution (K)
    DAV_Threshold_95_Percentile   — p95 of annual DAV distribution (K)
    DAV_Threshold_Average         — mean of p90 and p95 (K)

JOIN KEYS (shared with Tb CSV for MOD detection):
    Site, Snow Class, Latitude, Longitude, x, y

SKIP LOGIC
----------
- A cube is skipped if ALL year CSVs already exist for it.
- Individual year CSVs are skipped if they already exist.
- This makes reruns safe after partial failures.

USAGE
-----
Single job (all cubes sequentially):
    python run_DAVthresholds.py

SLURM array (one cube per task, recommended):
    sbatch --array=0-N run_DAVthresholds.slurm

DEPENDENCIES
------------
Custom:
    CETB_IO.py                 — read_Tb_whole, years_for, find_cube_offset
    CETB_algorithmsMB_Final.py — calc_DAV

Data:
    CETB SIR cubes             — /pl/active/PMESDR/CETB_cubesv2/
    SnowClass_Global_1km.nc    — NSIDC-0768 global snow classification

NOTES
-----
- DAV is computed once over the full time axis (all years) using calc_DAV(TB).
  This is more efficient than computing it year by year.
- calc_DAV inserts a 0 at time step 0 (no predecessor), so the first observation
  of each year that immediately follows a year boundary may have a spurious DAV=0.
  This is handled naturally by the quantile computation (0 values are included).
- Ocean and Maritime pixels are kept in the CSV but thresholds are set to None
  (blanked). This preserves spatial coverage for GeoTIFF production.
- Snow class is loaded once (cached) from the 4GB global NetCDF.
- The pixel_cache dict is reset each year — it has no effect since every pixel
  is unique per year. It is retained for structural consistency with earlier versions.

Author: Mahboubeh (Nava) Boueshagh
================================================================================
"""

import os
import csv
import math
import logging
import glob
from datetime import datetime

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from pyproj import Transformer


# ============================================================
# LOGGING
# Log messages go to the console and SLURM .err file.
# INFO    = normal progress (cube start/end, year progress)
# WARNING = non-fatal issues (no valid sites, missing years)
# ERROR   = cube failed — processing continues with next cube
# ============================================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)


# ============================================================
# USER / PATHS
# ThresholdDir receives all output CSVs (one per year per cube).
# scriptDir must contain CETB_IO.py and CETB_algorithmsMB_Final.py.
# ============================================================
USER = "Mahboubeh"

if USER == "Mahboubeh":
    dataDir      = "/pl/active/PMESDR/CETB_cubesv2"       # CETB cube storage (read-only)
    scriptDir    = "/projects/mabo7890/Melt_Codes/scripts"  # custom module location
    ThresholdDir = "/scratch/alpine/mabo7890/DAV_Thresholds" # output CSVs
    os.makedirs(ThresholdDir, exist_ok=True)
else:
    raise ValueError(f"Unknown user = {USER}")

os.chdir(scriptDir)
log.info("Working directory set to: %s", scriptDir)


# ============================================================
# CUSTOM IMPORTS
# Must be in scriptDir (set by os.chdir above).
# calc_DAV: computes |TB[t] - TB[t-1]| over the full 3D array.
# ============================================================
from CETB_IO import read_Tb_whole, years_for, find_cube_offset
from CETB_algorithmsMB_Final import calc_DAV


# ============================================================
# SENSOR / SATELLITE CONFIG
# To switch satellites:
#   F13/SSMI  (1995–2009): sat="F13", sensor="SSMI"
#   F16/SSMIS (2005–2019): sat="F16", sensor="SSMIS"
# ============================================================
sat_GRD    = "F16"    # satellite ID for GRD files
sat_SIR    = "F16"    # satellite ID for SIR files
sensor_GRD = "SSMIS"  # sensor name for GRD files
sensor_SIR = "SSMIS"  # sensor name for SIR files

channel_GRD = "37V"   # 37 GHz vertical polarization
channel_SIR = "37V"
alg_GRD     = "GRD"   # gridded resolution product
alg_SIR     = "SIR"   # scatterometer image reconstruction (higher res)
hemName     = "N"      # Northern Hemisphere

# Snow class raster — global NSIDC-0768 NetCDF (4GB, loaded once and cached)
snowclass_raster_path = "/projects/mabo7890/Melt_Codes/SnowClass_Global_1km.nc"

# Years — years_for() returns all valid years for this satellite.
# subYears[0:15] = first 15 years (good for testing).
# Replace with Years for the full record.
Years     = years_for(sat_GRD)
subYears  = Years[0:15]   # ← change to Years for full run
year_list = list(subYears)
log.info("Years to process: %s → %s (%d years)", year_list[0], year_list[-1], len(year_list))

# Provider / version — auto-set from sensor, do not edit
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

# Snow class category mapping — Sturm & Liston (2021) NSIDC-0768
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

# Pixels in BLANK_SNOW_CLASSES are kept in the CSV but get None thresholds.
# This preserves spatial completeness for GeoTIFF production while correctly
# marking non-snow-covered surfaces as having no DAV threshold.
BLANK_SNOW_CLASSES = {"Ocean", "Maritime"}


# ============================================================
# SLURM ARRAY HELPERS
# In array mode (--array=0-N), each task processes one cube.
# SLURM_ARRAY_TASK_ID is normalized to 0-based using TASK_MIN.
# In single-job mode (no SLURM env vars), all cubes run sequentially.
# ============================================================
def get_slurm_task_info():
    """
    Returns (task_id_0based, n_tasks).
    Works for both SLURM array jobs and single interactive runs.
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
        task0  = tid_i - tmin_i  # normalize to 0-based

        if tcount is not None:
            n_tasks = int(tcount)
        elif tmax is not None:
            n_tasks = int(tmax) - tmin_i + 1
        else:
            n_tasks = 1

        return max(0, task0), max(1, n_tasks)
    except Exception:
        return 0, 1


# ============================================================
# CUBE DISCOVERY
# Scans cube storage for all N25_* directories matching the
# configured sensor. A cube is skipped only if ALL its year
# CSVs already exist — partial reruns pick up missing years.
# ============================================================
def discover_cubes():
    """
    Auto-discover all N25_* cube directories for the configured sensor.
    Returns list of (region, SiteLabel) tuples sorted alphabetically.

    SiteLabel naming convention (mirrors Tb script):
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
        region = os.path.basename(cube_dir).replace("cubes_", "")

        suffix = region.replace("N25_", "")
        if suffix[0].isalpha():
            site_label = f"Cube{suffix[1:]}{suffix[0]}"
        else:
            site_label = f"Cube{suffix}"

        # Skip only if ALL year CSVs already exist for this cube.
        # If even one year is missing the cube is re-queued; the
        # per-year skip inside process_cube handles existing years.
        all_done = all(
            os.path.exists(os.path.join(
                ThresholdDir,
                f"DAVThresholds_{yr}_{sensor_SIR}_{site_label}_{region}V2Data.csv"
            ))
            for yr in year_list
        )
        if all_done:
            log.info("SKIP %s — all %d year CSVs already exist", region, len(year_list))
            skipped.append(region)
            continue

        # Verify .nc files exist (e.g. N25_77 had 0 SIR files before Molly reproduced it)
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
# SPATIAL HELPERS — identical to Tb script
# ============================================================
def get_full_cube_extent(datadir_SIR, prefix_SIR):
    """
    Return (env_rows_cols, lat2d, lon2d) for the ENTIRE cube.
    Uses try/finally to guarantee the NetCDF file is closed
    even if an error occurs during array extraction.
    env_rows_cols = (0, ny, 0, nx).
    lat2d, lon2d  = 2D pixel centre coordinate arrays (ny, nx).
    """
    sir_files = sorted(glob.glob(os.path.join(datadir_SIR, f"{prefix_SIR}*.nc")))
    if not sir_files:
        raise FileNotFoundError(
            f"No .nc files for prefix '{prefix_SIR}' in {datadir_SIR}"
        )
    log.info("Reference file: %s", os.path.basename(sir_files[0]))
    ds = xr.open_dataset(sir_files[0])
    try:
        ny, nx = ds.sizes["y"], ds.sizes["x"]
        env_rows_cols = (0, ny, 0, nx)
        lat2d = ds["latitude"].values
        lon2d = ds["longitude"].values
    finally:
        ds.close()
    log.info("Full cube: ny=%d nx=%d -> %d pixels total", ny, nx, ny * nx)
    return env_rows_cols, lat2d, lon2d


def build_sites_dataframe_for_box(lat2d, lon2d, env_rows_cols):
    """
    Build a flat DataFrame of all pixels in the cube envelope.
    Columns: lat_start, lon_start, y (row), x (col), Site ("row,col").
    Identical to Tb script — same join keys for downstream merging.
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
# SNOW CLASS — cached global load (identical to Tb script)
# Loaded once per job run, not once per cube or year.
# ============================================================
_snow_array = None

def attach_snow_class(sites_df: pd.DataFrame, raster_path: str) -> pd.DataFrame:
    """
    Attach a SnowClasses column using the global snow classification.
    NetCDF: vectorized xarray nearest-neighbour lookup (fast, global coverage).
    GeoTIFF: rasterio pixel sampling (fallback, regional coverage only).
    Pixels outside raster extent return NaN — dropped downstream.
    """
    global _snow_array

    if raster_path.endswith(".nc"):
        if _snow_array is None:
            log.info("Loading snow class NetCDF (one-time, ~4GB)...")
            _snow_ds    = xr.open_dataset(raster_path)
            _snow_array = _snow_ds["SnowClass"]
            log.info("Snow class cached.")

        lats = xr.DataArray(sites_df["lat_start"].values, dims="points")
        lons = xr.DataArray(sites_df["lon_start"].values, dims="points")
        vals = _snow_array.sel(lat=lats, lon=lons, method="nearest").values
    else:
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
# DOY HELPER
# Converts cal_date (which read_Tb_whole may return as pandas
# Timestamps, datetime objects, numpy datetime64, or YYYYMMDD
# integers) to day-of-year (1–366).
# Robust fallback chain handles all observed formats.
# ============================================================
def cal_date_to_doy(cal_date_array):
    """
    Convert cal_date array to DOY (1–366).
    Tries pd.Timestamp first (handles Timestamps, datetime, datetime64),
    then falls back to YYYYMMDD integer parsing.
    Returns -1 for any date that cannot be parsed.
    Logs the first element so the format is visible in the job log.
    """
    cal_date_arr = np.asarray(cal_date_array)
    doy = np.zeros(len(cal_date_arr), dtype=np.int32)

    if len(cal_date_arr) > 0:
        log.info("cal_date sample[0] = %s (type=%s)", cal_date_arr[0], type(cal_date_arr[0]))

    for i, d in enumerate(cal_date_arr):
        try:
            ts = pd.Timestamp(d)
            doy[i] = ts.day_of_year
        except Exception:
            try:
                dt = datetime.strptime(str(int(d)), "%Y%m%d")
                doy[i] = dt.timetuple().tm_yday
            except Exception:
                doy[i] = -1
    return doy


# ============================================================
# PROCESS ONE CUBE
# Steps:
#   1. Build pixel table and attach snow class
#   2. Read TB cube ONCE for all years
#   3. Compute DAV ONCE over the full time axis (vectorized)
#   4. Loop over years:
#        - Skip year if CSV already exists
#        - Extract annual DAV series per pixel
#        - Compute p90, p95, average
#        - Write per-year CSV
# ============================================================
def process_cube(region, SiteLabel):
    """
    Run the DAV threshold pipeline for one cube across all years.
    Returns list of paths to saved CSVs (one per year processed).
    Returns [] if no valid sites or all years already done.
    """
    cubeType_SIR = f"{channel_SIR}-{alg_SIR}"
    datadir_SIR  = f"{dataDir}/{sat_SIR}_{sensor_SIR}/{hemName}/nc_cubes/cubes_{region}/"
    prefix_SIR   = (f"CETB.cubefile.{region}.{sat_SIR}_{sensor_SIR}"
                    f"-{channel_SIR}-{alg_SIR}-{provider}-{version}")

    # Validate cube and read row/col offset relative to global EASE-Grid
    find_cube_offset(region, cubeDir=datadir_SIR, cubeType=cubeType_SIR, verbose=False)

    # Step 1: Spatial setup — identical to Tb script
    env_rows_cols, lat2d, lon2d = get_full_cube_extent(datadir_SIR, prefix_SIR)
    log.info("Cube rows/cols = %s", env_rows_cols)

    sites_df = build_sites_dataframe_for_box(lat2d, lon2d, env_rows_cols)
    sites_df = attach_snow_class(sites_df, snowclass_raster_path)
    sites_df = sites_df.dropna(subset=["SnowClasses"]).reset_index(drop=True)
    log.info("Valid sites after snow class filter: %d", len(sites_df))

    if len(sites_df) == 0:
        log.warning("No valid sites for %s — all pixels are Ocean/unclassified. Skipping.", region)
        return []

    # Step 2: Read TB cube ONCE for all years
    # Shape: (time, ny, nx). Reading all years at once avoids repeated disk I/O.
    # Memory: ~14GB for a 640×640 cube × 15 years × daily obs.
    log.info("Reading TB data for years %s–%s...", year_list[0], year_list[-1])
    data_SIR = read_Tb_whole(
        datadir_SIR, prefix_SIR, year_list,
        env_rows_cols[0], env_rows_cols[1],
        env_rows_cols[2], env_rows_cols[3],
    )
    TB       = np.asarray(data_SIR["TB"])        # (t, ny, nx)
    cal_year = np.asarray(data_SIR["cal_year"])  # (t,) — year for each time step

    ny_tb, nx_tb = TB.shape[1], TB.shape[2]
    y0, x0       = env_rows_cols[0], env_rows_cols[2]
    log.info("TB loaded: %d time steps, %d rows, %d cols", TB.shape[0], ny_tb, nx_tb)

    # Step 3: Compute DAV ONCE over the full time axis
    # calc_DAV(TB) = |TB[t] - TB[t-1]| for all t, vectorized over (ny, nx).
    # First time step gets DAV=0 (no predecessor — inserted by calc_DAV).
    # Result shape: (t, ny, nx) — same as TB.
    log.info("Computing DAV over full time axis...")
    DAV = np.asarray(calc_DAV(TB))  # (t, ny, nx)
    if DAV.shape != TB.shape:
        log.warning("DAV shape %s != TB shape %s — unexpected.", DAV.shape, TB.shape)
    log.info("DAV computed: shape=%s", DAV.shape)

    saved_paths = []
    n_sites     = len(sites_df)

    # Step 4: Loop over years — one CSV per year per cube
    for year in year_list:

        # Per-year skip: allows reruns to pick up only missing years
        csv_filename = (
            f"DAVThresholds_{year}_"
            f"{sensor_SIR}_{SiteLabel}_{region}V2Data.csv"
        )
        csv_path = os.path.join(ThresholdDir, csv_filename)
        if os.path.exists(csv_path):
            log.info("SKIP year=%d — CSV already exists", year)
            continue

        year_mask = (cal_year == year)
        if not np.any(year_mask):
            log.warning("No time steps for year=%d in cube %s — skipping.", year, region)
            continue

        log.info("--- year=%d: %d time steps ---", year, year_mask.sum())

        results_rows  = []
        pixel_cache   = {}   # cache p90/p95/avg by pixel key within this year
        blanked_count = 0    # pixels kept with None thresholds (Ocean/Maritime)
        outside_count = 0    # pixels outside the TB subset bounds (should be 0)
        all_nan_count = 0    # pixels with no valid DAV observations for this year

        for k, row in enumerate(sites_df.itertuples(index=False), start=1):
            Site_local = row.Site        # "row,col"
            snow_class = row.SnowClasses # e.g. "Tundra"
            lat_site   = float(row.lat_start)
            lon_site   = float(row.lon_start)
            y_abs      = int(row.y)
            x_abs      = int(row.x)

            # Convert absolute cube indices to local TB subset indices
            yi = y_abs - y0
            xi = x_abs - x0

            if yi < 0 or yi >= ny_tb or xi < 0 or xi >= nx_tb:
                # Should not occur for full-cube reads — guarded for safety
                outside_count += 1
                continue

            # Ocean/Maritime: keep in CSV for spatial completeness, blank thresholds
            if snow_class in BLANK_SNOW_CLASSES:
                blanked_count += 1
                results_rows.append([
                    Site_local, snow_class, lat_site, lon_site,
                    x_abs, y_abs, None, None, None,
                ])
                continue

            # pixel_cache avoids recomputing the same pixel if it appears twice
            # (should not happen for a regular grid, but retained for safety)
            SIR_key = f"{y_abs},{x_abs}"

            if SIR_key in pixel_cache:
                p90, p95, DAV_avg = pixel_cache[SIR_key]
            else:
                # Extract annual DAV series for this pixel
                dav_year = DAV[:, yi, xi][year_mask]  # 1D array for this year

                if np.all(np.isnan(dav_year)):
                    # No valid observations for this pixel-year (e.g. polar night gaps)
                    all_nan_count += 1
                    results_rows.append([
                        Site_local, snow_class, lat_site, lon_site,
                        x_abs, y_abs, None, None, None,
                    ])
                    continue

                # DAV thresholds: p90 and p95 of the annual distribution
                # DAV_avg = (p90 + p95) / 2 — used as the primary detection threshold
                p90     = float(np.nanquantile(dav_year, 0.90))
                p95     = float(np.nanquantile(dav_year, 0.95))
                DAV_avg = (p90 + p95) / 2.0
                pixel_cache[SIR_key] = (p90, p95, DAV_avg)

            results_rows.append([
                Site_local, snow_class, lat_site, lon_site,
                x_abs, y_abs, p90, p95, DAV_avg,
            ])

            if k % 500 == 0:
                log.info("  year=%d: %d/%d sites done", year, k, n_sites)

        log.info(
            "Done %s year=%d | blanked=%d | outside=%d | all-NaN=%d | cached=%d",
            region, year, blanked_count, outside_count, all_nan_count, len(pixel_cache)
        )

        # Write per-year CSV
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([
                "Site", "Snow Class", "Latitude", "Longitude", "x", "y",
                "DAV_Threshold_90_Percentile",
                "DAV_Threshold_95_Percentile",
                "DAV_Threshold_Average",
            ])
            w.writerows(results_rows)

        log.info("Saved: %s (%d rows)", csv_path, len(results_rows))
        saved_paths.append(csv_path)

    return saved_paths


# ============================================================
# MAIN
# Discovers cubes, assigns one cube per SLURM task (or runs all
# sequentially), calls process_cube() per cube.
# Errors on one cube are caught — processing continues with next.
# ============================================================
def main():
    log.info("=" * 60)
    log.info("DAV Threshold Workflow — %s %s %s", sat_SIR, sensor_SIR, channel_SIR)
    log.info("Years: %s → %s (%d years)", year_list[0], year_list[-1], len(year_list))
    log.info("Output: %s", ThresholdDir)
    log.info("=" * 60)

    cubes = discover_cubes()
    if not cubes:
        log.warning("No cubes to process. Either all done or none found.")
        return

    # SLURM array: one cube per task | Single job: all cubes sequentially
    task0, n_tasks = get_slurm_task_info()

    if n_tasks > 1:
        if task0 >= len(cubes):
            log.warning(
                "Task %d has no cube assigned (only %d cubes available). Exiting.",
                task0, len(cubes)
            )
            return
        cubes_this_task = [cubes[task0]]
        log.info("SLURM array task %d/%d → cube: %s", task0 + 1, n_tasks, cubes[task0][0])
    else:
        cubes_this_task = cubes
        log.info("Single job mode — processing all %d cubes sequentially", len(cubes))

    for region, SiteLabel in cubes_this_task:
        log.info("=" * 60)
        log.info("START  cube: region=%s  SiteLabel=%s", region, SiteLabel)
        try:
            saved = process_cube(region, SiteLabel)
            log.info("DONE   cube: %s → %d year CSV(s) saved", region, len(saved))
        except Exception as e:
            log.error("FAILED cube: %s — %s", region, e)
            import traceback
            traceback.print_exc()
            continue  # move to next cube even if this one failed

    log.info("=" * 60)
    log.info("All assigned cubes processed.")


if __name__ == "__main__":
    main()