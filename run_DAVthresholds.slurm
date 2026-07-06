#!/bin/bash
# ==============================================================================
# SLURM Array Job — DAV Threshold Computation
# ==============================================================================
#
# PURPOSE
# -------
# Submits one SLURM task per CETB cube for parallel DAV threshold computation.
# Each task runs run_DAVthresholds_informative.py, which auto-discovers all
# N25_* cubes and assigns exactly one cube to this task based on
# SLURM_ARRAY_TASK_ID.
#
# Unlike the Tb script (one CSV per cube), this script produces one CSV
# per year per cube. A 640×640 cube × 15 years = 15 output CSVs per task.
#
# USAGE
# -----
# First, check how many cubes need processing:
#     ls -d /pl/active/PMESDR/CETB_cubesv2/F16_SSMIS/N/nc_cubes/cubes_N25_* | wc -l
#
# Then update --array=0-(N-1) below and submit:
#     sbatch run_DAVthresholds.slurm
#
# Monitor progress:
#     squeue -u mabo7890
#     tail -f logs/dav_thresh_JOBID_0.err
#     tail -n 5 logs/dav_thresh_JOBID_*.err
#
# Check for failures:
#     grep "FAILED\|ERROR" logs/dav_thresh_JOBID_*.err
#
# Check how many year CSVs were produced per cube:
#     ls /scratch/alpine/mabo7890/DAV_Thresholds/DAVThresholds_*_N25_44*.csv | wc -l
#
# RERUN SAFETY
# ------------
# Safe to resubmit unchanged after partial failures. The Python script uses
# two levels of skip logic:
#   - Cube level : skipped if ALL year CSVs already exist
#   - Year level : individual year CSVs skipped if they already exist
# Only missing year/cube combinations are reprocessed.
#
# NOTES
# -----
# - Partition and QoS will be renamed on August 5, 2026:
#     amilan → acpu
#     normal → cpu-normal
#   Update those lines before that date.
#
# - 128G memory is required because:
#     TB array  : 640×640 × 15yr × daily obs ≈ 14GB
#     DAV array : same shape as TB ≈ 14GB (computed from TB, held in memory)
#     Snow class: 4GB global NetCDF (cached once)
#     Total     : ~35-40GB peak; 128G gives comfortable headroom
#
# - cpus-per-task=2: the Python script is single-threaded (serial per-pixel
#   loop per year). 2 CPUs allows numpy internal parallelism for the
#   vectorized calc_DAV() call over the full (time, ny, nx) array.
#
# - Each task writes its own set of CSVs. No inter-task communication or
#   file locking needed — tasks are fully independent.
#
# - Time limit of 8 hours: DAV computation is faster than Tb (no histogram
#   fitting), but there are 15 year CSVs per cube vs 1 for Tb, so total
#   I/O is comparable. 8 hours gives safe headroom.
#
# OUTPUT
# ------
# CSVs: /scratch/alpine/mabo7890/DAV_Thresholds/DAVThresholds_{year}_*.csv
#       (one file per year per cube — 15 files per task for a 15-year record)
# Logs: /projects/mabo7890/Melt_Codes/scripts/logs/dav_thresh_%A_%a.{out,err}
#   %A = array job ID (same for all tasks in this submission)
#   %a = task index (0, 1, 2, ... one per cube)
# ==============================================================================

# ------------------ SLURM DIRECTIVES ------------------

#SBATCH --job-name=dav_thresh           # name shown in squeue

#SBATCH --partition=acpu                # CPU partition
                                        # (was 'amilan' — renamed Aug 5, 2026)
#SBATCH --qos=cpu-normal                # QoS for acpu partition
                                        # (was 'normal' — renamed Aug 5, 2026)

#SBATCH --time=8:00:00                  # wall time limit per task
                                        # DAV loop over 15 years ~ 4-6 hrs
                                        # 8 hrs gives safe headroom

#SBATCH --nodes=1                       # always 1 — script is single-node
#SBATCH --ntasks=1                      # always 1 — no MPI
#SBATCH --cpus-per-task=2              # 2 CPUs for numpy vectorized DAV compute
#SBATCH --mem=128G                      # TB (~14GB) + DAV (~14GB) + snow class
                                        # (4GB) + overhead = ~35-40GB peak;
                                        # 128G gives comfortable headroom

#SBATCH --array=0-23                    # ← UPDATE THIS: 0-(N_cubes-1)
                                        # Check N with:
                                        # ls -d /pl/active/PMESDR/CETB_cubesv2/\
                                        #   F16_SSMIS/N/nc_cubes/cubes_N25_* | wc -l
                                        # Example: 24 cubes → --array=0-23

#SBATCH --output=/projects/mabo7890/Melt_Codes/scripts/logs/dav_thresh_%A_%a.out
#SBATCH --error=/projects/mabo7890/Melt_Codes/scripts/logs/dav_thresh_%A_%a.err
                                        # %A = array job ID, %a = task index
                                        # Each task gets its own log file pair

#SBATCH --mail-type=END,FAIL            # email only on completion or failure
                                        # (not ALL — avoids one email per task)
#SBATCH --mail-user=your_email@colorado.edu


# ------------------ ENVIRONMENT SETUP ------------------

# Remove any modules loaded by default — ensures clean, reproducible environment
module purge

# Use the cetb conda environment Python directly (no 'conda activate' in SLURM)
PYTHON=/projects/mabo7890/software/anaconda/envs/cetb/bin/python

# PROJ and GDAL data directories — required by pyproj and rasterio.
# Without these, coordinate transformations fail with:
#   RuntimeError: PROJ: proj_create_from_database: Open of .../share/proj failed
export PROJ_LIB=/projects/mabo7890/software/anaconda/envs/cetb/share/proj
export GDAL_DATA=/projects/mabo7890/software/anaconda/envs/cetb/share/gdal

# Print job context for debugging (appears in .out log)
echo "======================================================"
echo "  DAV Threshold — SLURM Array Job"
echo "======================================================"
echo "  Array Job ID  : ${SLURM_ARRAY_JOB_ID}"
echo "  Task ID       : ${SLURM_ARRAY_TASK_ID}"
echo "  Node          : $(hostname)"
echo "  Started       : $(date)"
echo "  Python        : ${PYTHON}"
echo "  PROJ_LIB      : ${PROJ_LIB}"
echo "======================================================"


# ------------------ RUN SCRIPT ------------------

cd /projects/mabo7890/Melt_Codes/scripts

# The Python script reads SLURM_ARRAY_TASK_ID to determine which cube
# to process. No arguments needed — cube assignment is automatic.
# Output: 15 per-year CSVs written to DAV_Thresholds/ for this cube.
$PYTHON run_DAVthreshold.py

EXIT_CODE=$?

echo "======================================================"
echo "  Task ${SLURM_ARRAY_TASK_ID} finished at $(date)"
echo "  Exit code: ${EXIT_CODE}"
echo "======================================================"

# Exit with the Python script's exit code so SLURM marks
# the task as FAILED if the script crashed (not silently SUCCESS)
exit ${EXIT_CODE}
