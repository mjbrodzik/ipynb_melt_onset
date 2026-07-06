#!/bin/bash
# ==============================================================================
# SLURM Array Job — Tb Threshold Computation
# ==============================================================================
#
# PURPOSE
# -------
# Submits one SLURM task per CETB cube for parallel Tb threshold computation.
# Each task runs run_Tbthresholds_informative.py, which auto-discovers all
# N25_* cubes and assigns exactly one cube to this task based on
# SLURM_ARRAY_TASK_ID.
#
# USAGE
# -----
# First, check how many cubes are available:
#     ls -d /pl/active/PMESDR/CETB_cubesv2/F16_SSMIS/N/nc_cubes/cubes_N25_* | wc -l
#
# Then update --array=0-(N-1) below and submit:
#     sbatch run_Tbthresholds.slurm
#
# Monitor progress:
#     squeue -u mabo7890
#     tail -f logs/tb_thresh_JOBID_0.out
#     tail -n 5 logs/tb_thresh_JOBID_*.err
#
# Check for failures:
#     grep "FAILED\|ERROR\|OOM" logs/tb_thresh_JOBID_*.err
#
# NOTES
# -----
# - Partition and QoS will be renamed on August 5, 2026:
#     amilan → acpu
#     normal → cpu-normal
#   Update those lines before that date.
#
# - 128G memory is required because the full TB array for a 640×640 cube
#   over 15 years of daily data is ~14GB. Other processing overhead
#   (snow class 4GB NetCDF, DataFrames) brings total to ~20-30GB.
#   128G gives comfortable headroom and avoids OOM kills.
#
# - cpus-per-task=2: the Python script is single-threaded (serial per-pixel loop).
#   2 CPUs allows numpy/xarray internal parallelism (OpenBLAS) for array ops.
#   More CPUs would not speed up the serial loop.
#
# - Each task writes its own CSV to ThresholdDir. No inter-task communication
#   or file locking is needed — tasks are fully independent.
#
# - If a task fails (OOM, timeout), resubmit this script unchanged.
#   The Python script skips cubes whose output CSV already exists.
#
# OUTPUT
# ------
# CSVs: /scratch/alpine/mabo7890/Tb_Thresholds_NewCubes/TbThresholds_*.csv
# Logs: /projects/mabo7890/Melt_Codes/scripts/logs/tb_thresh_%A_%a.{out,err}
#   %A = array job ID (same for all tasks)
#   %a = task index (0, 1, 2, ... one per cube)
# ==============================================================================

# ------------------ SLURM DIRECTIVES ------------------

#SBATCH --job-name=tb_thresh            # name shown in squeue

#SBATCH --partition=acpu                # CPU partition
                                        # (was 'amilan' — renamed Aug 5, 2026)
#SBATCH --qos=cpu-normal                # QoS for acpu partition
                                        # (was 'normal' — renamed Aug 5, 2026)

#SBATCH --time=8:00:00                  # wall time limit per task
                                        # 640x640 cube x 15 years ~ 4-6 hrs
                                        # 8 hrs gives safe headroom

#SBATCH --nodes=1                       # always 1 — script is single-node
#SBATCH --ntasks=1                      # always 1 — no MPI
#SBATCH --cpus-per-task=2              # 2 CPUs for numpy/xarray internal threads
#SBATCH --mem=128G                      # ~14GB TB array + 4GB snow class + overhead

#SBATCH --array=0-23                    # ← UPDATE THIS: 0-(N_cubes-1)
                                        # Check N with:
                                        # ls -d /pl/active/PMESDR/CETB_cubesv2/\
                                        #   F16_SSMIS/N/nc_cubes/cubes_N25_* | wc -l
                                        # Example: 24 cubes → --array=0-23

#SBATCH --output=/projects/mabo7890/Melt_Codes/scripts/logs/tb_thresh_%A_%a.out
#SBATCH --error=/projects/mabo7890/Melt_Codes/scripts/logs/tb_thresh_%A_%a.err
                                        # %A = array job ID, %a = task index
                                        # Each task gets its own log file

#SBATCH --mail-type=END,FAIL            # email only on completion or failure
                                        # (not on every task start — avoids spam)
#SBATCH --mail-user=your_email@colorado.edu


# ------------------ ENVIRONMENT SETUP ------------------

# Remove any modules loaded by default — ensures clean, reproducible environment
module purge

# Use the cetb conda environment Python directly (no 'conda activate' needed in SLURM)
PYTHON=/projects/mabo7890/software/anaconda/envs/cetb/bin/python

# PROJ and GDAL data directories — required by pyproj and rasterio.
# Without these, coordinate transformations fail with:
#   RuntimeError: PROJ: proj_create_from_database: Open of .../share/proj failed
export PROJ_LIB=/projects/mabo7890/software/anaconda/envs/cetb/share/proj
export GDAL_DATA=/projects/mabo7890/software/anaconda/envs/cetb/share/gdal

# Print job context for debugging (appears in .out log)
echo "======================================================"
echo "  Tb Threshold — SLURM Array Job"
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
$PYTHON run_Tbthreshold.py

EXIT_CODE=$?

echo "======================================================"
echo "  Task ${SLURM_ARRAY_TASK_ID} finished at $(date)"
echo "  Exit code: ${EXIT_CODE}"
echo "======================================================"

# Exit with the Python script's exit code so SLURM marks
# the task as FAILED if the script crashed
exit ${EXIT_CODE}
