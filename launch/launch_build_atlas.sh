#!/bin/bash
#SBATCH --job-name=MethID_atlas
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=3-00:00:00
#SBATCH --partition=cancergen,shared
#SBATCH --output=MethID_atlas_%j.log

#-----------------------------------------------------------
# MethID reference map build (build_atlas.nf)
#
# Run it from an atlas directory that holds this atlas's config:
#
#   cd /path/to/my_atlas               # contains atlas_params.json
#   sbatch /path/to/MethID/launch/launch_build_atlas.sh [atlas_params.json]
#
# Relative paths in atlas_params.json (outdir, ...) resolve against the atlas
# directory, and the atlas files, logs and Nextflow's resume history land there.
#-----------------------------------------------------------

# Run directory = where sbatch was called.
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1

# Repo = where this script lives. Under sbatch the script runs from SLURM's
# spool copy, so ask SLURM for the path it was submitted with. Override with
# METHID_DIR if needed.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  script_path=$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -1)
else
  script_path="${BASH_SOURCE[0]}"
fi
REPO="${METHID_DIR:-$(cd "$(dirname "$script_path")/.." && pwd)}"
if [[ ! -f "$REPO/build_atlas.nf" ]]; then
  echo "ERROR: could not find the MethID repo (no build_atlas.nf in '$REPO'); set METHID_DIR=/path/to/MethID" >&2
  exit 1
fi

# First argument = params file (default atlas_params.json) unless it is a Nextflow option.
if [[ $# -gt 0 && "$1" != -* ]]; then PARAMS="$1"; shift; else PARAMS="atlas_params.json"; fi
if [[ ! -f "$PARAMS" ]]; then
  echo "ERROR: $PARAMS not found in $(pwd). Copy $REPO/config/atlas_params.json here and edit it." >&2
  exit 1
fi

mkdir -p logs
export NXF_LOG_FILE="logs/nextflow.log"

module load nextflow/24.10.5
module load singularity/3.11.4

# Anything after the params file is passed to Nextflow.
# -work-dir (not NXF_WORK_DIR, which Nextflow ignores) puts intermediate
# files on scratch, one folder per run directory.
nextflow run "$REPO/build_atlas.nf" -params-file "$PARAMS" \
    -work-dir "${MYSCRATCH:-$PWD}/nextflow_work/$(basename "$PWD")" \
    -resume \
    "$@"
