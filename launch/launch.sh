#!/bin/bash
#SBATCH --job-name=MethID
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-00:00:00
#SBATCH --partition=cancergen,shared
#SBATCH --output=MethID_%j.log

#-----------------------------------------------------------
# MethID: deconvolution (main.nf)
#
# Run it from a run directory that holds this run's config:
#
#   cd /path/to/my_run                 # contains params.json (+ sample sheet)
#   sbatch /path/to/MethID/launch/launch.sh [params.json]
#
# Relative paths in params.json (methid_map, input_dir, sample_sheet, outdir, ...)
# resolve against the run directory, and results, logs and Nextflow's resume
# history land there too. Tasks run in the methid.sif container (see nextflow.config).
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
if [[ ! -f "$REPO/main.nf" ]]; then
  echo "ERROR: could not find the MethID repo (no main.nf in '$REPO'); set METHID_DIR=/path/to/MethID" >&2
  exit 1
fi

# First argument = params file (default params.json) unless it is a Nextflow option.
if [[ $# -gt 0 && "$1" != -* ]]; then PARAMS="$1"; shift; else PARAMS="params.json"; fi
if [[ ! -f "$PARAMS" ]]; then
  echo "ERROR: $PARAMS not found in $(pwd). Copy $REPO/config/params.json here and edit it." >&2
  exit 1
fi

mkdir -p logs
export NXF_LOG_FILE="logs/nextflow.log"

module load nextflow/24.10.5
module load singularity/3.11.4

# Anything after the params file is passed to Nextflow.
# -work-dir (not NXF_WORK_DIR, which Nextflow ignores) puts intermediate
# files on scratch, one folder per run directory.
nextflow run "$REPO/main.nf" -params-file "$PARAMS" \
    -work-dir "${MYSCRATCH:-$PWD}/nextflow_work/$(basename "$PWD")" \
    -resume \
    "$@"
