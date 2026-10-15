#!/bin/bash
#SBATCH --job-name=methid_build_map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=3-00:00:00
#SBATCH --partition=cancergen,shared
#SBATCH --output=logs/build_map_%j.txt

#-----------------------------------------------------------
# MethID reference map building (build_atlas.nf)
#
# Settings live in params_build_map.json; edit that file, not this one.
# Submit from the repo directory:
#   sbatch launch_build_map.sh [params_file] [extra nextflow args, e.g. -resume]
#-----------------------------------------------------------

PARAMS_FILE=${1:-params_build_map.json}
shift $(( $# > 0 ? 1 : 0 ))

REPO_DIR=${METHID_DIR:-$SLURM_SUBMIT_DIR}
if [[ ! -f "${REPO_DIR}/build_atlas.nf" ]]; then
    echo "build_atlas.nf not found in ${REPO_DIR}. Submit from the MethID repo directory or set METHID_DIR." >&2
    exit 1
fi

mkdir -p logs
export NXF_LOG_FILE="logs/nextflow_build_map.log"
export NXF_WORK="${MYSCRATCH}/nextflow_work/methid_build_map"

module load nextflow/24.10.5
module load singularity/3.11.4

nextflow run "${REPO_DIR}/build_atlas.nf" -params-file "${PARAMS_FILE}" "$@"
