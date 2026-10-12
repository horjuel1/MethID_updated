#!/bin/bash
#SBATCH --job-name=NF_deconvolution  
#SBATCH --nodes=1                       
#SBATCH --cpus-per-task=4              
#SBATCH --mem-per-cpu=8G               
#SBATCH --time=8:00:00                  
#SBATCH --partition=cancergen,shared    
#SBATCH --output=logs/LAUNCH_NF_deconv_o%j.txt

#-----------------------------------------------------------
# MethID Deconvolution Pipeline
#
# Input parameters are set in params.json (methid_map, input_dir,
# ref_genome, outdir, min_cov). Edit that file rather than this script.
#-----------------------------------------------------------

#Receive the first argument as the path to the params.json file
PARAMS_FILE=${1:-params.json}

#Set up log directories in full output directory
mkdir -p "logs"
export NXF_LOG_FILE="logs/nextflow.log"
export NXF_WORK_DIR="${MYSCRATCH}/nextflow_work"

# Load necessary modules for Nextflow and Singularity
module load nextflow/24.10.5
module load singularity/3.11.4

# Launch the Nextflow pipeline
nextflow run main.nf -params-file ${PARAMS_FILE}