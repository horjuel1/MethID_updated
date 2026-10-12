#!/bin/bash
#SBATCH --job-name=NF_deconvolution  
#SBATCH --nodes=1                       
#SBATCH --cpus-per-task=2              
#SBATCH --mem-per-cpu=130G               
#SBATCH --time=2-00:00:00                  
#SBATCH --partition=cancergen,shared    
#SBATCH --output=/dcs11/scharpf/data/horjuela/logs/nextflow_deconvolution_%j.out

#-----------------------------------------------------------
# User Parameters
#-----------------------------------------------------------
methID_map_path="/dcs11/scharpf/data/horjuela/horjuela/methID_maps/repeat_evolution_meta_final/fstat_map/fstat_map_95.tsv"
sample_sheet_path="/dcs11/scharpf/data/horjuela/3_2026_tumor_evolution/sample_sheet/6_16_2026_tumor_evo_adjusted.clean.tsv"
output_directory="/dcs11/scharpf/data/horjuela/3_2026_tumor_evolution/9_15_26_methid_final"
reference_genome="/dcl01/scharpf/data/pipeline-hub/pipeline-resources/genome/hg19/hg19.fa"

#-----------------------------------------------------------
# MethID Deconvolution Pipeline
#-----------------------------------------------------------

#Set up log directories in full output directory
mkdir -p "${output_directory}"
mkdir -p "${output_directory}/logs"
export NXF_LOG_FILE="${output_directory}/logs/nextflow.log"
export NXF_WORK_DIR="${output_directory}/nextflow_work"

# Load necessary modules for Nextflow and Singularity
module load nextflow/24.10.5
module load singularity/3.11.4

# Launch the Nextflow pipeline
nextflow run /dcs11/scharpf/data/horjuela/scripts/MethID-1.0.2/main.nf \
    --methid_map ${methID_map_path} \
    --sample_sheet ${sample_sheet_path} \
    --ref_genome ${reference_genome} \
    --outdir ${output_directory} \
    --min_cov 2 

