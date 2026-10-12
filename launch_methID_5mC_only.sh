#!/bin/bash
#SBATCH --job-name=NF_deconvolution  
#SBATCH --nodes=1                       
#SBATCH --cpus-per-task=4              
#SBATCH --mem-per-cpu=8G               
#SBATCH --time=8:00:00                  
#SBATCH --partition=cancergen,shared    
#SBATCH --output=/dcs07/scharpf/data/horjuela/dcs10_horjuela/Methylation_Atlas_Dev/publications/revisions_Orjuela2026/logs/LAUNCH_NF_deconv_o%j.txt

#-----------------------------------------------------------
# User Parameters
#-----------------------------------------------------------
methID_map_path="/dcs07/scharpf/data/horjuela/dcs10_horjuela/Methylation_Atlas_Dev/publications/Orjuela2026/data/methid_references/liver_blood.tsv"
sample_sheet_path="/dcs07/scharpf/data/horjuela/dcs10_horjuela/Methylation_Atlas_Dev/publications/Orjuela2026/data/MethID_Full_Sample_Sheet.tsv"
output_directory="/dcs07/scharpf/data/horjuela/dcs10_horjuela/Methylation_Atlas_Dev/publications/revisions_Orjuela2026/data/5mC_only"
reference_genome="/dcl01/scharpf/data/pipeline-hub/pipeline-resources/genome/hg19/hg19.fa"

#-----------------------------------------------------------
# MethID Deconvolution Pipeline
#-----------------------------------------------------------

#Set up log directories in full output directory
mkdir -p "${output_directory}"
mkdir -p "${output_directory}/logs"
export NXF_LOG_FILE="${output_directory}/logs/nextflow.log"
export NXF_WORK_DIR="${MYSCRATCH}/nextflow_work"

# Load necessary modules for Nextflow and Singularity
module load nextflow/24.10.5
module load singularity/3.11.4

# Launch the Nextflow pipeline
nextflow run /dcs07/scharpf/data/horjuela/dcs10_horjuela/Methylation_Atlas_Dev/publications/revisions_Orjuela2026/code/MethID-1.0.2/main_5mC.nf \
    --methid_map ${methID_map_path} \
    --sample_sheet ${sample_sheet_path} \
    --ref_genome ${reference_genome} \
    --outdir ${output_directory} 

