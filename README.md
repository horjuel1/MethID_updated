*MethID Deconvolution*

MethID deconvolution is a method for estimating the cell contributions for a given cell-free DNA sample, using low coverage, whole genome methylation data. It uses a non-negative least squares algorithm for identifying the closest combination of cell-specific methylation profiles to match input methylation data. 

Update- 1.0.2: Update modkit from Modkit V0.4.4 to Modkit V0.6.1 and combines 5mC and 5hmC to match WGBS data-based map 

To use this pipeline:

1. create a sample sheet (.tsv), which must have columns 'id', 'label', and 'bam_path'. Absolute bam_path is necessary. 

2. select a methID map (from /dcs07/Methylation_Atlas_Dev/methID_maps) and reference genome. Genome must match the genome used to align .bam files

3. launch the `launch.sh` nextflow script. 
