*MethID Deconvolution*

MethID deconvolution is a method for estimating the cell contributions for a given cell-free DNA sample, using low coverage, whole genome methylation data. It uses a non-negative least squares algorithm of many positions of atleast 2x coverage for identifying the closest combination of cell-specific methylation profiles to match input methylation data. 

This Nextflow pipeline is intended to be used on a high performance computing cluster using the job scheduler SLURM. 

*SETUP*

## Setup

### Prerequisites
- **Local machine:** Docker (Desktop or Engine)
- **Compute environment:** Singularity or Apptainer

Docker and Singularity usually live on different machines — you build the image
locally, then move it to the cluster.

### 1. Clone the repo

```bash
git clone <REPO_URL>
cd <REPO_NAME>
```

### 2. Build the Docker image (local machine)

```bash
docker buildx build --platform linux/amd64 -t methid --load .
```

`--platform linux/amd64` is required. Prebuilt R package binaries are only
published for amd64, and the pipeline runs on x86 clusters. On Apple Silicon
this builds under emulation — expect it to be slow, and don't remove the flag.

### 3. Export the image to a tarball (local machine)

```bash
docker save -o methid.tar methid:latest
```

### 4. Transfer to your compute environment

```bash
# gzip is optional but the tarball is large
gzip methid.tar
scp methid.tar.gz user@cluster:/path/to/parent-of-repo/
```

### 5. Convert to a Singularity image (cluster)

```bash
gunzip methid.tar.gz     # docker-archive needs an uncompressed tar
singularity build methid.sif docker-archive://methid.tar
```

Use `apptainer` in place of `singularity` on newer systems — the commands are
otherwise identical. If the build fails on permissions, try `--fakeroot`.

### 6. Place the container

`methid.sif` must sit in the **parent directory of this repo**:

```
parent-directory/
├── methid.sif
└── <REPO_NAME>/
```

### 7. Verify

```bash
singularity exec methid.sif Rscript -e 'library(GenomicRanges); cat("ok\n")'
singularity exec methid.sif modkit --help > /dev/null && echo "modkit ok"
```

*To use this SLURM pipeline*

1. place your modified-base BAM files (and their `.bai` indexes) in an input folder. Each sample is named after its BAM filename; no sample sheet is required.
   Alternatively, list samples in a sample sheet (see below).

2. select a methID reference and reference genome. MethID reference is of format:
`chr/start/start+2/ref_cell_1_beta/ref_cell_2_beta/etc`
MethID reference must refer to the genome used to align .bam files. 
Included in the repo is the MethID reference used to generate figure 5C in Orjuela et al 2026. 

3. edit `params.json` to set the input paths:
   - `methid_map`  : path to the methID reference map (.tsv)
   - `input_dir`   : folder containing the input BAMs (+ `.bai` indexes)
   - `ref_genome`  : reference genome FASTA (must match the BAM alignment genome)
   - `outdir`      : output directory
   - `min_cov`     : minimum coverage for a site to be considered
   - `input_file_type` : `bam` (default), `bedmethyl`, `shortbed`, or `mod_pileup`
   - `mod_bases`   : BAM input only. `5mC_5hmC` (default) counts 5mC + 5hmC together as methylated,
                     matching WGBS-based maps; `5mC` counts 5mC alone
   - `sample_sheet`: optional, used instead of `input_dir` (see below)

*Sample sheet (optional)*

Instead of `input_dir`, you can pass a tab-separated `sample_sheet` with columns `id`, `label`, and `path`
(`bam_path` also works, so MethID 1.0.2 sheets can be reused). `path` must be absolute and point to a file
of the selected `input_file_type`. Samples are named `<id>_<label>`.

```
id	label	path
P001	baseline	/abs/path/P001.bam
P001	week4	/abs/path/P001_wk4.bam
```

Existing shortbeds (what the old `main_deconv_only.nf` did) can be deconvolved again with
`"input_file_type": "shortbed"` and `input_dir` set to `<old outdir>/shortbeds`.

*5-base / 5CmC data*

For Illumina 5-base or Watchmaker data processed with the 5CmC pipeline, start from Rastair's
`mod_pileup.txt` instead of a BAM:

```json
{
    "methid_map": "map.tsv",
    "input_file_type": "mod_pileup",
    "input_dir": "path/to/5CmC/output",
    "outdir": "results_dir",
    "min_cov": 2
}
```

- `input_dir` is the folder holding the per-sample 5CmC directories (the central
  `5CmC/<version>/<project>/` folder, or the folder you launched 5CmC from, which holds symlinks to them).
  Each `<sample>_..._hg19_wm_idt_spike_v1.0__5cmc-v1.1/mod_pileup.txt` is picked up, and its sample is named
  after the part before the CG-Align/5CmC suffix (`SampleABC`).
- Rows where the sample's genotype isn't the reference C (`C/C` on `+`, `G/G` on `-`) are removed,
  because C>T SNPs look like methylation calls. Then both strands are summed per CpG, the same way
  `modkit pileup --combine-strands` combines them.
- No modkit step and no `ref_genome` are needed. The CG-Align reference must be the same genome build as the MethID map.

4. submit the pipeline with `sbatch launch.sh`. Each sample is deconvolved in parallel; combined results are written to `<outdir>/deconvolution/deconvolved_samples.csv`. 
