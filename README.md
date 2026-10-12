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

*Building a reference atlas (build_atlas.nf)*

`build_atlas.nf` is the same atlas pipeline as smMethID's. Run it once per atlas:

```bash
nextflow run build_atlas.nf -params-file config/atlas_params.json
```

Copy `config/atlas_params.example.json` and `config/atlas_source.example.json` and edit both. Give the reference
WGBS beds either as `config` (a JSON of sources and cell types, read as `<base_dir>/<cell_type>/<bed_subdir>`) or as
`mapped_bed_dirs` (a comma-separated list of per-cell-type folders). It writes to `outdir`:

- `atlas.tsv`: replicates pooled (sum of methylated reads / sum of reads). Usable directly as `methid_map`.
- `atlas_beta_binomial_mean.tsv` / `_var.tsv` / `_dispersion.tsv`: the beta-binomial fit. The mean file is
  also usable as `methid_map`, and the var file is what `marker_ranking: standardized` needs.
- `atlas_pair_separability.csv` and `atlas_confusability_map*.png`: how confusable each pair of cell types is,
  on all CpGs and again on each `separability_marker_sets` panel (default `5000,20000`).
- `atlas_marker_balance.tsv`: per cell type, how many CpGs set it apart from its closest other cell type, and how
  strong its 5,000th / 20,000th / 50,000th best marker is. Use it to choose `markers_per_cell_type`.

*Marker selection*

Marker selection works the same way as in smMethID and runs at deconvolution time, so one atlas serves any panel
size. Set it in `params.json`:

| Field | Default | Description |
|---|---|---|
| `markers_per_cell_type` | `0` | `0` = every atlas CpG. `N` = for each cell type, the N CpGs where it differs most from its closest other cell type, pooled, so every cell type is represented equally. CpGs with a missing beta for any cell type are dropped first. |
| `marker_ranking` | `gap` | `gap`: difference to the closest other cell type. `standardized`: that difference divided by the combined replicate SD of the two cell types, so CpGs where replicates disagree rank low. Needs `var_atlas_file`; use an atlas built with `min_coverage` of 10 or more. |
| `var_atlas_file` | unset | `atlas_beta_binomial_var.tsv` from `build_atlas.nf`, for `standardized` ranking. |

```json
"methid_map": "results_atlas/atlas_beta_binomial_mean.tsv",
"markers_per_cell_type": 5000,
"marker_ranking": "standardized",
"var_atlas_file": "results_atlas/atlas_beta_binomial_var.tsv"
```

The selected CpGs are identical to smMethID's for the same atlas and settings. Each sample's log lists the
markers kept per cell type, and `atlas_length` in the coverage table becomes the panel size. MethID only uses
panel CpGs that a sample covers at least `min_cov` times, so check `overlap_sites`: smMethID's per-read model
gets by with fewer markers than this site-level fit does on low-coverage cfDNA.
