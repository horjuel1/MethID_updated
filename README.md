# MethID

MethID estimates the cell-type contributions to a cell-free DNA sample from low-coverage, whole-genome
methylation data. It fits a non-negative least squares model of reference cell-type methylation profiles to the
sample's methylation at every CpG covered at least `min_cov` times.

It runs as two Nextflow pipelines on a SLURM cluster, each with one params file and one launch script:

| Step | Pipeline | Params template | Launch script |
|---|---|---|---|
| 1. Build a reference map (once per map) | `build_atlas.nf` | `config/atlas_params.json` | `launch/launch_build_atlas.sh` |
| 2. Deconvolve samples | `main.nf` | `config/params.json` | `launch/launch.sh` |

Step 1 is optional. `Orjuela2026.tsv.gz`, the hg19 map used for figure 5C in Orjuela et al. 2026, works as is.

```mermaid
flowchart TD
    subgraph BUILD["Step 1: build_atlas.nf"]
        direction TB
        REF[/"Reference WGBS beds<br/>sources: base_dir / cell_type / bed_hg19"/]
        POOL["BUILD_ATLAS<br/>pool replicates per cell type"]
        BB["ATLAS_BETA_BINOMIAL_FIT<br/>mean + variance per CpG"]
        SEP["ATLAS_SEPARABILITY<br/>how confusable each pair of cell types is"]
        ATLAS[("atlas.tsv")]
        MEAN[("atlas_beta_binomial_mean.tsv")]
        VAR[("atlas_beta_binomial_var.tsv")]
        QC[("separability heatmaps<br/>atlas_marker_balance.tsv")]
        REF --> POOL --> ATLAS --> BB
        BB --> MEAN
        BB --> VAR
        MEAN --> SEP --> QC
    end

    subgraph DECON["Step 2: main.nf"]
        direction TB
        BAM[/"platform: ont<br/>modBAMs"/]
        BEDM[/"start_from: bedmethyl<br/>modkit output"/]
        PILE[/"platform: 5base<br/>5CmC mod_pileup.txt"/]
        SHORTIN[/"start_from: shortbed<br/>earlier run"/]
        MODKIT["MODKIT_PILEUP<br/>mod_bases: 5mC_5hmC or 5mC"]
        B2S["BEDMETHYL_TO_SHORTBED"]
        P2S["MOD_PILEUP_TO_SHORTBED<br/>drop SNPs, merge strands"]
        SHORT[("per-sample shortbed<br/>chr start end n_meth n_cov")]
        DEC["DECONVOLVE, one task per sample<br/>CpGs with at least min_cov reads<br/>optional balanced marker selection<br/>non-negative least squares"]
        COMB["COMBINE"]
        OUT[("deconvolved_samples.csv<br/>coverage table")]
        BAM --> MODKIT --> B2S
        BEDM --> B2S
        B2S --> SHORT
        PILE --> P2S --> SHORT
        SHORTIN --> SHORT
        SHORT --> DEC --> COMB --> OUT
    end

    ORJ[("Orjuela2026.tsv.gz<br/>bundled hg19 map")]
    ATLAS -. "methid_map" .-> DEC
    MEAN -. "methid_map" .-> DEC
    ORJ -. "methid_map" .-> DEC
    VAR -. "var_atlas_file<br/>standardized ranking" .-> DEC
    QC -. "choose<br/>markers_per_cell_type" .-> DEC
```

```
MethID/
├── main.nf                          # pipeline: deconvolve samples
├── build_atlas.nf                   # pipeline: build a reference map
├── nextflow.config                  # SLURM resources and the container
├── Dockerfile                       # container recipe
├── Orjuela2026.tsv.gz               # bundled hg19 reference map
├── bin/                             # R scripts the pipelines call
├── config/
│   ├── params.json                  # template: settings for main.nf
│   ├── atlas_params.json            # template: settings for build_atlas.nf
│   ├── sample_sheet.example.csv     # template: optional sample sheet
│   └── metadata_sheet.example.csv   # template: optional metadata for the summaries
└── launch/
    ├── launch.sh                    # sbatch launcher for main.nf
    └── launch_build_atlas.sh        # sbatch launcher for build_atlas.nf
```

The repo itself stays untouched. Each run lives in its own directory holding its params file (copied from
`config/`); you launch from there, and its results, logs and Nextflow resume history land there too.

## Setup

### 1. Clone the repo on the cluster

```bash
git clone git@github.com:horjuel1/MethID_updated.git
```

### 2. Build the container (on a local machine with Docker)

The container holds R 4.4 with every package the pipelines use, samtools and modkit 0.6.1.

```bash
docker buildx build --platform linux/amd64 -t methid --load .
```

`--platform linux/amd64` is required: the cluster is x86 and the R package binaries are amd64 only. On Apple
Silicon the build runs under emulation and is slower.

Then convert it to a Singularity image. Either convert on your machine using Singularity in Docker (the same
version as the cluster's module):

```bash
docker save methid:latest -o methid.tar
docker run --rm --privileged --platform linux/amd64 -v "$PWD:/work" -w /work \
    quay.io/singularity/singularity:v3.11.4 build methid.sif docker-archive:///work/methid.tar
```

or copy `methid.tar` to the cluster and run `singularity build methid.sif docker-archive://methid.tar` there
(add `--fakeroot` if it fails on permissions).

### 3. Place the container

`nextflow.config` expects `methid.sif` in the folder above the repo. To keep it elsewhere, change `container`
in `nextflow.config` to its full path.

```
parent-directory/
├── methid.sif
└── MethID_updated/
```

Check it:

```bash
singularity exec methid.sif modkit --version
singularity exec methid.sif Rscript -e 'library(GenomicRanges); library(limSolve); cat("ok\n")'
```

Rebuild the container only when `Dockerfile` or `bin/requirements.r` change. The pipelines run the scripts in
`bin/` from the repo, not from inside the container.

## Running

Build a map (once per map; optional):

```bash
mkdir my_atlas && cd my_atlas
cp /path/to/MethID_updated/config/atlas_params.json .      # edit it
sbatch /path/to/MethID_updated/launch/launch_build_atlas.sh
```

Deconvolve a batch of samples:

```bash
mkdir my_run && cd my_run
cp /path/to/MethID_updated/config/params.json .            # edit it (+ a sample sheet, if used)
sbatch /path/to/MethID_updated/launch/launch.sh
```

Relative paths in the params file (`methid_map`, `input_dir`, `sample_sheet`, `outdir`, `base_dir`) resolve
against the run directory. The SLURM log (`MethID_<jobid>.log`) and `logs/nextflow.log` land there, and
Nextflow's work files go to `$MYSCRATCH/nextflow_work/<run directory name>`. Runs resume by default, so
resubmitting picks up where a failed run stopped.

A launcher takes another params file as its first argument and passes anything after it to Nextflow:

```bash
sbatch /path/to/MethID_updated/launch/launch.sh params_other.json -with-report report.html
```

The launchers find the repo from the path you submit them with; set `METHID_DIR=/path/to/MethID_updated` if
that ever fails.

## 1. Building a reference map (`config/atlas_params.json`)

`build_atlas.nf` is the same atlas pipeline as smMethID's. It pools the reference WGBS replicates per cell type,
fits a beta-binomial model, and measures how separable the cell types are.

| Field | Default | Description |
|---|---|---|
| `sources` | — | The reference beds: `{ "<source>": { "base_dir": "...", "cell_types": [...] } }`. Each cell type is read from `<base_dir>/<cell_type>/<bed_subdir>/*.bed` and named `<source>_<cell_type>`. |
| `bed_subdir` | `bed_hg19` | Folder under each cell type holding its beds (`chr start end meth cov`). |
| `samples_to_exclude` | `false` | Comma-separated bed file names to leave out. |
| `outdir` | `results_atlas` | Output folder. |
| `min_coverage` | `10` | Minimum pooled reads for a CpG to be reported. Keep it at 10 or more for `standardized` marker ranking. |
| `cpus` | `4` | CPUs for pooling. |
| `bb_min_df_for_dispersion`, `bb_min_phi`, `bb_max_phi` | `30`, `0.0001`, `0.999` | Beta-binomial fit settings. |
| `separability_marker_sets` | `5000,20000` | Panel sizes (markers per cell type) to also compute separability for; `""` to skip. |

Instead of `sources`, you can give `mapped_bed_dirs` (a comma-separated list of per-cell-type bed folders) or
`config` (a separate JSON file in the `sources` format). Use only one of the three.

Outputs in `outdir`:

- `atlas.tsv`: replicates pooled (sum of methylated reads / sum of reads), one column per cell type. Usable as `methid_map`.
- `atlas_beta_binomial_mean.tsv`: the beta-binomial mean per CpG. Also usable as `methid_map`.
- `atlas_beta_binomial_var.tsv`: the per-CpG variance, for `marker_ranking: standardized`.
- `atlas_pair_separability.csv` and `atlas_confusability_map*.png`: how hard each pair of cell types is to tell
  apart, on all CpGs and on each `separability_marker_sets` panel.
- `atlas_marker_balance.tsv`: per cell type, how many CpGs set it apart from its closest other cell type, and how
  strong its 5,000th / 20,000th / 50,000th best marker is. Use it to choose `markers_per_cell_type`.
- `atlas_sources.json`: the `sources` used, for the record.

## 2. Deconvolving samples (`config/params.json`)

| Field | Default | Description |
|---|---|---|
| `methid_map` | — | Reference map: `chr start end <cell type betas...>`, betas from 0 to 1, same genome build as the samples. |
| `platform` | — | `ont` (ONT modBAMs) or `5base` (5CmC output). Required. |
| `start_from` | `false` | Restart from an intermediate instead of the platform's raw input: `shortbed` (either platform) or `bedmethyl` (`ont` only). |
| `input_dir` | — | Folder with the inputs. Use this or `sample_sheet`. |
| `sample_sheet` | `false` | CSV listing the inputs (see below). Use this or `input_dir`. |
| `metadata_sheet` | `false` | Optional CSV adding `label` and `exposure` to each sample; turns on the summaries (see below). |
| `ref_genome` | — | Reference FASTA. `ont` only. |
| `mod_bases` | `5mC_5hmC` | `ont` only. `5mC_5hmC` counts 5mC and 5hmC together as methylated, matching WGBS-based maps; `5mC` counts 5mC alone. |
| `markers_per_cell_type` | `0` | Marker selection (see below). `0` uses every CpG in the map. |
| `marker_ranking` | `gap` | `gap` or `standardized`. |
| `var_atlas_file` | `false` | `atlas_beta_binomial_var.tsv`, needed for `standardized`. |
| `outdir` | `results` | Output folder. |
| `min_cov` | `2` | Minimum reads at a CpG for it to be used. |

### Inputs

| `platform` | `start_from` | `input_dir` holds | Sample name |
|---|---|---|---|
| `ont` | `false` | `*.bam` with `.bai` indexes (MM/ML tags) | BAM file name |
| `5base` | `false` | 5CmC sample folders, each with `mod_pileup.txt` | folder name before the CG-Align suffix |
| `ont` | `bedmethyl` | `*.bed` from `modkit pileup` | file name |
| either | `shortbed` | `*_short.bed` from an earlier MethID run (`<outdir>/shortbeds`) | file name before `_short` |

**5-base / 5CmC (`platform: 5base`).** Point `input_dir` at the 5CmC project folder, for example
`/dcs10/scharpf/data/5b-project-hub/2026-03_5-Base_Healthy_Ref`. Every
`<sample>_bt2-bwa_umi_true_hg19_wm_idt_spike_v0.1/mod_pileup.txt` is picked up and named `<sample>`
(`CGPLH339P_5B_lib2`). Positions where the sample's genotype isn't the reference C (`C/C` on `+`, `G/G` on `-`)
are dropped, since C>T SNPs look like methylation. The two strands are then summed per CpG, as
`modkit pileup --combine-strands` does. No `ref_genome` is needed.

**Sample sheet.** Same format as smMethID: a CSV (tab-separated if the file ends in `.tsv` or `.txt`) with a
`CGID` column, which names each sample, and one path column that depends on the input:

| Input | Path column | Points to |
|---|---|---|
| `platform: 5base` | `sample_dir` | the 5CmC sample folder (holding `mod_pileup.txt`) |
| `platform: ont` | `bam_path` | the BAM (with its `.bai` next to it) |
| `start_from: shortbed` / `bedmethyl` | `path` | the shortbed or bedmethyl file |

```
CGID,sample_dir
CGPLH339P,/dcs10/scharpf/data/5b-project-hub/2026-03_5-Base_Healthy_Ref/CGPLH339P_5B_lib2_bt2-bwa_umi_true_hg19_wm_idt_spike_v0.1
```

### Marker selection

This works the same way as in smMethID and runs at deconvolution time, so one map serves any panel size. With
`markers_per_cell_type: N`, for each cell type MethID keeps the N CpGs where it differs most from its closest other
cell type, and pools them, so every cell type is represented equally. CpGs with a missing beta for any cell type
are dropped first.

- `gap` ranks by the difference in beta to the closest other cell type.
- `standardized` divides that difference by the combined replicate SD of the two cell types, so CpGs where
  replicates disagree rank low. It needs `var_atlas_file` and a map built with `min_coverage` of 10 or more.

```json
"methid_map": "/path/to/results_atlas/atlas_beta_binomial_mean.tsv",
"markers_per_cell_type": 5000,
"marker_ranking": "standardized",
"var_atlas_file": "/path/to/results_atlas/atlas_beta_binomial_var.tsv"
```

For the same map and settings, the selected CpGs are identical to smMethID's. MethID only uses panel CpGs a sample
covers at least `min_cov` times, so with low-coverage cfDNA use larger panels than smMethID needs, and check
`overlap_sites` in the coverage table.

### Summaries by label and exposure (optional)

Set `metadata_sheet` to a CSV (tab-separated if named `.tsv` or `.txt`) with the sample sheet's columns plus
`label` and `exposure`. You can use your sample sheet with the two columns added:

```
CGID,sample_dir,label,exposure
CGPLH339P,/dcs10/.../CGPLH339P_5B_lib2_bt2-bwa_umi_true_hg19_wm_idt_spike_v0.1,healthy,unexposed
```

Samples are matched by `CGID` to the sample names in `deconvolved_samples.csv`. With a sample sheet those are the
`CGID`s; with `input_dir` they are the names derived from the files, such as `CGPLH339P_5B_lib2`. Empty `label`
or `exposure` values are grouped as `NA`. The columns are checked when the run starts. Samples that don't match
are left out of the summaries and listed in `summary_unmatched_samples.csv`.

Written to `<outdir>/summary/`:

- `summary_by_label.csv`: per label and cell type, the number of samples and the mean, SD, median, 25th and 75th
  percentiles, minimum and maximum percentage.
- `summary_by_label_exposure.csv`: the same per label, exposure and cell type.
- `summary_sample_counts.csv`: samples per label and exposure.
- `deconvolved_samples_with_metadata.csv`: every sample's percentages with its label and exposure, one row per
  sample and cell type, for your own plots and tests.
- `summary_unmatched_samples.csv`: samples that were deconvolved but not in the metadata sheet, or the reverse.

### Outputs

- `<outdir>/deconvolution/deconvolved_samples.csv`: cell-type percentages, one column per sample.
- `<outdir>/deconvolution/coverage_<date>.tsv`: per sample, CpGs in the input (`sample_sites`), CpGs used
  (`overlap_sites`) and the map or panel size (`atlas_length`).
- `<outdir>/deconvolution/per_sample/`: the same, one file per sample.
- `<outdir>/shortbeds/`: per-sample methylated and total read counts per CpG, reusable with `start_from: shortbed`.
- `<outdir>/bedmethyl/`: modkit output (`ont` only).
- `<outdir>/summary/`: summaries by label and exposure (with `metadata_sheet` only).
