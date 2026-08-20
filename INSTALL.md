# Installation & Reproducibility Notes

## Tested environment

The downstream track was installed and run successfully on:

- **Machine:** Apple MacBook Pro, M2 Pro (Apple Silicon, arm64)
- **OS:** macOS
- **Package manager:** Miniforge (conda 26.x + mamba)
- **R:** installed via conda (r-base from conda-forge)
- **Date verified:** February 2026

Toy-data run completed end to end in ~3 minutes, recovering 137 differentially
accessible regions (71 gained, 66 lost) and 2 NMF accessibility programs — the
signal that the synthetic generator embeds, confirming the pipeline works.

## Install (downstream track)

```bash
# 1. Install Miniforge (Apple Silicon example)
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh      # accept license, allow conda init

# 2. Create the environment
mamba env create -f environment.yml
conda activate atacseq-downstream

# 3. Install the NMF R package separately (see note below)
R -e 'install.packages("NMF", repos="https://cloud.r-project.org")'
```

## Run

```bash
Rscript make_toy_data.R
Rscript R/atacseq_downstream_pipeline.R config/config.yaml
```

Outputs are written to `results/downstream/` (qc/, dars/, nmf/).

## Known gotchas (and fixes)

- **NMF is installed separately, not via conda.** The conda `r-nmf` package can
  conflict with the R version pulled in by the Bioconductor stack. Installing
  NMF from CRAN inside R after the env is built avoids the conflict.
- **Avoid hard version pins in `environment.yml`.** Pinning individual
  Bioconductor packages (e.g. `deseq2=1.42`) causes unsolvable dependency
  chains. Let conda/mamba resolve a mutually consistent R + Bioconductor set.
- **On toy data, annotation / GO / chromVAR steps skip themselves** because no
  genome is configured in `config.yaml`. This is expected; QC, DARs, heatmap,
  and NMF still run.

## Footprinting track (TOBIAS)

TOBIAS is kept in a **separate environment** because its dependencies can clash
with the large R stack. It also requires real BAM files, not the toy count
matrix. See `workflow/footprinting.smk` and the `footprinting:` block in
`config/config.yaml`.
