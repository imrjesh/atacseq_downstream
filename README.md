# ATAC-seq Downstream Analysis Pipeline

A modular, reproducible pipeline for downstream bulk ATAC-seq analysis, taking
a peak × sample count matrix through the full interpretive workflow:
quality control and filtering, differential accessibility testing (DARs),
unsupervised discovery of accessibility programs (NMF), peak-to-gene
annotation, functional enrichment, and transcription-factor activity inference
via two complementary and cross-validated methods — chromVAR motif
accessibility and TOBIAS footprinting. Config-driven, environment-pinned, and
runnable end-to-end on an included synthetic dataset.

## Two tracks

| Track | Input | Runs on toy data? | Outputs |
|-------|-------|-------------------|---------|
| Downstream | count matrix + peaks BED + metadata | yes | DARs, heatmap, NMF, annotation, GO, chromVAR |
| Footprinting | BAMs + genome FASTA + motifs | needs real BAMs | TOBIAS footprints, differential TF binding |

> Note: ATAC-seq measures **chromatin accessibility**, so the differential
> test yields **differentially accessible regions (DARs)**, not differential
> gene expression. Peak calling (MACS2) is an upstream step; a count matrix
> already implies peaks exist.

## Quickstart (downstream track, toy data)

```bash
conda env create -f environment.yml
conda activate atacseq-downstream

Rscript make_toy_data.R                                   # generate toy dataset
Rscript R/atacseq_downstream_pipeline.R config/config.yaml
```

Results land in `results/downstream/`. The toy data has 10% true DARs and two
latent programs, so you should recover significant DARs and see the two
programs in `nmf/`.

## Footprinting track (real data)

Requires duplicate-marked, filtered, per-condition merged BAMs, a genome FASTA,
and a JASPAR motif file. Fill the `footprinting:` block in `config/config.yaml`,
then:

```bash
snakemake -s workflow/footprinting.smk --use-conda --cores 8
```

## Integrating the two TF readouts

chromVAR (`motifs/differential_tf_activity.csv`) and TOBIAS
(`bindetect/bindetect_results.txt`) are independent estimates of TF activity.
TFs flagged by **both** are the highest-confidence regulatory drivers and top
the integrated prioritization.

## License

MIT
