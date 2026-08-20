# =============================================================================
# workflow/footprinting.smk
#
# TOBIAS TF-footprinting module for ATAC-seq.
#
# WHERE THIS FITS:
#   This is an UPSTREAM/PARALLEL module. It operates on BAM files (base-pair
#   resolution Tn5 cut signal), NOT on the peak x sample count matrix used by
#   the downstream R pipeline. It complements chromVAR:
#       chromVAR  -> per-sample motif accessibility (from the count matrix)
#       TOBIAS    -> per-condition TF footprints / bound-vs-unbound (from BAMs)
#
# TOBIAS stages:
#   1. ATACorrect   - correct Tn5 insertion sequence bias in the BAM
#   2. ScoreBigwig  - compute continuous footprint scores across peaks
#   3. BINDetect    - per-motif differential binding between conditions
#
# Requirements (install via conda; see environment.yml):
#   tobias, samtools, bedtools
#
# Inputs (declared in config/config.yaml under a `footprinting:` block):
#   footprinting:
#     conditions:                       # one merged BAM per condition
#       tumor:  data/bam/tumor.merged.bam
#       normal: data/bam/normal.merged.bam
#     genome_fasta: refs/hg38.fa        # genome FASTA (indexed .fai alongside)
#     peaks:        data/peaks/consensus_peaks.bed   # consensus peak set (BED)
#     motifs:       refs/JASPAR2022_CORE_vertebrates.jaspar  # motif database
#     blacklist:    refs/hg38.blacklist.bed          # optional ENCODE blacklist
#     outdir:       results/footprinting
#     threads:      8
#
# NOTE on BAMs: TOBIAS expects duplicate-marked, properly filtered ATAC BAMs
# (mitochondrial reads removed, blacklist-filtered). Merge replicates per
# condition upstream (samtools merge) so each condition has one representative
# signal track. That merge step is included below.
#
# Usage:
#   snakemake -s workflow/footprinting.smk --use-conda --cores 8
# =============================================================================

import os

configfile: "config/config.yaml"
FP = config["footprinting"]

CONDITIONS  = list(FP["conditions"].keys())
GENOME      = FP["genome_fasta"]
PEAKS       = FP["peaks"]
MOTIFS      = FP["motifs"]
OUTDIR      = FP["outdir"]
THREADS     = FP.get("threads", 8)
BLACKLIST   = FP.get("blacklist", None)

# ---------------------------------------------------------------------------
# Target rule: BINDetect output is the final deliverable (differential binding)
# ---------------------------------------------------------------------------
rule all:
    input:
        os.path.join(OUTDIR, "bindetect", "bindetect_results.txt"),
        expand(os.path.join(OUTDIR, "footprints", "{cond}_footprints.bw"),
               cond=CONDITIONS)

# ---------------------------------------------------------------------------
# 0. (Optional) merge replicate BAMs per condition.
#    If your config already points at merged BAMs, this rule is a no-op copy.
# ---------------------------------------------------------------------------
rule merged_bam:
    input:
        lambda wc: FP["conditions"][wc.cond]
    output:
        bam = os.path.join(OUTDIR, "bam", "{cond}.bam"),
        bai = os.path.join(OUTDIR, "bam", "{cond}.bam.bai")
    threads: THREADS
    shell:
        # If input is already a single merged, indexed BAM this just links+indexes.
        "mkdir -p $(dirname {output.bam}) && "
        "cp {input} {output.bam} && "
        "samtools index -@ {threads} {output.bam}"

# ---------------------------------------------------------------------------
# 1. ATACorrect: remove Tn5 cut-site sequence bias. Produces a bias-corrected
#    signal bigwig that footprint scoring depends on.
# ---------------------------------------------------------------------------
rule atacorrect:
    input:
        bam    = os.path.join(OUTDIR, "bam", "{cond}.bam"),
        genome = GENOME,
        peaks  = PEAKS
    output:
        corrected = os.path.join(OUTDIR, "atacorrect",
                                 "{cond}_corrected.bw")
    params:
        outdir    = os.path.join(OUTDIR, "atacorrect"),
        blacklist = f"--blacklist {BLACKLIST}" if BLACKLIST else ""
    threads: THREADS
    shell:
        "TOBIAS ATACorrect "
        "--bam {input.bam} "
        "--genome {input.genome} "
        "--peaks {input.peaks} "
        "{params.blacklist} "
        "--outdir {params.outdir} "
        "--prefix {wildcards.cond} "
        "--cores {threads}"

# ---------------------------------------------------------------------------
# 2. ScoreBigwig: turn the corrected signal into continuous footprint scores
#    across the peak regions.
# ---------------------------------------------------------------------------
rule footprint_scores:
    input:
        signal = os.path.join(OUTDIR, "atacorrect", "{cond}_corrected.bw"),
        peaks  = PEAKS
    output:
        footprints = os.path.join(OUTDIR, "footprints", "{cond}_footprints.bw")
    threads: THREADS
    shell:
        "TOBIAS ScoreBigwig "
        "--signal {input.signal} "
        "--regions {input.peaks} "
        "--output {output.footprints} "
        "--cores {threads}"

# ---------------------------------------------------------------------------
# 3. BINDetect: for every motif, estimate bound vs unbound sites in each
#    condition and test for DIFFERENTIAL binding between conditions.
#    This is the payoff: a ranked table of TFs with changed footprint occupancy.
# ---------------------------------------------------------------------------
rule bindetect:
    input:
        footprints = expand(os.path.join(OUTDIR, "footprints",
                                         "{cond}_footprints.bw"),
                            cond=CONDITIONS),
        motifs = MOTIFS,
        genome = GENOME,
        peaks  = PEAKS
    output:
        results = os.path.join(OUTDIR, "bindetect", "bindetect_results.txt")
    params:
        outdir = os.path.join(OUTDIR, "bindetect"),
        conds  = " ".join(CONDITIONS)
    threads: THREADS
    shell:
        "TOBIAS BINDetect "
        "--motifs {input.motifs} "
        "--signals {input.footprints} "
        "--genome {input.genome} "
        "--peaks {input.peaks} "
        "--cond_names {params.conds} "
        "--outdir {params.outdir} "
        "--cores {threads}"

# ---------------------------------------------------------------------------
# INTERPRETATION (documented for repo readers):
#   bindetect_results.txt ranks each TF motif by a differential binding score
#   between conditions. Positive = more footprint occupancy in condition B,
#   negative = more in condition A. Cross-reference these TFs against the
#   chromVAR differential_tf_activity.csv from the downstream R pipeline:
#   TFs flagged by BOTH methods are your highest-confidence regulatory drivers,
#   and belong at the top of the integrated prioritization shortlist.
# ---------------------------------------------------------------------------
