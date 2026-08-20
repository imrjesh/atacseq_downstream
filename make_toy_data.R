#!/usr/bin/env Rscript
# =============================================================================
# make_toy_data.R
#
# Generates a small, realistic SYNTHETIC ATAC-seq dataset so the downstream
# pipeline runs end-to-end on `git clone` with no external data.
#
# Produces (into data/example/):
#   counts.tsv     peak x sample raw count matrix (rownames = chr:start-end)
#   peaks.bed      matching peak coordinates (BED, 0-based half-open)
#   metadata.csv   sample sheet with `sample` and `condition` columns
#
# Design: negative-binomial counts (realistic ATAC overdispersion), a subset of
# peaks made truly differential between the two conditions, plus two latent
# "accessibility programs" so NMF has structure to recover.
#
# Usage:
#   Rscript make_toy_data.R            # defaults below
#   Rscript make_toy_data.R 3000 12    # 3000 peaks, 12 samples (6 vs 6)
# =============================================================================

set.seed(42)  # reproducible toy data

args      <- commandArgs(trailingOnly = TRUE)
n_peaks   <- ifelse(length(args) >= 1, as.integer(args[1]), 2000L)
n_samples <- ifelse(length(args) >= 2, as.integer(args[2]), 8L)
if (n_samples %% 2 != 0) stop("Use an even number of samples (split evenly across 2 conditions).")

outdir <- file.path("data", "example")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- sample metadata: half control, half treatment --------------------------
half   <- n_samples / 2
sample <- sprintf("S%02d", seq_len(n_samples))
cond   <- rep(c("control", "treatment"), each = half)
meta   <- data.frame(sample = sample, condition = cond, stringsAsFactors = FALSE)

# ---- peak coordinates on a single toy chromosome ----------------------------
peak_w  <- 500L
gap     <- 1500L
starts  <- seq(10000L, by = peak_w + gap, length.out = n_peaks)
ends    <- starts + peak_w
peak_id <- sprintf("chr1:%d-%d", starts, ends)

# ---- baseline accessibility per peak (log-normal spread of mean depth) -------
base_mu <- rlnorm(n_peaks, meanlog = 3.0, sdlog = 1.0)  # ~ tens to hundreds

# ---- inject true differential peaks -----------------------------------------
n_diff       <- round(0.10 * n_peaks)            # 10% of peaks are DARs
diff_idx     <- sample(seq_len(n_peaks), n_diff)
gained_idx   <- diff_idx[seq_len(n_diff %/% 2)]  # up in treatment
lost_idx     <- setdiff(diff_idx, gained_idx)    # down in treatment
fc_gained    <- runif(length(gained_idx), 2.0, 5.0)
fc_lost      <- runif(length(lost_idx),   2.0, 5.0)

# ---- two latent programs so NMF recovers structure --------------------------
# Program peaks get correlated boosts within a subset of samples.
prog1_peaks  <- sample(seq_len(n_peaks), round(0.05 * n_peaks))
prog2_peaks  <- sample(setdiff(seq_len(n_peaks), prog1_peaks), round(0.05 * n_peaks))
prog1_samps  <- sample(seq_len(n_samples), half)          # arbitrary sample subset
prog2_samps  <- setdiff(seq_len(n_samples), prog1_samps)

# ---- simulate counts (negative binomial; size controls overdispersion) ------
size_disp <- 8  # smaller = more overdispersed
counts <- matrix(0L, nrow = n_peaks, ncol = n_samples,
                 dimnames = list(peak_id, sample))

for (j in seq_len(n_samples)) {
  mu <- base_mu
  is_treat <- meta$condition[j] == "treatment"
  if (is_treat) {
    mu[gained_idx] <- mu[gained_idx] * fc_gained
    mu[lost_idx]   <- mu[lost_idx]   / fc_lost
  }
  if (j %in% prog1_samps) mu[prog1_peaks] <- mu[prog1_peaks] * 2.0
  if (j %in% prog2_samps) mu[prog2_peaks] <- mu[prog2_peaks] * 2.0

  # per-sample library-size factor (realistic depth variation)
  lib <- runif(1, 0.7, 1.3)
  counts[, j] <- rnbinom(n_peaks, mu = mu * lib, size = size_disp)
}

# ---- write outputs ----------------------------------------------------------
# count matrix (tab-delimited, peak rownames in first column)
write.table(data.frame(peak = rownames(counts), counts, check.names = FALSE),
            file = file.path(outdir, "counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# BED (0-based start to match the coordinate convention read by the pipeline)
bed <- data.frame(chr = "chr1", start = starts - 1L, end = ends)
write.table(bed, file = file.path(outdir, "peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# metadata
write.csv(meta, file = file.path(outdir, "metadata.csv"), row.names = FALSE)

# ---- report -----------------------------------------------------------------
cat(sprintf("Wrote toy dataset to %s/\n", outdir))
cat(sprintf("  counts.tsv    %d peaks x %d samples\n", n_peaks, n_samples))
cat(sprintf("  peaks.bed     %d regions\n", n_peaks))
cat(sprintf("  metadata.csv  %d samples (%d control / %d treatment)\n",
            n_samples, half, half))
cat(sprintf("  true DARs:    %d gained, %d lost (should be recovered downstream)\n",
            length(gained_idx), length(lost_idx)))
cat(sprintf("  latent NMF programs: 2 (should appear in nmf/ outputs)\n"))
