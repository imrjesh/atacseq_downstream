#!/usr/bin/env Rscript
# =============================================================================
# atacseq_downstream_pipeline.R
#
# Bulk ATAC-seq DOWNSTREAM analysis pipeline: from a peak x sample count matrix
# to differentially accessible regions (DARs), accessibility programs (NMF),
# peak annotation, functional enrichment, TF motif activity, and an integrated
# prioritization of candidate regulatory regions / transcription factors.
#
# SCOPE NOTE:
#   This pipeline starts FROM a count matrix. Peak calling (MACS2), alignment,
#   and matrix construction are UPSTREAM steps documented in `workflow/upstream/`.
#   A count matrix already implies peaks exist, so peak calling is not repeated
#   here.
#
# TERMINOLOGY NOTE:
#   ATAC-seq measures chromatin accessibility, not gene expression. The
#   differential test yields DIFFERENTIALLY ACCESSIBLE REGIONS (DARs), not DGE.
#
# Usage:
#   Rscript atacseq_downstream_pipeline.R config.yaml
#
# Author: Rajesh Kumar|Ph.D.   |   License: MIT
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
  library(SummarizedExperiment)
  library(matrixStats)
  library(ComplexHeatmap)
  library(circlize)
  library(NMF)
  library(ChIPseeker)
  library(clusterProfiler)
  library(chromVAR)
  library(motifmatchr)
  library(TFBSTools)
  library(GenomicRanges)
})

# -----------------------------------------------------------------------------
# 0. CONFIG
# -----------------------------------------------------------------------------
# All run-specific parameters live in an external config.yaml so the code stays
# generic. Example config.yaml is documented in the repo README. Expected keys:
#
#   counts:        path to peak x sample count matrix (rows = peaks, cols = samples)
#   peaks:         path to peak BED (chr,start,end matching count matrix rownames)
#   metadata:      path to sample metadata CSV (must contain 'sample' + 'condition')
#   outdir:        output directory
#   condition:     name of the metadata column to contrast on (e.g. "condition")
#   reference:     reference level of that column (e.g. "control")
#   genome:        BSgenome package name (e.g. "BSgenome.Hsapiens.UCSC.hg38")
#   txdb:          TxDb package name (e.g. "TxDb.Hsapiens.UCSC.hg38.knownGene")
#   orgdb:         OrgDb package name (e.g. "org.Hs.eg.db")
#   species:       JASPAR species tax id (e.g. 9606 for human)
#   min_count:     min total count per peak to keep (default 10)
#   min_samples:   peak must exceed min_count in >= this many samples (default 3)
#   padj:          adjusted p-value cutoff for DARs (default 0.05)
#   lfc:           |log2FC| cutoff for DARs (default 1)
#   n_top_var:     number of most-variable peaks used for NMF/heatmaps (default 3000)
#   nmf_ranks:     vector of ranks to survey for NMF (e.g. [2,3,4,5,6])
#   seed:          random seed for reproducibility (default 42)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript atacseq_downstream_pipeline.R config.yaml")
cfg <- yaml::read_yaml(args[1])

# Fill defaults for optional keys
defaults <- list(min_count = 10, min_samples = 3, padj = 0.05, lfc = 1,
                 n_top_var = 3000, nmf_ranks = 2:6, seed = 42)
for (k in names(defaults)) if (is.null(cfg[[k]])) cfg[[k]] <- defaults[[k]]

set.seed(cfg$seed)
dir.create(cfg$outdir, showWarnings = FALSE, recursive = TRUE)
subdirs <- c("qc", "dars", "nmf", "annotation", "enrichment", "motifs", "prioritization")
for (s in subdirs) dir.create(file.path(cfg$outdir, s), showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# -----------------------------------------------------------------------------
# 1. LOAD INPUTS
# -----------------------------------------------------------------------------
load_inputs <- function(cfg) {
  log_msg("Loading count matrix, peaks, and metadata")
  counts <- as.matrix(read.delim(cfg$counts, row.names = 1, check.names = FALSE))
  meta   <- read.csv(cfg$metadata, stringsAsFactors = FALSE)
  rownames(meta) <- meta$sample

  # Align matrix columns to metadata rows
  common <- intersect(colnames(counts), rownames(meta))
  if (length(common) < 2) stop("Fewer than 2 samples shared between counts and metadata.")
  counts <- counts[, common, drop = FALSE]
  meta   <- meta[common, , drop = FALSE]

  # Parse peak coordinates into a GRanges (rownames expected as chr:start-end or a BED)
  if (!is.null(cfg$peaks) && file.exists(cfg$peaks)) {
    bed <- read.delim(cfg$peaks, header = FALSE)
    peaks <- GRanges(bed[[1]], IRanges(bed[[2]] + 1, bed[[3]]))  # BED is 0-based
    names(peaks) <- rownames(counts)
  } else {
    peaks <- GRanges(rownames(counts))  # assumes "chr:start-end" rownames
  }
  list(counts = counts, meta = meta, peaks = peaks)
}

# -----------------------------------------------------------------------------
# 2. QC & FILTERING
# -----------------------------------------------------------------------------
qc_and_filter <- function(counts, meta, cfg) {
  log_msg("QC and peak filtering")

  # Library size QC
  libsize <- colSums(counts)
  pdf(file.path(cfg$outdir, "qc", "library_sizes.pdf"), width = 7, height = 5)
  barplot(libsize / 1e6, las = 2, ylab = "Library size (millions)", main = "Fragments in peaks")
  dev.off()

  # Keep peaks with enough signal in enough samples (removes noise peaks)
  keep <- rowSums(counts >= cfg$min_count) >= cfg$min_samples
  log_msg(sprintf("Retained %d / %d peaks after filtering", sum(keep), nrow(counts)))
  counts[keep, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 3. NORMALIZATION, VST, AND EXPLORATORY PCA
# -----------------------------------------------------------------------------
build_dds_and_vst <- function(counts, meta, cfg) {
  log_msg("Building DESeq2 dataset and variance-stabilizing transform")

  meta[[cfg$condition]] <- relevel(factor(meta[[cfg$condition]]), ref = cfg$reference)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta,
                                design = as.formula(paste0("~", cfg$condition)))
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)

  # PCA on VST
  pca <- plotPCA(vsd, intgroup = cfg$condition, returnData = TRUE)
  pct <- round(100 * attr(pca, "percentVar"))
  pdf(file.path(cfg$outdir, "qc", "pca.pdf"), width = 6, height = 5)
  with(pca, {
    plot(PC1, PC2, col = as.integer(factor(group)), pch = 19,
         xlab = sprintf("PC1 (%d%%)", pct[1]), ylab = sprintf("PC2 (%d%%)", pct[2]),
         main = "PCA (VST-normalized accessibility)")
    text(PC1, PC2, labels = name, pos = 3, cex = 0.6)
  })
  dev.off()

  # Sample-sample correlation heatmap
  cor_mat <- cor(assay(vsd))
  pdf(file.path(cfg$outdir, "qc", "sample_correlation.pdf"), width = 6, height = 5)
  draw(Heatmap(cor_mat, name = "Pearson", column_title = "Sample correlation (VST)"))
  dev.off()

  list(dds = dds, vsd = vsd)
}

# -----------------------------------------------------------------------------
# 4. DIFFERENTIAL ACCESSIBILITY (DARs)
# -----------------------------------------------------------------------------
call_dars <- function(dds, cfg) {
  log_msg("Testing for differentially accessible regions (DARs)")

  res <- results(dds, alpha = cfg$padj)
  res <- lfcShrink(dds, coef = resultsNames(dds)[length(resultsNames(dds))],
                   res = res, type = "apeglm")
  res_df <- as.data.frame(res)
  res_df$peak <- rownames(res_df)
  res_df$significant <- with(res_df,
                             !is.na(padj) & padj < cfg$padj & abs(log2FoldChange) >= cfg$lfc)
  res_df$direction <- ifelse(res_df$log2FoldChange > 0, "gained", "lost")

  write.csv(res_df, file.path(cfg$outdir, "dars", "differential_accessibility.csv"),
            row.names = FALSE)

  # Volcano plot
  pdf(file.path(cfg$outdir, "dars", "volcano.pdf"), width = 6, height = 5)
  with(res_df, {
    plot(log2FoldChange, -log10(padj), pch = 20, cex = 0.4,
         col = ifelse(significant, ifelse(direction == "gained", "firebrick", "steelblue"), "grey70"),
         xlab = "log2 fold-change (accessibility)", ylab = "-log10 adjusted p",
         main = "Differentially accessible regions")
    abline(h = -log10(cfg$padj), lty = 2); abline(v = c(-cfg$lfc, cfg$lfc), lty = 2)
  })
  dev.off()

  log_msg(sprintf("Found %d significant DARs (%d gained, %d lost)",
                  sum(res_df$significant, na.rm = TRUE),
                  sum(res_df$significant & res_df$direction == "gained", na.rm = TRUE),
                  sum(res_df$significant & res_df$direction == "lost", na.rm = TRUE)))
  res_df
}

# -----------------------------------------------------------------------------
# 5. HEATMAP OF TOP DARs
# -----------------------------------------------------------------------------
dar_heatmap <- function(vsd, res_df, meta, cfg) {
  log_msg("Drawing DAR heatmap")
  sig <- res_df$peak[which(res_df$significant)]
  if (length(sig) < 2) { log_msg("Too few DARs for a heatmap; skipping"); return(invisible()) }

  mat <- assay(vsd)[sig, , drop = FALSE]
  mat <- t(scale(t(mat)))  # row z-score
  ann <- HeatmapAnnotation(condition = meta[[cfg$condition]])
  col_fun <- colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick"))

  pdf(file.path(cfg$outdir, "dars", "dar_heatmap.pdf"), width = 7, height = 8)
  draw(Heatmap(mat, name = "z-score", top_annotation = ann, col = col_fun,
               show_row_names = FALSE, column_title = "Top DARs (row-scaled VST)"))
  dev.off()
}

# -----------------------------------------------------------------------------
# 6. NMF: UNSUPERVISED ACCESSIBILITY PROGRAMS
# -----------------------------------------------------------------------------
# NMF factorizes the non-negative (normalized) count matrix into "programs"
# (metafeatures) x sample-loadings, revealing coordinated accessibility modules
# beyond a single pairwise contrast.
run_nmf <- function(dds, cfg) {
  log_msg("Running NMF on top-variable peaks")

  norm <- counts(dds, normalized = TRUE)
  vars <- rowVars(norm)
  top  <- head(order(vars, decreasing = TRUE), cfg$n_top_var)
  mat  <- norm[top, , drop = FALSE]
  mat  <- mat[rowSums(mat) > 0, , drop = FALSE]  # NMF requires non-negative, non-empty rows

  # Rank survey (cophenetic / dispersion) to choose factor number
  ranks <- as.integer(cfg$nmf_ranks)
  estim <- nmf(mat, rank = ranks, nrun = 20, seed = cfg$seed, .options = "p")
  pdf(file.path(cfg$outdir, "nmf", "rank_survey.pdf"), width = 7, height = 5)
  plot(estim)
  dev.off()

  # Fit at the best rank by cophenetic coefficient
  best_rank <- ranks[which.max(sapply(estim$fit, function(f) cophcor(f)))]
  log_msg(sprintf("Selected NMF rank = %d", best_rank))
  fit <- nmf(mat, rank = best_rank, nrun = 50, seed = cfg$seed)

  # Sample-to-program assignment from the coefficient (H) matrix
  H <- coef(fit)
  assign_df <- data.frame(sample = colnames(H),
                          program = paste0("P", apply(H, 2, which.max)))
  write.csv(assign_df, file.path(cfg$outdir, "nmf", "sample_program_assignment.csv"),
            row.names = FALSE)

  # Program x peak feature scores (W matrix) — top peaks per program
  W <- basis(fit)
  colnames(W) <- paste0("P", seq_len(ncol(W)))
  feat <- lapply(seq_len(ncol(W)), function(i) {
    ord <- order(W[, i], decreasing = TRUE)
    data.frame(program = colnames(W)[i], peak = rownames(W)[head(ord, 100)],
               score = W[head(ord, 100), i])
  })
  feat <- do.call(rbind, feat)
  write.csv(feat, file.path(cfg$outdir, "nmf", "program_top_peaks.csv"), row.names = FALSE)

  pdf(file.path(cfg$outdir, "nmf", "consensus_heatmap.pdf"), width = 7, height = 6)
  consensusmap(fit)
  dev.off()

  list(fit = fit, rank = best_rank, features = feat, assignment = assign_df)
}

# -----------------------------------------------------------------------------
# 7. PEAK ANNOTATION
# -----------------------------------------------------------------------------
annotate_peaks <- function(peaks, res_df, cfg) {
  log_msg("Annotating peaks to genomic features and nearest genes")
  if (is.null(cfg$txdb)) { log_msg("No TxDb specified; skipping annotation"); return(NULL) }

  suppressPackageStartupMessages(library(cfg$txdb, character.only = TRUE))
  txdb <- get(cfg$txdb)

  sig_peaks <- peaks[names(peaks) %in% res_df$peak[which(res_df$significant)]]
  if (length(sig_peaks) == 0) { log_msg("No significant peaks to annotate"); return(NULL) }

  anno <- annotatePeak(sig_peaks, TxDb = txdb, tssRegion = c(-3000, 3000),
                       annoDb = cfg$orgdb, verbose = FALSE)
  pdf(file.path(cfg$outdir, "annotation", "feature_distribution.pdf"), width = 7, height = 4)
  plotAnnoBar(anno)
  dev.off()

  anno_df <- as.data.frame(anno)
  write.csv(anno_df, file.path(cfg$outdir, "annotation", "dar_annotation.csv"), row.names = FALSE)
  anno_df
}

# -----------------------------------------------------------------------------
# 8. FUNCTIONAL ENRICHMENT OF DAR-LINKED GENES
# -----------------------------------------------------------------------------
functional_enrichment <- function(anno_df, cfg) {
  log_msg("GO enrichment of DAR-linked genes")
  if (is.null(anno_df) || is.null(cfg$orgdb)) return(invisible())
  suppressPackageStartupMessages(library(cfg$orgdb, character.only = TRUE))

  genes <- unique(na.omit(anno_df$geneId))
  if (length(genes) < 10) { log_msg("Too few genes for enrichment; skipping"); return(invisible()) }

  ego <- enrichGO(gene = genes, OrgDb = cfg$orgdb, keyType = "ENTREZID",
                  ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file.path(cfg$outdir, "enrichment", "go_bp.csv"), row.names = FALSE)
    pdf(file.path(cfg$outdir, "enrichment", "go_dotplot.pdf"), width = 8, height = 6)
    print(dotplot(ego, showCategory = 20))
    dev.off()
  }
}

# -----------------------------------------------------------------------------
# 9. TF MOTIF ACTIVITY (chromVAR + JASPAR)
# -----------------------------------------------------------------------------
# For bulk ATAC we quantify per-sample deviations in accessibility across peaks
# containing each TF motif — a robust readout of differential TF activity.
motif_activity <- function(counts, peaks, meta, cfg) {
  log_msg("Computing TF motif accessibility deviations (chromVAR)")
  if (is.null(cfg$genome)) { log_msg("No BSgenome specified; skipping motif analysis"); return(NULL) }
  suppressPackageStartupMessages(library(cfg$genome, character.only = TRUE))
  genome <- get(cfg$genome)

  peaks <- peaks[names(peaks) %in% rownames(counts)]
  peaks <- resize(peaks, width = 500, fix = "center")  # standardize widths for GC bias
  cmat  <- counts[names(peaks), , drop = FALSE]

  se <- SummarizedExperiment(assays = list(counts = cmat),
                             rowRanges = peaks, colData = meta)
  se <- addGCBias(se, genome = genome)
  se <- filterPeaks(se, non_overlapping = TRUE)

  motifs   <- getJasparMotifs(species = cfg$species)
  mot_ix   <- matchMotifs(motifs, se, genome = genome)
  dev      <- computeDeviations(object = se, annotations = mot_ix)
  var_tf   <- computeVariability(dev)

  write.csv(var_tf, file.path(cfg$outdir, "motifs", "tf_variability.csv"), row.names = FALSE)

  # Differential motif activity between conditions (t-test on deviation z-scores)
  z <- deviationScores(dev)
  grp <- meta[[cfg$condition]]
  levs <- levels(factor(grp))
  if (length(levs) == 2) {
    tt <- apply(z, 1, function(v) {
      tryCatch(t.test(v[grp == levs[2]], v[grp == levs[1]])$p.value, error = function(e) NA)
    })
    eff <- rowMeans(z[, grp == levs[2], drop = FALSE], na.rm = TRUE) -
           rowMeans(z[, grp == levs[1], drop = FALSE], na.rm = TRUE)
    tf_df <- data.frame(motif = rownames(z),
                        tf = TFBSTools::name(motifs)[match(rownames(z), names(motifs))],
                        delta_activity = eff, p = tt,
                        padj = p.adjust(tt, "BH"))
    tf_df <- tf_df[order(tf_df$padj), ]
    write.csv(tf_df, file.path(cfg$outdir, "motifs", "differential_tf_activity.csv"),
              row.names = FALSE)
    return(tf_df)
  }
  NULL
}

# -----------------------------------------------------------------------------
# 10. INTEGRATED PRIORITIZATION
# -----------------------------------------------------------------------------
# Combine evidence into a ranked shortlist of candidate regulatory regions:
#   - DAR strength (signed -log10 padj * log2FC)
#   - proximity to a gene (annotation)
#   - presence in an NMF program
# Produces a transparent, weighted score reviewers can trace.
prioritize <- function(res_df, anno_df, nmf_res, cfg) {
  log_msg("Prioritizing candidate regulatory regions")

  pr <- res_df[which(res_df$significant), c("peak", "log2FoldChange", "padj", "direction")]
  if (nrow(pr) == 0) { log_msg("No DARs to prioritize"); return(invisible()) }

  pr$dar_score <- abs(pr$log2FoldChange) * -log10(pr$padj + 1e-300)

  if (!is.null(anno_df)) {
    key <- paste0(anno_df$seqnames, ":", anno_df$start, "-", anno_df$end)
    pr$linked_gene <- anno_df$SYMBOL[match(pr$peak, key)]
    pr$dist_to_tss <- anno_df$distanceToTSS[match(pr$peak, key)]
  }
  if (!is.null(nmf_res)) {
    pr$nmf_program <- nmf_res$features$program[match(pr$peak, nmf_res$features$peak)]
  }

  # Normalize DAR score to 0-1 and add a small bonus for TSS-proximal, program-tagged peaks
  pr$score <- pr$dar_score / max(pr$dar_score, na.rm = TRUE)
  if (!is.null(pr$dist_to_tss)) pr$score <- pr$score + 0.1 * (abs(pr$dist_to_tss) < 50000)
  if (!is.null(pr$nmf_program)) pr$score <- pr$score + 0.1 * (!is.na(pr$nmf_program))
  pr <- pr[order(pr$score, decreasing = TRUE), ]

  write.csv(pr, file.path(cfg$outdir, "prioritization", "prioritized_regions.csv"),
            row.names = FALSE)
  log_msg(sprintf("Wrote %d prioritized regions", nrow(pr)))
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
main <- function(cfg) {
  inp <- load_inputs(cfg)
  counts_f <- qc_and_filter(inp$counts, inp$meta, cfg)
  dv <- build_dds_and_vst(counts_f, inp$meta, cfg)
  res_df <- call_dars(dv$dds, cfg)
  dar_heatmap(dv$vsd, res_df, inp$meta, cfg)
  nmf_res <- tryCatch(run_nmf(dv$dds, cfg), error = function(e) { log_msg("NMF failed:", conditionMessage(e)); NULL })
  anno_df <- annotate_peaks(inp$peaks, res_df, cfg)
  functional_enrichment(anno_df, cfg)
  tf_df   <- tryCatch(motif_activity(counts_f, inp$peaks, inp$meta, cfg),
                      error = function(e) { log_msg("Motif step failed:", conditionMessage(e)); NULL })
  prioritize(res_df, anno_df, nmf_res, cfg)

  # Reproducibility: always record the environment
  writeLines(capture.output(sessionInfo()),
             file.path(cfg$outdir, "sessionInfo.txt"))
  log_msg("Pipeline complete. Outputs in:", cfg$outdir)
}

if (sys.nframe() == 0) main(cfg)
