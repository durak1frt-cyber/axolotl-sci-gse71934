# Reads CEL files with the GEO-supplied custom CDF, plots raw-intensity QC,
# runs RMA, then plots normalised QC (RLE, PCA, sample correlation clustering).
source("R/config.R")
suppressPackageStartupMessages({
  library(affy)
  library(makecdfenv)
  library(pheatmap)
})

samples <- read_samples()
cel_paths <- file.path(paths$raw, "cel", samples$cel_file)

cdf_env <- make.cdf.env("AMBY_002a520748F.cdf", cdf.path = paths$raw)
raw <- ReadAffy(filenames = cel_paths, sampleNames = samples$gsm, cdfname = geo$cdf_name)
# affy resolves a non-Bioconductor CDF by looking up its name in the global environment.
assign(cdfName(raw), cdf_env, envir = globalenv())

write.csv(data.frame(samples[, c("gsm", "group", "replicate")], scan_date = protocolData(raw)$ScanDate),
          file.path(paths$results, "array_scan_dates.csv"), row.names = FALSE)

long_intensities <- function(mat, n = 50000) {
  set.seed(1)
  rows <- sample(nrow(mat), min(n, nrow(mat)))
  data.frame(
    gsm   = rep(colnames(mat), each = length(rows)),
    value = as.vector(mat[rows, ])
  ) |> merge(samples[, c("gsm", "group", "label")], by = "gsm")
}

intensity_plots <- function(df, stage, value_label) {
  df$label <- factor(df$label, levels = samples$label)
  box <- ggplot(df, aes(label, value, fill = group)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
    scale_fill_manual(values = group_colours, labels = group_labels) +
    labs(title = paste0("Per-array intensity distribution (", stage, ")"),
         x = NULL, y = value_label, fill = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dens <- ggplot(df, aes(value, group = gsm, colour = group)) +
    geom_density(linewidth = 0.5) +
    scale_colour_manual(values = group_colours, labels = group_labels) +
    labs(title = paste0("Intensity density (", stage, ")"), x = value_label, y = "Density", colour = NULL)
  list(box = box, dens = dens)
}

raw_plots <- intensity_plots(long_intensities(log2(pm(raw))), "raw PM probes", "log2 intensity")
save_png(raw_plots$box, "qc_boxplot_raw", width = 9, height = 5)
save_png(raw_plots$dens, "qc_density_raw", width = 8, height = 5)

eset <- rma(raw)
expr <- exprs(eset)
stopifnot(identical(colnames(expr), samples$gsm))

norm_plots <- intensity_plots(long_intensities(expr), "after RMA", "log2 expression")
save_png(norm_plots$box, "qc_boxplot_rma", width = 9, height = 5)
save_png(norm_plots$dens, "qc_density_rma", width = 8, height = 5)

rle <- long_intensities(expr - matrixStats::rowMedians(expr))
rle$label <- factor(rle$label, levels = samples$label)
save_png(
  ggplot(rle, aes(label, value, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = group_colours, labels = group_labels) +
    labs(title = "Relative log expression (RLE)", subtitle = "Each probeset minus its median across arrays",
         x = NULL, y = "RLE", fill = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "qc_rle", width = 9, height = 5)

pca <- prcomp(t(expr), center = TRUE, scale. = FALSE)
var_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- cbind(samples, pca$x[, 1:3])
save_png(
  ggplot(pca_df, aes(PC1, PC2, colour = group)) +
    geom_point(size = 3.5) +
    ggrepel::geom_text_repel(aes(label = replicate), size = 3, show.legend = FALSE) +
    scale_colour_manual(values = group_colours, labels = group_labels) +
    labs(title = "PCA of RMA-normalised expression",
         subtitle = sprintf("All %d probesets; labels are replicate numbers", nrow(expr)),
         x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]), colour = NULL),
  "qc_pca", width = 7, height = 6)
write.csv(data.frame(pc = paste0("PC", seq_along(var_explained)), percent_variance = var_explained),
          file.path(paths$results, "pca_variance_explained.csv"), row.names = FALSE)

cor_mat <- cor(expr, method = "pearson")
dimnames(cor_mat) <- list(samples$label, samples$label)
annotation_col <- data.frame(Group = group_labels[as.character(samples$group)], row.names = samples$label)
pheatmap(cor_mat,
         clustering_distance_rows = as.dist(1 - cor_mat),
         clustering_distance_cols = as.dist(1 - cor_mat),
         clustering_method = "average",
         annotation_col = annotation_col, annotation_row = annotation_col,
         annotation_colors = list(Group = setNames(group_colours, group_labels[names(group_colours)])),
         main = "Sample-to-sample Pearson correlation (average-linkage clustering)",
         filename = file.path(paths$figures, "qc_sample_correlation.png"),
         width = 9, height = 7.5)

saveRDS(eset, file.path(paths$processed, "eset_rma.rds"))
message("RMA complete: ", nrow(expr), " probesets x ", ncol(expr), " arrays")
