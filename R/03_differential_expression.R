# limma differential expression for every contrast in contrasts_registry,
# with volcano plots, DE counts, a heatmap of top time-course genes and a FOS check.
source("R/config.R")
suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(pheatmap)
})

samples <- read_samples()
eset <- readRDS(file.path(paths$processed, "eset_rma.rds"))
expr <- exprs(eset)
annotation <- read.csv(file.path(paths$processed, "probe_annotation.csv"), stringsAsFactors = FALSE)

annot <- annotation[match(rownames(expr), annotation$probe_id), ]
stopifnot(!anyNA(annot$probe_id))
display_name <- ifelse(is.na(annot$symbol), rownames(expr), annot$symbol)

design <- model.matrix(~ 0 + group, data = samples)
colnames(design) <- levels(samples$group)
# ivermectin_vs_uninjured feeds only the ivermectin summary at the end, not the per-contrast outputs.
contrast_matrix <- makeContrasts(contrasts = c(contrasts_registry$expression, "ivermectin_dpi1 - uninjured"),
                                 levels = design)
colnames(contrast_matrix) <- c(contrasts_registry$id, "ivermectin_vs_uninjured")

fit <- eBayes(contrasts.fit(lmFit(expr, design), contrast_matrix))

classify <- function(tt) {
  ifelse(tt$adj.P.Val < thresholds$adj_p & tt$logFC >= thresholds$log2fc, "Up",
         ifelse(tt$adj.P.Val < thresholds$adj_p & tt$logFC <= -thresholds$log2fc, "Down", "Not significant"))
}

de_tables <- lapply(setNames(contrasts_registry$id, contrasts_registry$id), function(id) {
  tt <- topTable(fit, coef = id, number = Inf, sort.by = "none")
  tt <- data.frame(probe_id = rownames(tt), symbol = annot$symbol, gene_name = annot$gene_name,
                   tt[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], row.names = NULL)
  tt$status <- classify(tt)
  tt <- tt[order(tt$P.Value), ]
  write.csv(tt, file.path(paths$processed, paste0("de_", id, ".csv")), row.names = FALSE)
  write.csv(head(tt[tt$status != "Not significant", ], 50),
            file.path(paths$results, paste0("top50_de_", id, ".csv")), row.names = FALSE)
  tt
})
saveRDS(de_tables, file.path(paths$processed, "de_tables.rds"))

de_counts <- do.call(rbind, lapply(contrasts_registry$id, function(id) {
  tt <- de_tables[[id]]
  data.frame(contrast = id,
             label = contrasts_registry$label[contrasts_registry$id == id],
             up = sum(tt$status == "Up"), down = sum(tt$status == "Down"),
             total = sum(tt$status != "Not significant"),
             unique_symbols = length(unique(na.omit(tt$symbol[tt$status != "Not significant"]))))
}))
write.csv(de_counts, file.path(paths$results, "de_counts.csv"), row.names = FALSE)
print(de_counts)

counts_long <- rbind(
  data.frame(label = de_counts$label, direction = "Up", n = de_counts$up),
  data.frame(label = de_counts$label, direction = "Down", n = de_counts$down))
counts_long$label <- factor(counts_long$label, levels = contrasts_registry$label)
save_png(
  ggplot(counts_long, aes(label, n, fill = direction)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = n), position = position_dodge(width = 0.75), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c(Up = "#C8463D", Down = "#2F6EA6")) +
    labs(title = "Differentially expressed probesets per contrast",
         subtitle = sprintf("limma, BH-adjusted p < %.2f and |log2FC| >= %g", thresholds$adj_p, thresholds$log2fc),
         x = NULL, y = "Probesets", fill = NULL),
  "de_counts", width = 8, height = 5)

volcano_df <- do.call(rbind, lapply(contrasts_registry$id, function(id) {
  tt <- de_tables[[id]]
  tt$label <- contrasts_registry$label[contrasts_registry$id == id]
  tt$display <- ifelse(is.na(tt$symbol), tt$probe_id, tt$symbol)
  tt$show_label <- cumsum(tt$status != "Not significant") <= 12 & tt$status != "Not significant"
  tt
}))
volcano_df$label <- factor(volcano_df$label, levels = contrasts_registry$label)
volcano_df$status <- factor(volcano_df$status, levels = c("Up", "Down", "Not significant"))
save_png(
  ggplot(volcano_df, aes(logFC, -log10(adj.P.Val), colour = status)) +
    geom_point(data = subset(volcano_df, status == "Not significant"), size = 0.6, alpha = 0.4) +
    geom_point(data = subset(volcano_df, status != "Not significant"), size = 0.9, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1) * thresholds$log2fc, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(thresholds$adj_p), linetype = "dashed", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = subset(volcano_df, show_label), aes(label = display),
                             size = 2.8, colour = "black", max.overlaps = 30, min.segment.length = 0) +
    scale_colour_manual(values = c(Up = "#C8463D", Down = "#2F6EA6", "Not significant" = "grey70")) +
    facet_wrap(~ label, ncol = 2) +
    labs(title = "Volcano plots",
         subtitle = "Labels: 12 most significant DE probesets per contrast (presumptive human ortholog symbol)",
         x = "log2 fold change", y = "-log10 BH-adjusted p-value", colour = NULL),
  "volcano_plots", width = 11, height = 9.5)

time_course_ids <- contrasts_registry$id[contrasts_registry$time_course]
f_table <- topTable(fit, coef = time_course_ids, number = Inf, sort.by = "F")
top_probes <- head(rownames(f_table)[f_table$adj.P.Val < thresholds$adj_p], 50)
z <- t(scale(t(expr[top_probes, ])))
rownames(z) <- make.unique(display_name[match(top_probes, rownames(expr))])
colnames(z) <- samples$label
annotation_col <- data.frame(Group = group_labels[as.character(samples$group)], row.names = samples$label)
pheatmap(z,
         cluster_cols = FALSE,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         annotation_col = annotation_col,
         annotation_colors = list(Group = setNames(group_colours, group_labels[names(group_colours)])),
         color = colorRampPalette(c("#2F6EA6", "white", "#C8463D"))(101),
         breaks = seq(-2.5, 2.5, length.out = 102),
         gaps_col = cumsum(table(samples$group))[-nlevels(samples$group)],
         fontsize_row = 7,
         main = "Top 50 time-course probesets (moderated F-test across 1/3/7 dpi vs uninjured), row z-scores",
         filename = file.path(paths$figures, "heatmap_top_de.png"),
         width = 10, height = 10)

fos_probes <- which(annot$symbol == "FOS")
if (length(fos_probes)) {
  fos <- data.frame(
    probe_id = rep(rownames(expr)[fos_probes], each = ncol(expr)),
    value    = as.vector(t(expr[fos_probes, , drop = FALSE])),
    group    = rep(samples$group, length(fos_probes)))
  save_png(
    ggplot(fos, aes(group, value, colour = group)) +
      geom_jitter(width = 0.12, size = 2.5) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "black", linewidth = 0.3) +
      scale_x_discrete(labels = group_labels) +
      scale_colour_manual(values = group_colours, guide = "none") +
      facet_wrap(~ probe_id) +
      labs(title = "FOS (c-Fos) probesets",
           subtitle = "The original study reports c-Fos as an early injury-response gene",
           x = NULL, y = "RMA log2 expression") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)),
    "fos_expression", width = 4 + 2.5 * length(fos_probes), height = 5)
  fos_stats <- do.call(rbind, lapply(contrasts_registry$id, function(id) {
    tt <- de_tables[[id]]
    cbind(contrast = id, tt[tt$probe_id %in% rownames(expr)[fos_probes],
                            c("probe_id", "logFC", "adj.P.Val", "status")])
  }))
  write.csv(fos_stats, file.path(paths$results, "fos_de_stats.csv"), row.names = FALSE)
}

injury <- de_tables$dpi1_vs_uninjured
injury_de <- injury$probe_id[injury$status != "Not significant"]
ivermectin_vs_uninjured <- topTable(fit, coef = "ivermectin_vs_uninjured", number = Inf)
ivermectin_summary <- data.frame(
  metric = c("DE probesets, 1 dpi vs uninjured",
             "DE probesets, ivermectin 1 dpi vs uninjured (same thresholds)",
             "Median |log2FC| of 1 dpi DE probesets, 1 dpi vs uninjured",
             "Median |log2FC| of the same probesets, ivermectin 1 dpi vs uninjured"),
  value = c(length(injury_de),
            sum(ivermectin_vs_uninjured$adj.P.Val < thresholds$adj_p &
                  abs(ivermectin_vs_uninjured$logFC) >= thresholds$log2fc),
            median(abs(injury$logFC[injury$probe_id %in% injury_de])),
            median(abs(ivermectin_vs_uninjured[injury_de, "logFC"]))))
write.csv(ivermectin_summary, file.path(paths$results, "ivermectin_vs_injury_response.csv"), row.names = FALSE)
print(ivermectin_summary)
