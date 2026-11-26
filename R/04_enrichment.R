# Maps axolotl probesets to human Entrez IDs through the platform's presumptive
# ortholog symbols, collapses to one row per gene, and tests GO Biological Process
# over-representation with limma::goana against an array-specific universe.
source("R/config.R")
suppressPackageStartupMessages({
  library(limma)
  library(GO.db)
  library(org.Hs.eg.db)
})

de_tables <- readRDS(file.path(paths$processed, "de_tables.rds"))
annotation <- read.csv(file.path(paths$processed, "probe_annotation.csv"), stringsAsFactors = FALSE)

symbols <- unique(na.omit(annotation$symbol))
by_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, symbols, "ENTREZID", "SYMBOL", multiVals = "first")
by_alias <- AnnotationDbi::mapIds(org.Hs.eg.db, symbols[is.na(by_symbol)], "ENTREZID", "ALIAS", multiVals = "list")
by_alias <- vapply(by_alias, function(ids) if (length(ids) == 1 && !is.na(ids)) ids else NA_character_, "")
symbol_to_entrez <- c(by_symbol[!is.na(by_symbol)], by_alias)

annotation$entrez <- unname(symbol_to_entrez[annotation$symbol])
mapping_summary <- data.frame(
  step = c("Probesets on array", "with a presumptive human ortholog symbol",
           "symbol matched to a current HGNC symbol", "symbol matched only as an unambiguous alias",
           "probesets mapped to an Entrez gene", "distinct Entrez genes (enrichment universe)"),
  n = c(nrow(annotation), sum(!is.na(annotation$symbol)),
        sum(annotation$symbol %in% names(by_symbol)[!is.na(by_symbol)]),
        sum(annotation$symbol %in% names(by_alias)[!is.na(by_alias)]),
        sum(!is.na(annotation$entrez)), length(unique(na.omit(annotation$entrez)))))
write.csv(mapping_summary, file.path(paths$results, "ortholog_mapping_summary.csv"), row.names = FALSE)
print(mapping_summary)
universe <- unique(na.omit(annotation$entrez))

gene_level <- function(tt) {
  tt$entrez <- annotation$entrez[match(tt$probe_id, annotation$probe_id)]
  tt <- tt[!is.na(tt$entrez), ]
  tt <- tt[order(tt$P.Value), ]
  tt[!duplicated(tt$entrez), ]
}

run_goana <- function(id) {
  genes <- gene_level(de_tables[[id]])
  sets <- list(Up = genes$entrez[genes$status == "Up"], Down = genes$entrez[genes$status == "Down"])
  go <- goana(sets, universe = universe, species = "Hs")
  go <- go[go$Ont == "BP" & go$N >= 10 & go$N <= 500, ]
  go$FDR.Up <- p.adjust(go$P.Up, "BH")
  go$FDR.Down <- p.adjust(go$P.Down, "BH")
  out <- rbind(
    data.frame(go_id = rownames(go), term = go$Term, direction = "Up", genes_in_term = go$N,
               de_genes = go$Up, de_list_size = length(sets$Up), p_value = go$P.Up, fdr = go$FDR.Up),
    data.frame(go_id = rownames(go), term = go$Term, direction = "Down", genes_in_term = go$N,
               de_genes = go$Down, de_list_size = length(sets$Down), p_value = go$P.Down, fdr = go$FDR.Down))
  out$contrast <- id
  out[order(out$direction, out$p_value), ]
}

enrichment <- do.call(rbind, lapply(contrasts_registry$id, run_goana))
significant <- enrichment[enrichment$fdr < thresholds$adj_p, ]
write.csv(significant, file.path(paths$results, "go_bp_enrichment_fdr05.csv"), row.names = FALSE)
print(table(significant$contrast, significant$direction))

top_terms <- do.call(rbind, lapply(split(significant, list(significant$contrast, significant$direction), drop = TRUE),
                                   function(d) head(d[order(d$p_value), ], 6)))
plot_df <- enrichment[enrichment$go_id %in% top_terms$go_id, ]
plot_df <- merge(plot_df, contrasts_registry[, c("id", "label")], by.x = "contrast", by.y = "id")
plot_df$label <- factor(plot_df$label, levels = contrasts_registry$label)
plot_df$gene_ratio <- plot_df$de_genes / plot_df$de_list_size
plot_df$term <- factor(plot_df$term, levels = unique(top_terms$term[order(top_terms$direction, -top_terms$p_value)]))
plot_df$significant <- plot_df$fdr < thresholds$adj_p

save_png(
  ggplot(plot_df, aes(label, term)) +
    geom_point(aes(size = gene_ratio, colour = -log10(fdr), alpha = significant)) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
    scale_colour_gradient(low = "#9EC5E0", high = "#B2182B", name = "-log10 FDR") +
    scale_size_continuous(name = "DE genes / list size", range = c(1, 6)) +
    facet_wrap(~ direction, labeller = labeller(direction = c(Up = "Up-regulated genes", Down = "Down-regulated genes"))) +
    labs(title = "GO Biological Process over-representation (human orthologs)",
         subtitle = "Top 6 terms per contrast and direction; faded points have FDR >= 0.05",
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right"),
  "go_enrichment", width = 13, height = max(5, 0.28 * nlevels(plot_df$term) + 2.5))
