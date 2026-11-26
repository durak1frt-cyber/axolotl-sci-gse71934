suppressPackageStartupMessages({
  library(ggplot2)
})

paths <- list(
  raw       = "data/raw",
  processed = "data/processed",
  results   = "results",
  figures   = "figures"
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

geo <- list(
  series   = "GSE71934",
  platform = "GPL15153",
  cdf_name = "AMBY_002a520748F",
  raw_tar  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71934/suppl/GSE71934_RAW.tar",
  cdf_gz   = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15153/suppl/GPL15153_AMBY_002a520748F.cdf.gz"
)

group_levels <- c("uninjured", "dpi1", "dpi3", "dpi7", "ivermectin_dpi1")
group_labels <- c(
  uninjured       = "Uninjured",
  dpi1            = "1 dpi",
  dpi3            = "3 dpi",
  dpi7            = "7 dpi",
  ivermectin_dpi1 = "1 dpi + ivermectin"
)
group_colours <- c(
  uninjured       = "#4D4D4D",
  dpi1            = "#F2A541",
  dpi3            = "#E0603A",
  dpi7            = "#A8324A",
  ivermectin_dpi1 = "#3B7EA1"
)

contrasts_registry <- data.frame(
  id         = c("dpi1_vs_uninjured", "dpi3_vs_uninjured", "dpi7_vs_uninjured", "ivermectin_vs_dpi1"),
  expression = c("dpi1 - uninjured", "dpi3 - uninjured", "dpi7 - uninjured", "ivermectin_dpi1 - dpi1"),
  label      = c("1 dpi vs uninjured", "3 dpi vs uninjured", "7 dpi vs uninjured", "Ivermectin 1 dpi vs 1 dpi"),
  time_course = c(TRUE, TRUE, TRUE, FALSE)
)

thresholds <- list(adj_p = 0.05, log2fc = 1)

read_samples <- function() {
  s <- read.csv("metadata/samples.csv", stringsAsFactors = FALSE)
  s$group <- factor(s$group, levels = group_levels)
  s$label <- paste0(group_labels[as.character(s$group)], " #", s$replicate)
  s
}

save_png <- function(plot, name, width = 8, height = 6) {
  ggsave(file.path(paths$figures, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 150, bg = "white")
}

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold"),
                  legend.position = "bottom"))
