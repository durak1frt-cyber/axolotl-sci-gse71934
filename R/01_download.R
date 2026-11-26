# Fetches raw CEL files, the custom AMBY_002 CDF and the GPL15153 annotation from GEO.
# Re-running skips anything already on disk.
source("R/config.R")
options(timeout = 600)

fetch <- function(url, dest) {
  if (!file.exists(dest)) {
    message("Downloading ", basename(dest))
    download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

cel_dir <- file.path(paths$raw, "cel")
tar_path <- fetch(geo$raw_tar, file.path(paths$raw, "GSE71934_RAW.tar"))
if (!dir.exists(cel_dir)) untar(tar_path, exdir = cel_dir)

cdf_gz <- fetch(geo$cdf_gz, file.path(paths$raw, "AMBY_002a520748F.cdf.gz"))
# makecdfenv recognises the binary (XDA) CDF format only in an uncompressed file.
cdf_path <- sub("\\.gz$", "", cdf_gz)
if (!file.exists(cdf_path)) R.utils::gunzip(cdf_gz, destname = cdf_path, remove = FALSE)

gpl <- GEOquery::getGEO(geo$platform, destdir = paths$raw)
annotation <- GEOquery::Table(gpl)[, c("ID", "Contig_ID", "GeneSymbol", "GeneName")]
names(annotation) <- c("probe_id", "contig_id", "symbol", "gene_name")
annotation$symbol[annotation$symbol == ""] <- NA
write.csv(annotation, file.path(paths$processed, "probe_annotation.csv"), row.names = FALSE)

samples <- read_samples()
missing <- setdiff(samples$cel_file, list.files(cel_dir))
if (length(missing)) stop("CEL files listed in metadata/samples.csv but not in the GEO archive: ",
                          paste(missing, collapse = ", "))
message("Download complete: ", nrow(samples), " CEL files, ", nrow(annotation), " annotated probesets")
