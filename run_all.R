# Runs the whole pipeline from the repository root: Rscript run_all.R
for (step in c("R/01_download.R", "R/02_preprocess_qc.R",
               "R/03_differential_expression.R", "R/04_enrichment.R")) {
  message("\n==> ", step)
  source(step, local = new.env())
}
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
