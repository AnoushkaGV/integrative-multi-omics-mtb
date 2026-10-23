# =============================================================================
# Runs the five analysis scripts in order and saves a log of the console output.
# Run from the repository root:
#   source("scripts/run_all.R")
# Tables and figures are written to output_revised/<analysis>/{tables,figures}.
# =============================================================================
scripts <- c("01_meta_mrna.R",     # Fig 2a-c, Supp Table S1 (must run first)
             "02_cibersort.R",     # Fig 2d, Supp Table S8
             "03_mirna.R",         # Fig 4a-c, Supp Tables S2-S3
             "04_trfs.R",          # Fig 5a-f, Supp Table S4
             "05_predictions.R")   # Fig 3, Table 2, Supp Tables S5-S7, Supp Figs S1-S2

dir.create("output_revised", showWarnings = FALSE)
log_file <- file.path("output_revised", paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt"))
con <- file(log_file, open = "wt")
sink(con, split = TRUE); sink(con, type = "message")

for (s in scripts) {
  cat("\n\n#################### ", s, " ####################\n", sep = "")
  t0 <- Sys.time()
  ok <- tryCatch({
    source(file.path("scripts", s), echo = FALSE, print.eval = TRUE); TRUE
  }, error = function(e) {
    # print errors to the console as well as the log
    cat("\n*** ERROR in ", s, ": ", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  if (!ok) { sink(type = "message"); sink(); close(con)
             stop("Stopped at ", s, " -- see the message above and ", log_file, call. = FALSE) }
  graphics.off()
  cat(sprintf("\n%s finished in %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(list = setdiff(ls(), c("scripts", "s", "log_file", "con")))
}

cat("\n\n#################### sessionInfo ####################\n")
print(sessionInfo())
sink(type = "message"); sink(); close(con)
cat("Done. Log saved to", log_file, "\n")
