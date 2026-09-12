options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "02e-build-tcga-source-modules-v2-no-cli"

# Archival clean-replay version.
# This file fixes only the final R string-concatenation bug in 02e v1.
# It is NOT needed to continue the current analysis because the original
# 02e v1 first run already completed module detection and wrote the scientific
# artifacts before its serialization-only error.
#
# To preserve first-run provenance, use:
#   02e1_finalize_tcga_source_modules_v1.R
#
# This v2 file intentionally exits unless ALLOW_FULL_REPLAY is manually set TRUE.

ALLOW_FULL_REPLAY <- FALSE

if (!ALLOW_FULL_REPLAY) {
  stop(
    "02e v1 already completed the first frozen module-detection computation. ",
    "Do not recompute it for the current analysis. Run ",
    "`02e1_finalize_tcga_source_modules_v1.R` instead. ",
    "Set ALLOW_FULL_REPLAY <- TRUE only for a future independent regression replay."
  )
}

# If an independent replay is required later, copy the scientific WGCNA call
# from 02e v1 unchanged and replace the faulty freeze_rule expression with:
#
# freeze_rule = paste0(
#   "Every non-grey module returned by the pre-frozen WGCNA contract is retained. ",
#   "No module is selected or discarded based on size beyond the pre-frozen detection rule, ",
#   "biology, target coverage, or target preservation."
# )
#
# The current project intentionally does not execute a second module-detection
# run before the first-run artifacts are finalized and hashed.
