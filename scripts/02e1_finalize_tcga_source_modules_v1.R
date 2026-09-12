options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "02e1-finalize-tcga-source-modules-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"
MODULE_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_tcga_source_modules_v1")
CONTRACT_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_wgcna_module_contract_v2")

NETWORK_RDS <- file.path(MODULE_ROOT, "tcga_source_wgcna_network_v1.rds")
MEMBERSHIP_TSV <- file.path(MODULE_ROOT, "tcga_source_module_membership_frozen_v1.tsv")
SUMMARY_TSV <- file.path(MODULE_ROOT, "tcga_source_module_summary_v1.tsv")
EIGENGENES_TSV <- file.path(MODULE_ROOT, "tcga_source_module_eigengenes_v1.tsv")
EIGENCOR_TSV <- file.path(MODULE_ROOT, "tcga_source_module_eigengene_correlations_v1.tsv")
TOM_RDATA <- file.path(MODULE_ROOT, "tcga_source_TOM_v1-block.1.RData")
CONTRACT_JSON <- file.path(CONTRACT_ROOT, "tcga_wgcna_module_contract_v2.json")

FINAL_JSON <- file.path(MODULE_ROOT, "tcga_source_modules_frozen_v1.json")
FINAL_MD <- file.path(MODULE_ROOT, "tcga_source_modules_frozen_v1.md")
MANIFEST_TSV <- file.path(MODULE_ROOT, "tcga_source_modules_sha256_manifest_v1.tsv")

cat(strrep("=", 124), "\n", sep = "")
cat("Paper 4 / TCBB - finalize already-computed frozen TCGA source modules after serialization-only failure\n")
cat(strrep("=", 124), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("\nProvenance:\n")
cat("  02e v1 completed WGCNA module detection and wrote all scientific artifacts.\n")
cat("  It then failed only while constructing the final R result list because of an R string-concatenation syntax error.\n")
cat("  This script DOES NOT recompute modules. It validates and finalizes the already-written first-run artifacts.\n")
cat("\nScientific guard:\n")
cat("  SCAN-B opened:                                      NO\n")
cat("  METABRIC opened:                                    NO\n")
cat("  GSE239948 opened:                                   NO\n")
cat("  Clinical outcomes loaded:                           NO\n")
cat("  Target preservation statistics calculated:          NO\n")
cat("  WGCNA modules recomputed here:                       NO\n")
cat(strrep("=", 124), "\n\n", sep = "")

required <- c(
  NETWORK_RDS, MEMBERSHIP_TSV, SUMMARY_TSV, EIGENGENES_TSV,
  EIGENCOR_TSV, TOM_RDATA, CONTRACT_JSON
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing partial-run artifact(s):\n", paste(missing, collapse = "\n"))
}

if (!requireNamespace("WGCNA", quietly = TRUE)) stop("WGCNA is not installed.")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed.")
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is not installed.")

suppressPackageStartupMessages(library(WGCNA))

cat("[1/4] Loading and validating first-run WGCNA object ...\n")
net <- readRDS(NETWORK_RDS)

if (length(net$colors) != 10000) {
  stop("Expected 10,000 module labels; found ", length(net$colors))
}
if (length(net$blocks) != 10000) {
  stop("Expected 10,000 block labels; found ", length(net$blocks))
}
if (length(unique(net$blocks)) != 1) {
  stop("Expected one WGCNA block; found ", length(unique(net$blocks)))
}

labels <- as.integer(net$colors)
colors <- WGCNA::labels2colors(labels)

non_grey_labels <- sort(unique(labels[labels != 0]))
n_modules <- length(non_grey_labels)
grey_n <- sum(labels == 0)
grey_fraction <- grey_n / length(labels)

if (n_modules != 12) {
  stop("Expected 12 non-grey modules from first run; found ", n_modules)
}
if (grey_n != 54) {
  stop("Expected 54 grey genes from first run; found ", grey_n)
}

cat("  non-grey modules: ", n_modules, "\n", sep = "")
cat("  grey genes:       ", grey_n, " / 10000 (",
    sprintf("%.2f", 100 * grey_fraction), "%)\n", sep = "")

cat("\n[2/4] Replaying saved membership/summary consistency checks ...\n")
membership <- read.delim(
  MEMBERSHIP_TSV, sep = "\t", header = TRUE, check.names = FALSE,
  quote = "", comment.char = ""
)
sizes <- read.delim(
  SUMMARY_TSV, sep = "\t", header = TRUE, check.names = FALSE,
  quote = "", comment.char = ""
)

if (nrow(membership) != 10000) stop("Membership table is not 10,000 rows.")
if (!identical(as.integer(membership$wgcna_label), labels)) {
  bad <- which(as.integer(membership$wgcna_label) != labels)
  stop("Saved membership does not replay network object. First mismatch row: ", bad[1])
}
if (!identical(as.character(membership$wgcna_color), as.character(colors))) {
  bad <- which(as.character(membership$wgcna_color) != as.character(colors))
  stop("Saved colors do not replay network object. First mismatch row: ", bad[1])
}

calc_sizes <- sort(table(labels), decreasing = FALSE)
saved_sizes <- sizes$n_genes[match(as.integer(names(calc_sizes)), sizes$wgcna_label)]
if (any(is.na(saved_sizes)) || !identical(as.integer(saved_sizes), as.integer(calc_sizes))) {
  stop("Saved module-size summary does not replay network object.")
}

MEs <- as.data.frame(net$MEs, check.names = FALSE)
me_cor <- cor(MEs, use = "pairwise.complete.obs", method = "pearson")
diag(me_cor) <- NA_real_
max_positive_me_cor <- if (ncol(MEs) > 1) max(me_cor, na.rm = TRUE) else NA_real_
min_me_cor <- if (ncol(MEs) > 1) min(me_cor, na.rm = TRUE) else NA_real_

cat("  eigengenes: ", ncol(MEs), "\n", sep = "")
cat("  max off-diagonal ME correlation: ",
    signif(max_positive_me_cor, 6), "\n", sep = "")

cat("\n[3/4] Writing final frozen result metadata ...\n")

non_grey_sizes <- sizes[sizes$wgcna_label != 0, , drop = FALSE]
module_size_list <- as.list(as.integer(non_grey_sizes$n_genes))
names(module_size_list) <- non_grey_sizes$program_id

result <- list(
  result_id = "paper4-tcbb-tcga-source-modules-v1",
  script_version_finalizer = SCRIPT_VERSION,
  original_compute_script = "02e-build-tcga-source-modules-v1-no-cli",
  status = "FROZEN_SOURCE_MODULES_FINALIZED_FROM_FIRST_RUN",
  provenance = list(
    first_run_scientific_computation_completed = TRUE,
    first_run_artifacts_written_before_failure = TRUE,
    first_run_failure_stage = "final result-list construction / serialization only",
    first_run_failure_reason = "R syntax error from adjacent string literals in freeze_rule",
    modules_recomputed_by_finalizer = FALSE
  ),
  scientific_guard = list(
    scanb_opened = FALSE,
    metabric_opened = FALSE,
    gse239948_opened = FALSE,
    clinical_outcomes_loaded = FALSE,
    target_coverage_or_preservation_used = FALSE,
    posthoc_module_parameter_tuning = FALSE
  ),
  source = list(
    samples = 1082L,
    genes = 10000L
  ),
  contract = CONTRACT_JSON,
  wgcna = list(
    version = as.character(packageVersion("WGCNA")),
    blocks = length(unique(net$blocks)),
    non_grey_modules = n_modules,
    grey_genes = grey_n,
    grey_fraction = grey_fraction,
    module_sizes = module_size_list
  ),
  source_eigengene_structure = list(
    n_eigengenes = ncol(MEs),
    max_off_diagonal_pearson_correlation = max_positive_me_cor,
    min_off_diagonal_pearson_correlation = min_me_cor
  ),
  freeze_rule = paste0(
    "Every non-grey module returned by the pre-frozen WGCNA contract is retained. ",
    "No module is selected or discarded based on size beyond the pre-frozen detection rule, ",
    "biology, target coverage, or target preservation."
  )
)

jsonlite::write_json(
  result,
  FINAL_JSON,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

writeLines(
  c(
    "# Paper 4 / TCBB — Frozen TCGA-BRCA Source Modules v1",
    "",
    "Status: **FROZEN AND FINALIZED FROM THE ORIGINAL FIRST MODULE-DETECTION RUN**",
    "",
    "The original `02e` run completed WGCNA module detection, validation, and artifact writing.",
    "It failed only afterward while constructing a final R metadata object because adjacent",
    "string literals were used instead of explicit concatenation. No scientific computation",
    "failed, and this finalizer does not recompute the modules.",
    "",
    paste0("- Source samples: **1082**"),
    paste0("- Source genes: **10000**"),
    paste0("- Non-grey modules: **", n_modules, "**"),
    paste0("- Grey genes: **", grey_n, " (", sprintf("%.2f", 100 * grey_fraction), "%)**"),
    paste0("- Source eigengenes returned: **", ncol(MEs), "**"),
    paste0("- Maximum off-diagonal source eigengene correlation: **",
           signif(max_positive_me_cor, 6), "**"),
    "",
    "All non-grey modules are retained automatically; no biological or target-based",
    "module selection has occurred."
  ),
  FINAL_MD
)

cat("\n[4/4] Creating SHA-256 manifest ...\n")
manifest_files <- c(
  CONTRACT_JSON,
  NETWORK_RDS,
  MEMBERSHIP_TSV,
  SUMMARY_TSV,
  EIGENGENES_TSV,
  EIGENCOR_TSV,
  TOM_RDATA,
  FINAL_JSON,
  FINAL_MD
)

sha <- vapply(
  manifest_files,
  function(f) digest::digest(file = f, algo = "sha256", serialize = FALSE),
  character(1)
)

manifest <- data.frame(
  file = manifest_files,
  bytes = as.numeric(file.info(manifest_files)$size),
  sha256 = sha,
  stringsAsFactors = FALSE
)
write.table(
  manifest,
  MANIFEST_TSV,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

cat("\n", strrep("=", 124), "\n", sep = "")
cat("02e1 FROZEN TCGA SOURCE MODULE FINALIZATION: PASS\n")
cat(strrep("=", 124), "\n", sep = "")
cat("Scientific modules recomputed: NO\n")
cat("First-run non-grey modules:    ", n_modules, "\n", sep = "")
cat("First-run grey genes:          ", grey_n, " (",
    sprintf("%.2f", 100 * grey_fraction), "%)\n", sep = "")
cat("First-run eigengenes:          ", ncol(MEs), "\n", sep = "")
cat("Max source ME correlation:     ", signif(max_positive_me_cor, 6), "\n", sep = "")
cat("\nFinal outputs:\n")
cat("  ", FINAL_JSON, "\n", sep = "")
cat("  ", FINAL_MD, "\n", sep = "")
cat("  ", MANIFEST_TSV, "\n", sep = "")
cat("\nNo target dataset was opened.\n")
cat("No target preservation statistic was calculated.\n")
cat(strrep("=", 124), "\n", sep = "")
