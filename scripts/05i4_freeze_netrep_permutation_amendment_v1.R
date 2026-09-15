options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "05i4-freeze-netrep-permutation-amendment-v1-no-cli"
DATA_ROOT <- "D:/paper4_tcbb_data"

BASELINE_CONTRACT <- file.path(
  DATA_ROOT,
  "paper4_tcbb_postprimary_comparator_benchmark_contract_v1",
  "postprimary_comparator_benchmark_contract_v1.json"
)
INTERPRETATION_CONTRACT <- file.path(
  DATA_ROOT,
  "paper4_tcbb_comparator_interpretation_headtohead_contract_v1",
  "comparator_interpretation_headtohead_contract_v1.json"
)
OUT_DIR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_netrep_permutation_amendment_v1"
)
OUT_JSON <- file.path(
  OUT_DIR,
  "netrep_permutation_amendment_v1.json"
)

POSSIBLE_RESULT_DIRS <- file.path(
  DATA_ROOT,
  c(
    "paper4_tcbb_postprimary_comparator_benchmark_v1",
    "paper4_tcbb_postprimary_comparator_benchmark_v2",
    "paper4_tcbb_postprimary_comparator_benchmark_v3",
    "paper4_tcbb_postprimary_comparator_benchmark_v4"
  )
)

ORIGINAL_NETREP_PERM <- 1000L
AMENDED_NETREP_PERM <- 10000L
WGCNA_PERM_UNCHANGED <- 200L

cat(strrep("=", 154), "\n", sep = "")
cat("Paper 4 / TCBB - freeze pre-result NetRep permutation-count amendment\n")
cat(strrep("=", 154), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n\n", sep = "")

cat("Reason for amendment:\n")
cat("  NetRep documentation recommends at least 10,000 permutations so the\n")
cat("  permutation null distributions are representative. The earlier 05i2\n")
cat("  contract froze 1,000. This amendment is made BEFORE any NetRep/WGCNA\n")
cat("  scientific comparator result has been calculated or inspected.\n\n")

cat("Scientific guard:\n")
cat("  Comparator results observed:                          NO\n")
cat("  NetRep nPerm changed:                                 1,000 -> 10,000\n")
cat("  WGCNA nPermutations changed:                          NO (remains 200)\n")
cat("  Comparator statistics/null/interpretation changed:   NO\n")
cat("  Primary MTA classification changed:                   NO\n")
cat(strrep("=", 154), "\n\n", sep = "")

for (f in c(BASELINE_CONTRACT, INTERPRETATION_CONTRACT)) {
  if (!file.exists(f)) stop("Missing required contract: ", f)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required.")
}

base <- jsonlite::read_json(BASELINE_CONTRACT, simplifyVector = TRUE)
interp <- jsonlite::read_json(INTERPRETATION_CONTRACT, simplifyVector = TRUE)

if (as.character(base$status) !=
    "FROZEN_POST_PRIMARY_BEFORE_FIRST_COMPARATOR_STATISTIC") {
  stop("05i2 baseline comparator contract has unexpected status.")
}
if (as.character(interp$status) !=
    "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT") {
  stop("05i3 interpretation contract has unexpected status.")
}
if (as.integer(base$NetRep$nPerm) != ORIGINAL_NETREP_PERM) {
  stop("05i2 NetRep nPerm is not 1000; review provenance.")
}

scientific_patterns <- c(
  "netrep_modulePreservation.*\\.rds$",
  "netrep_modulePreservation.*\\.tsv$",
  "wgcna_modulePreservation.*\\.rds$",
  "wgcna_modulePreservation.*\\.tsv$",
  "integrated_mta_netrep_wgcna.*\\.tsv$",
  "netrep_all_targets.*\\.tsv$",
  "wgcna_modulePreservation_all_targets.*\\.tsv$",
  "mta_netrep_wgcna_integrated_all_targets.*\\.tsv$"
)

found <- character(0)
for (d in POSSIBLE_RESULT_DIRS) {
  if (!dir.exists(d)) next
  ff <- list.files(d, recursive = TRUE, full.names = TRUE)
  if (length(ff) == 0) next
  ff <- ff[file.info(ff)$isdir %in% FALSE]
  if (length(ff) == 0) next
  bn <- basename(ff)
  keep <- rep(FALSE, length(ff))
  for (pat in scientific_patterns) keep <- keep | grepl(pat, bn, ignore.case = TRUE)
  found <- c(found, ff[keep])
}
if (length(found) > 0) {
  stop(
    "Scientific comparator artifacts already exist before 05i4 amendment:\n",
    paste(found, collapse = "\n")
  )
}

amendment <- list(
  contract_id = "paper4-tcbb-netrep-permutation-amendment-v1",
  script_version = SCRIPT_VERSION,
  status = "FROZEN_BEFORE_FIRST_COMPARATOR_STATISTIC",
  supersedes_for_netrep_permutation_count = "05i2 nPerm=1000 only",
  provenance = list(
    comparator_results_observed_before_amendment = FALSE,
    reason = paste(
      "NetRep documentation recommends at least 10,000 permutations",
      "for representative null distributions."
    ),
    target_sets_changed = FALSE,
    module_sets_changed = FALSE,
    statistics_changed = FALSE,
    null_model_changed = FALSE,
    alternative_changed = FALSE,
    interpretation_rules_changed = FALSE
  ),
  NetRep = list(
    previous_nPerm = ORIGINAL_NETREP_PERM,
    amended_nPerm = AMENDED_NETREP_PERM,
    null = "overlap",
    alternative = "greater",
    all_seven_statistics_retained = TRUE,
    nThreads = 1L
  ),
  WGCNA = list(
    nPermutations = WGCNA_PERM_UNCHANGED,
    changed_by_this_amendment = FALSE
  ),
  primary_mta_classification_changed = FALSE
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  amendment,
  OUT_JSON,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

cat(strrep("=", 154), "\n", sep = "")
cat("05i4 NETREP PERMUTATION AMENDMENT: PASS\n")
cat(strrep("=", 154), "\n", sep = "")
cat("NetRep nPerm:                       10,000\n")
cat("WGCNA nPermutations:                200 (unchanged)\n")
cat("Comparator results seen beforehand:NO\n")
cat("Scientific comparator estimand:     unchanged\n")
cat("Primary MTA classification changed: NO\n")
cat("\nOutput: ", OUT_JSON, "\n", sep = "")
cat(strrep("=", 154), "\n", sep = "")
