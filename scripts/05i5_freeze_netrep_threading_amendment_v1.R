options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "05i5-freeze-netrep-threading-amendment-v1-no-cli"
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
PERM_AMENDMENT <- file.path(
  DATA_ROOT,
  "paper4_tcbb_netrep_permutation_amendment_v1",
  "netrep_permutation_amendment_v1.json"
)

OUT_DIR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_netrep_threading_amendment_v1"
)
OUT_JSON <- file.path(
  OUT_DIR,
  "netrep_threading_amendment_v1.json"
)

POSSIBLE_RESULT_DIRS <- file.path(
  DATA_ROOT,
  c(
    "paper4_tcbb_postprimary_comparator_benchmark_v1",
    "paper4_tcbb_postprimary_comparator_benchmark_v2",
    "paper4_tcbb_postprimary_comparator_benchmark_v3",
    "paper4_tcbb_postprimary_comparator_benchmark_v4",
    "paper4_tcbb_postprimary_comparator_benchmark_v5"
  )
)

PREVIOUS_THREADS <- 1L
AMENDED_THREADS <- 8L
NETREP_PERM <- 10000L
WGCNA_PERM <- 200L

cat(strrep("=", 154), "\n", sep = "")
cat("Paper 4 / TCBB - freeze pre-result NetRep threading amendment\n")
cat(strrep("=", 154), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n\n", sep = "")

cat("Reason for amendment:\n")
cat("  The 10,000-permutation NetRep run was observed to be computationally\n")
cat("  impractical at nThreads=1. The run was stopped before completion and\n")
cat("  before any NetRep/WGCNA scientific comparator result was produced or\n")
cat("  inspected. This amendment changes execution parallelism only.\n\n")

cat("Scientific guard:\n")
cat("  Comparator scientific results observed:             NO\n")
cat("  NetRep permutations changed:                         NO (10,000)\n")
cat("  NetRep null/statistics/modules changed:              NO\n")
cat("  NetRep nThreads changed:                             1 -> 8\n")
cat("  WGCNA settings changed:                              NO\n")
cat("  Interpretation contract changed:                     NO\n")
cat("  Primary MTA classification changed:                  NO\n")
cat(strrep("=", 154), "\n\n", sep = "")

for (f in c(BASELINE_CONTRACT, INTERPRETATION_CONTRACT, PERM_AMENDMENT)) {
  if (!file.exists(f)) stop("Missing required contract/amendment: ", f)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required.")
}

base <- jsonlite::read_json(BASELINE_CONTRACT, simplifyVector = TRUE)
interp <- jsonlite::read_json(INTERPRETATION_CONTRACT, simplifyVector = TRUE)
perm <- jsonlite::read_json(PERM_AMENDMENT, simplifyVector = TRUE)

if (as.character(base$status) !=
    "FROZEN_POST_PRIMARY_BEFORE_FIRST_COMPARATOR_STATISTIC") {
  stop("05i2 baseline comparator contract has unexpected status.")
}
if (as.character(interp$status) !=
    "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT") {
  stop("05i3 interpretation contract has unexpected status.")
}
if (as.character(perm$status) !=
    "FROZEN_BEFORE_FIRST_COMPARATOR_STATISTIC") {
  stop("05i4 permutation amendment has unexpected status.")
}
if (as.integer(perm$NetRep$amended_nPerm) != NETREP_PERM) {
  stop("05i4 NetRep permutation count is not 10,000.")
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
  for (pat in scientific_patterns) {
    keep <- keep | grepl(pat, bn, ignore.case = TRUE)
  }
  found <- c(found, ff[keep])
}

if (length(found) > 0) {
  stop(
    "Scientific comparator artifacts already exist before 05i5 threading amendment:\n",
    paste(found, collapse = "\n"),
    "\nDo not continue until provenance is reviewed."
  )
}

amendment <- list(
  contract_id = "paper4-tcbb-netrep-threading-amendment-v1",
  script_version = SCRIPT_VERSION,
  status = "FROZEN_BEFORE_FIRST_COMPARATOR_STATISTIC",

  supersedes_for_execution_parallelism =
    "05i2/05i4 NetRep nThreads=1 only",

  provenance = list(
    comparator_results_observed_before_amendment = FALSE,
    aborted_single_thread_run_completed = FALSE,
    reason = paste(
      "The 10,000-permutation single-thread run was computationally",
      "impractical and was stopped before any comparator result was",
      "produced or inspected."
    ),
    scientific_estimand_changed = FALSE,
    target_sets_changed = FALSE,
    module_sets_changed = FALSE,
    statistics_changed = FALSE,
    null_model_changed = FALSE,
    permutation_count_changed = FALSE,
    interpretation_rules_changed = FALSE
  ),

  NetRep = list(
    nPerm = NETREP_PERM,
    previous_nThreads = PREVIOUS_THREADS,
    amended_nThreads = AMENDED_THREADS,
    null = "overlap",
    alternative = "greater",
    all_seven_statistics_retained = TRUE,
    reproducibility_rule = paste(
      "The v5 run fixes nThreads=8 and the already-frozen target-specific",
      "set.seed schedule. The aborted nThreads=1 permutation stream is not",
      "used as a scientific reference because it produced no completed result."
    )
  ),

  WGCNA = list(
    nPermutations = WGCNA_PERM,
    parallelCalculation = FALSE,
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
cat("05i5 NETREP THREADING AMENDMENT: PASS\n")
cat(strrep("=", 154), "\n", sep = "")
cat("NetRep nPerm:                       10,000\n")
cat("NetRep nThreads:                    8\n")
cat("WGCNA nPermutations:                200 (unchanged)\n")
cat("Comparator results seen beforehand:NO\n")
cat("Scientific comparator estimand:     unchanged\n")
cat("Primary MTA classification changed: NO\n")
cat("\nOutput: ", OUT_JSON, "\n", sep = "")
cat(strrep("=", 154), "\n", sep = "")
