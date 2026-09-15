options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "05i5a-freeze-partial-comparator-continuation-v1-no-cli"
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

V4_ROOT <- file.path(
  DATA_ROOT,
  "paper4_tcbb_postprimary_comparator_benchmark_v4"
)
SCANB_V4_DIR <- file.path(V4_ROOT, "SCANB_GSE96058")
SCANB_NETREP_RDS <- file.path(
  SCANB_V4_DIR,
  "netrep_modulePreservation_v4.rds"
)
SCANB_NETREP_TSV <- file.path(
  SCANB_V4_DIR,
  "netrep_modulePreservation_long_v4.tsv"
)

FORBIDDEN_ALREADY_CALCULATED <- c(
  file.path(SCANB_V4_DIR, "wgcna_modulePreservation_v4.rds"),
  file.path(SCANB_V4_DIR, "wgcna_modulePreservation_summary_v4.tsv"),
  file.path(V4_ROOT, "METABRIC", "netrep_modulePreservation_v4.rds"),
  file.path(V4_ROOT, "METABRIC", "netrep_modulePreservation_long_v4.tsv"),
  file.path(V4_ROOT, "METABRIC", "wgcna_modulePreservation_v4.rds"),
  file.path(V4_ROOT, "METABRIC", "wgcna_modulePreservation_summary_v4.tsv")
)

OUT_DIR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_partial_comparator_continuation_contract_v1"
)
OUT_JSON <- file.path(
  OUT_DIR,
  "partial_comparator_continuation_contract_v1.json"
)

NETREP_THREADS_FOR_REMAINING_RUNS <- 8L
NETREP_PERM <- 10000L
WGCNA_PERM <- 200L

cat(strrep("=", 156), "\n", sep = "")
cat("Paper 4 / TCBB - freeze continuation after partial comparator completion\n")
cat(strrep("=", 156), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n\n", sep = "")

cat("Provenance situation:\n")
cat("  SCAN-B NetRep v4 artifacts exist:                   YES\n")
cat("  Their numerical contents inspected here:             NO\n")
cat("  They are frozen by file hash before further work:     YES\n")
cat("  SCAN-B NetRep will be recomputed:                     NO\n")
cat("  Remaining NetRep work may use 8 threads:              YES\n")
cat("  Scientific settings/interpretation changed:           NO\n")
cat(strrep("=", 156), "\n\n", sep = "")

required <- c(
  BASELINE_CONTRACT,
  INTERPRETATION_CONTRACT,
  PERM_AMENDMENT,
  SCANB_NETREP_RDS,
  SCANB_NETREP_TSV
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required file(s):\n", paste(missing, collapse = "\n"))
}

bad <- FORBIDDEN_ALREADY_CALCULATED[file.exists(FORBIDDEN_ALREADY_CALCULATED)]
if (length(bad) > 0) {
  stop(
    "Additional comparator artifacts already exist, so this exact partial-continuation contract is not valid:\n",
    paste(bad, collapse = "\n"),
    "\nStop and review provenance before continuing."
  )
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required.")
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("R package 'digest' is required for SHA-256 file freezing.")
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

# Hash exact bytes only; do not parse scientific values.
rds_sha256 <- digest::digest(
  file = SCANB_NETREP_RDS,
  algo = "sha256",
  serialize = FALSE
)
tsv_sha256 <- digest::digest(
  file = SCANB_NETREP_TSV,
  algo = "sha256",
  serialize = FALSE
)

rds_info <- file.info(SCANB_NETREP_RDS)
tsv_info <- file.info(SCANB_NETREP_TSV)

contract <- list(
  contract_id =
    "paper4-tcbb-partial-comparator-continuation-v1",
  script_version = SCRIPT_VERSION,
  status =
    "FROZEN_AFTER_SCANB_NETREP_BYTES_EXIST_BEFORE_NUMERICAL_INSPECTION",

  provenance = list(
    explanation = paste(
      "05i v4 completed and serialized the SCAN-B NetRep comparator",
      "before the attempted threading amendment. The existence guard",
      "correctly prevented pretending that no comparator result had",
      "been calculated. This continuation freezes the exact completed",
      "SCAN-B NetRep bytes without parsing their numerical contents."
    ),
    scanb_netrep_numerical_contents_inspected_by_this_script = FALSE,
    scanb_netrep_recomputed = FALSE,
    scientific_interpretation_changed = FALSE,
    baseline_comparator_settings_changed = FALSE
  ),

  authoritative_completed_artifact = list(
    target = "SCANB_GSE96058",
    method = "NetRep",
    source_runner =
      "05i-run-postprimary-comparator-benchmark-v4-no-cli",
    nPerm = NETREP_PERM,
    nThreads = 1L,
    rds = list(
      path = SCANB_NETREP_RDS,
      sha256 = rds_sha256,
      bytes = unname(as.numeric(rds_info$size))
    ),
    tsv = list(
      path = SCANB_NETREP_TSV,
      sha256 = tsv_sha256,
      bytes = unname(as.numeric(tsv_info$size))
    ),
    rule =
      "These exact bytes are authoritative and must be reused, never recomputed."
  ),

  remaining_execution = list(
    scanb_wgcna = list(
      required = TRUE,
      nPermutations = WGCNA_PERM,
      settings_changed = FALSE
    ),
    metabric_netrep = list(
      required = TRUE,
      nPerm = NETREP_PERM,
      nThreads = NETREP_THREADS_FOR_REMAINING_RUNS,
      reason_for_threading =
        "runtime optimization decided without opening SCAN-B NetRep numerical values",
      scientific_settings_changed = FALSE
    ),
    metabric_wgcna = list(
      required = TRUE,
      nPermutations = WGCNA_PERM,
      settings_changed = FALSE
    )
  ),

  output_rule =
    "Continue in a new v6 result directory; copy/reuse exact SCAN-B NetRep v4 artifacts after verifying hashes.",

  primary_mta_classification_changed = FALSE
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  contract,
  OUT_JSON,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

cat(strrep("=", 156), "\n", sep = "")
cat("05i5a PARTIAL COMPARATOR CONTINUATION CONTRACT: PASS\n")
cat(strrep("=", 156), "\n", sep = "")
cat("Authoritative completed result:       SCAN-B NetRep v4\n")
cat("SCAN-B NetRep recomputation:          FORBIDDEN\n")
cat("Frozen RDS SHA-256:                   ", rds_sha256, "\n", sep = "")
cat("Frozen TSV SHA-256:                   ", tsv_sha256, "\n", sep = "")
cat("Remaining NetRep threads:             8\n")
cat("Remaining NetRep permutations:        10,000\n")
cat("WGCNA permutations:                   200\n")
cat("Scientific interpretation changed:    NO\n")
cat("\nOutput: ", OUT_JSON, "\n", sep = "")
cat(strrep("=", 156), "\n", sep = "")
