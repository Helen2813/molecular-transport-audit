options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "02f-freeze-tcga-source-program-contract-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"
MODULE_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_tcga_source_modules_v1")

MODULE_JSON <- file.path(MODULE_ROOT, "tcga_source_modules_frozen_v1.json")
MEMBERSHIP_TSV <- file.path(MODULE_ROOT, "tcga_source_module_membership_frozen_v1.tsv")
SOURCE_MATRIX <- file.path(
  DATA_ROOT,
  "paper4_tcbb_wgcna_soft_threshold_v1",
  "tcga_source_datExpr_1082x10000_v1.rds"
)

OUT_DIR <- file.path(DATA_ROOT, "paper4_tcbb_source_program_contract_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# FROZEN REPRESENTATION CONTRACT
# =============================================================================

SOURCE_EXPRESSION <- "log2(RSEM + 1)"
WITHIN_COHORT_GENE_STANDARDIZATION <- "z-score across source samples using sample SD"
SOURCE_EDGE_CORRELATION <- "Pearson"
EDGE_VECTOR_ORDER <- "upper triangle in frozen source-gene order"

WEIGHT_DEFINITION <- "first right singular vector (PC1 loading) of standardized source module"
WEIGHT_NORM <- "unit L2 norm as returned by SVD"

ORIENTATION_REFERENCE <- "unweighted mean of standardized module genes across each source sample"
ORIENTATION_TOL <- 1e-12
ORIENTATION_FALLBACK <- paste0(
  "if abs(cor(PC1 score, unweighted standardized mean)) <= 1e-12, ",
  "orient so the lexicographically first gene among maximum-absolute-loading genes ",
  "has positive loading"
)

TARGET_ORIENTATION_SCORE <- paste0(
  "for later target use only: mean(sign(frozen source loading_g) * target z_g) ",
  "over evaluable frozen genes; target PC1 sign is flipped if its correlation ",
  "with this frozen signed score is negative"
)

SOURCE_COHERENCE_PRIMARY <- "mean absolute Pearson correlation across nonredundant within-program edges"

cat(strrep("=", 122), "\n", sep = "")
cat("Paper 4 / TCBB - freeze TCGA source-program representation contract\n")
cat(strrep("=", 122), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("\nScientific guard:\n")
cat("  SCAN-B opened:                                      NO\n")
cat("  METABRIC opened:                                    NO\n")
cat("  GSE239948 opened:                                   NO\n")
cat("  Clinical outcomes loaded:                           NO\n")
cat("  Target preservation statistics calculated:          NO\n")
cat("  Source PC1/loadings calculated by this script:       NO\n")
cat("  Operation: freeze representation definitions only\n")
cat(strrep("=", 122), "\n\n", sep = "")

for (f in c(MODULE_JSON, MEMBERSHIP_TSV, SOURCE_MATRIX)) {
  if (!file.exists(f)) stop("Missing required frozen source artifact: ", f)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed.")

module_meta <- jsonlite::read_json(MODULE_JSON, simplifyVector = TRUE)
if (as.integer(module_meta$wgcna$non_grey_modules) != 12L) {
  stop("Expected 12 frozen non-grey modules.")
}

membership <- read.delim(
  MEMBERSHIP_TSV, sep = "\t", header = TRUE, check.names = FALSE,
  quote = "", comment.char = ""
)
program_ids <- sort(unique(membership$program_id[membership$program_id != "GREY"]))
if (length(program_ids) != 12L) stop("Membership does not contain exactly 12 non-grey programs.")

datExpr <- readRDS(SOURCE_MATRIX)
if (nrow(datExpr) != 1082 || ncol(datExpr) != 10000) {
  stop("Unexpected frozen source matrix dimensions.")
}
rm(datExpr)
gc(verbose = FALSE)

contract <- list(
  contract_id = "paper4-tcbb-source-program-representation-v1",
  script_version = SCRIPT_VERSION,
  status = "FROZEN_BEFORE_SOURCE_PROGRAM_STATISTICS",
  scientific_guard = list(
    scanb_opened = FALSE,
    metabric_opened = FALSE,
    gse239948_opened = FALSE,
    clinical_outcomes_loaded = FALSE,
    target_preservation_statistics_calculated = FALSE,
    source_pc1_or_edge_statistics_calculated_by_this_script = FALSE
  ),
  frozen_programs = program_ids,
  representation = list(
    source_expression = SOURCE_EXPRESSION,
    within_cohort_gene_standardization = WITHIN_COHORT_GENE_STANDARDIZATION,
    source_weight_definition = WEIGHT_DEFINITION,
    source_weight_norm = WEIGHT_NORM,
    source_pc1_orientation_reference = ORIENTATION_REFERENCE,
    source_pc1_orientation_tolerance = ORIENTATION_TOL,
    source_pc1_orientation_fallback = ORIENTATION_FALLBACK,
    later_target_orientation_score = TARGET_ORIENTATION_SCORE
  ),
  source_structure = list(
    gene_gene_correlation = SOURCE_EDGE_CORRELATION,
    edge_vector_order = EDGE_VECTOR_ORDER,
    primary_coherence_metric_for_later_null_matching = SOURCE_COHERENCE_PRIMARY,
    additional_frozen_summaries = c(
      "median absolute edge correlation",
      "mean signed edge correlation",
      "SD of signed edge correlation",
      "5th/25th/50th/75th/95th percentiles",
      "fraction of positive edges",
      "PC1 variance explained"
    )
  ),
  persistence = list(
    freeze_gene_membership = TRUE,
    freeze_oriented_source_loadings = TRUE,
    freeze_source_edge_vectors = TRUE,
    freeze_source_pc1_scores = TRUE,
    freeze_source_signed_scores = TRUE
  ),
  later_null_constraint = paste0(
    "The primary source-coherence matching statistic is fixed here as mean absolute ",
    "source Pearson edge correlation. The exact panel-matching tolerance/algorithm ",
    "will be frozen separately using source-only information before any target ",
    "preservation statistic is calculated."
  )
)

json_out <- file.path(OUT_DIR, "tcga_source_program_representation_contract_v1.json")
jsonlite::write_json(contract, json_out, pretty = TRUE, auto_unbox = TRUE, digits = NA)

md_out <- file.path(OUT_DIR, "tcga_source_program_representation_contract_v1.md")
writeLines(
  c(
    "# Paper 4 / TCBB — TCGA Source-Program Representation Contract v1",
    "",
    "Status: **FROZEN BEFORE SOURCE PROGRAM STATISTICS**",
    "",
    "No target dataset or clinical outcome was used.",
    "",
    "## Frozen source weights",
    "",
    "- Within each frozen non-grey WGCNA module, source expression is standardized gene-wise across the 1,082 source samples.",
    "- Frozen weights are the PC1 loadings (first right singular vector) of that standardized source module.",
    "- The loading vector retains the SVD unit-L2 normalization.",
    "- Source PC1 sign is oriented to correlate positively with the unweighted mean standardized module score.",
    paste0("- Orientation tie tolerance: `", ORIENTATION_TOL, "`."),
    "- If the orientation correlation is within tolerance of zero, the maximum-absolute-loading gene provides a deterministic sign fallback.",
    "",
    "## Frozen source edge structure",
    "",
    "- Within-program gene-gene correlation: **Pearson**.",
    "- Edge vector: nonredundant upper triangle in frozen source-gene order.",
    "- Primary source-coherence quantity for later null matching: **mean absolute Pearson edge correlation**.",
    "- Additional summaries are retained but do not replace the primary coherence definition.",
    "",
    "## Later target orientation",
    "",
    "Target PC1 orientation remains outcome-blind and uses the frozen sign-score:",
    "`mean(sign(w_g) * Z_ig)` over evaluable frozen genes.",
    "",
    "The exact random-panel matching tolerance and algorithm are not chosen here;",
    "they will be frozen from source-only information before target preservation is computed."
  ),
  md_out
)

cat(strrep("=", 122), "\n", sep = "")
cat("02f SOURCE-PROGRAM REPRESENTATION CONTRACT: PASS\n")
cat(strrep("=", 122), "\n", sep = "")
cat("Frozen non-grey programs:     ", length(program_ids), "\n", sep = "")
cat("Frozen source weights:        PC1 loadings on gene-standardized source expression\n")
cat("Frozen source PC1 orientation:positive correlation with unweighted standardized module mean\n")
cat("Frozen source edge type:      Pearson\n")
cat("Frozen coherence metric:      mean absolute source edge correlation\n")
cat("Frozen target orientation:    sign(frozen source weights) score\n")
cat("\nNo source PC1/loadings were calculated by this script.\n")
cat("No target dataset was opened.\n")
cat("\nOutputs:\n")
cat("  ", json_out, "\n", sep = "")
cat("  ", md_out, "\n", sep = "")
cat(strrep("=", 122), "\n", sep = "")
