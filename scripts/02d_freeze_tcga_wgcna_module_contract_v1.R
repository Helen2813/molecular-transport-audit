options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "02d-freeze-tcga-wgcna-module-contract-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"
SFT_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_wgcna_soft_threshold_v1")
POWER_FILE <- file.path(SFT_ROOT, "tcga_wgcna_soft_power_frozen_v1.tsv")
SOURCE_MATRIX <- file.path(SFT_ROOT, "tcga_source_datExpr_1082x10000_v1.rds")

OUT_DIR <- file.path(DATA_ROOT, "paper4_tcbb_wgcna_module_contract_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PREDECLARED MODULE-CONSTRUCTION CONTRACT
# =============================================================================

EXPECTED_POWER <- 8

NETWORK_TYPE <- "signed"
COR_TYPE <- "pearson"
TOM_TYPE <- "signed"
TOM_DENOM <- "min"

# Force the 10,000-gene network into one block. This avoids projective-k-means
# block assignment and makes the source module definition easier to reproduce.
MAX_BLOCK_SIZE <- 12000L
RANDOM_SEED <- 54321L

# Dynamic Tree Cut / module merge settings.
DEEP_SPLIT <- 2L
DETECT_CUT_HEIGHT <- 0.995
MIN_MODULE_SIZE <- 30L
PAM_STAGE <- TRUE
PAM_RESPECTS_DENDRO <- FALSE

# Disable post-hoc kME reassignment and make trimming maximally permissive.
REASSIGN_THRESHOLD <- 0
MIN_KME_TO_STAY <- 0

# Merge modules whose eigengene dissimilarity is <= 0.25,
# i.e. eigengene correlation >= 0.75.
MERGE_CUT_HEIGHT <- 0.25

NUMERIC_LABELS <- TRUE
SAVE_TOMS <- TRUE
QUICK_COR <- 0
USE_COR_OPTIONS_THROUGHOUT <- TRUE

# MTA assessability guard, frozen before any target preservation analysis.
PRIMARY_MTA_MIN_EVALUABLE_GENES <- 30L
PRIMARY_MTA_MIN_TARGET_COVERAGE <- 0.80

cat(strrep("=", 120), "\n", sep = "")
cat("Paper 4 / TCBB - freeze TCGA-BRCA WGCNA module-construction contract\n")
cat(strrep("=", 120), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("\nScientific guard:\n")
cat("  SCAN-B opened:                                      NO\n")
cat("  METABRIC opened:                                    NO\n")
cat("  GSE239948 opened:                                   NO\n")
cat("  Clinical outcomes loaded:                           NO\n")
cat("  Target coverage/preservation inspected:             NO\n")
cat("  WGCNA modules calculated by this script:            NO\n")
cat("  Operation: freeze module-construction parameters only\n")
cat(strrep("=", 120), "\n\n", sep = "")

required <- c(POWER_FILE, SOURCE_MATRIX)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required files:\n", paste(missing, collapse = "\n"))
}

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("WGCNA is not installed.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is not installed.")
}

power_tbl <- read.delim(
  POWER_FILE,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

if (nrow(power_tbl) != 1) {
  stop("Expected exactly one frozen soft-power record.")
}

frozen_power <- as.numeric(power_tbl$chosen_power[1])

if (!isTRUE(all.equal(frozen_power, EXPECTED_POWER))) {
  stop(
    "Frozen soft power does not match predeclared module contract. Expected ",
    EXPECTED_POWER, ", got ", frozen_power
  )
}

# Verify the frozen matrix dimensions without calculating any network/module.
datExpr <- readRDS(SOURCE_MATRIX)
if (!is.matrix(datExpr) && !is.data.frame(datExpr)) {
  stop("Frozen source datExpr is not a matrix/data.frame.")
}
if (nrow(datExpr) != 1082 || ncol(datExpr) != 10000) {
  stop(
    "Unexpected frozen source matrix dimensions: ",
    nrow(datExpr), " x ", ncol(datExpr)
  )
}
if (any(!is.finite(as.matrix(datExpr)))) {
  stop("Frozen source matrix contains non-finite values.")
}
rm(datExpr)
gc(verbose = FALSE)

contract <- list(
  contract_id = "paper4-tcbb-wgcna-module-contract-v1",
  script_version = SCRIPT_VERSION,
  status = "FROZEN_BEFORE_MODULE_DETECTION",
  scientific_guard = list(
    scanb_opened = FALSE,
    metabric_opened = FALSE,
    gse239948_opened = FALSE,
    clinical_outcomes_loaded = FALSE,
    target_coverage_or_preservation_inspected = FALSE,
    modules_calculated_by_this_script = FALSE
  ),
  source_input = list(
    frozen_matrix = SOURCE_MATRIX,
    samples = 1082L,
    genes = 10000L
  ),
  network = list(
    power = frozen_power,
    networkType = NETWORK_TYPE,
    corType = COR_TYPE,
    TOMType = TOM_TYPE,
    TOMDenom = TOM_DENOM,
    maxBlockSize = MAX_BLOCK_SIZE,
    single_block_required = TRUE,
    randomSeed = RANDOM_SEED,
    quickCor = QUICK_COR,
    useCorOptionsThroughout = USE_COR_OPTIONS_THROUGHOUT
  ),
  module_detection = list(
    method = "WGCNA::blockwiseModules",
    deepSplit = DEEP_SPLIT,
    detectCutHeight = DETECT_CUT_HEIGHT,
    minModuleSize = MIN_MODULE_SIZE,
    pamStage = PAM_STAGE,
    pamRespectsDendro = PAM_RESPECTS_DENDRO,
    reassignThreshold = REASSIGN_THRESHOLD,
    minKMEtoStay = MIN_KME_TO_STAY,
    mergeCutHeight = MERGE_CUT_HEIGHT,
    numericLabels = NUMERIC_LABELS,
    saveTOMs = SAVE_TOMS
  ),
  assessability = list(
    minimum_evaluable_genes = PRIMARY_MTA_MIN_EVALUABLE_GENES,
    minimum_target_gene_coverage = PRIMARY_MTA_MIN_TARGET_COVERAGE,
    rule = paste0(
      "A module is primary-assessable in a target only if at least ",
      PRIMARY_MTA_MIN_EVALUABLE_GENES,
      " frozen genes are evaluable and evaluable coverage is at least ",
      100 * PRIMARY_MTA_MIN_TARGET_COVERAGE,
      "%. Otherwise it is reported descriptively as not assessable, not as failed."
    )
  ),
  software = list(
    R = R.version.string,
    WGCNA = as.character(packageVersion("WGCNA"))
  )
)

json_out <- file.path(OUT_DIR, "tcga_wgcna_module_contract_v1.json")
jsonlite::write_json(
  contract,
  path = json_out,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

md_out <- file.path(OUT_DIR, "tcga_wgcna_module_contract_v1.md")
writeLines(
  c(
    "# Paper 4 / TCBB — TCGA-BRCA WGCNA Module-Construction Contract v1",
    "",
    "Status: **FROZEN BEFORE MODULE DETECTION**",
    "",
    "No SCAN-B, METABRIC, GSE239948, target coverage, target preservation",
    "result, or clinical outcome was inspected.",
    "",
    "## Frozen network definition",
    "",
    paste0("- Soft power: **", frozen_power, "**"),
    "- Network type: **signed**",
    "- Correlation: **Pearson**",
    "- TOM type: **signed**",
    "- TOM denominator: **min**",
    paste0("- Maximum block size: **", MAX_BLOCK_SIZE,
           "** (forces the 10,000 genes into one block)"),
    paste0("- Random seed: **", RANDOM_SEED, "**"),
    "",
    "## Frozen module detection",
    "",
    paste0("- `deepSplit = ", DEEP_SPLIT, "`"),
    paste0("- `detectCutHeight = ", DETECT_CUT_HEIGHT, "`"),
    paste0("- `minModuleSize = ", MIN_MODULE_SIZE, "`"),
    paste0("- `pamStage = ", PAM_STAGE, "`"),
    paste0("- `pamRespectsDendro = ", PAM_RESPECTS_DENDRO, "`"),
    paste0("- `reassignThreshold = ", REASSIGN_THRESHOLD, "`"),
    paste0("- `minKMEtoStay = ", MIN_KME_TO_STAY, "`"),
    paste0("- `mergeCutHeight = ", MERGE_CUT_HEIGHT,
           "` (merge when eigengene correlation is at least 0.75)"),
    "- Numeric labels are used.",
    "- TOM is saved for exact downstream replay.",
    "",
    "## Frozen MTA assessability guard",
    "",
    paste0(
      "A target-module result is primary-assessable only when at least **",
      PRIMARY_MTA_MIN_EVALUABLE_GENES,
      " genes** remain and frozen-gene coverage is at least **",
      100 * PRIMARY_MTA_MIN_TARGET_COVERAGE,
      "%**."
    ),
    "",
    "Modules below either threshold will be labelled **not assessable** and",
    "reported descriptively; they will not be interpreted as negative preservation."
  ),
  con = md_out
)

# Machine-readable single-row parameter table.
params_out <- file.path(OUT_DIR, "tcga_wgcna_module_parameters_frozen_v1.tsv")
params <- data.frame(
  power = frozen_power,
  networkType = NETWORK_TYPE,
  corType = COR_TYPE,
  TOMType = TOM_TYPE,
  TOMDenom = TOM_DENOM,
  maxBlockSize = MAX_BLOCK_SIZE,
  randomSeed = RANDOM_SEED,
  deepSplit = DEEP_SPLIT,
  detectCutHeight = DETECT_CUT_HEIGHT,
  minModuleSize = MIN_MODULE_SIZE,
  pamStage = PAM_STAGE,
  pamRespectsDendro = PAM_RESPECTS_DENDRO,
  reassignThreshold = REASSIGN_THRESHOLD,
  minKMEtoStay = MIN_KME_TO_STAY,
  mergeCutHeight = MERGE_CUT_HEIGHT,
  numericLabels = NUMERIC_LABELS,
  saveTOMs = SAVE_TOMS,
  primary_mta_min_evaluable_genes = PRIMARY_MTA_MIN_EVALUABLE_GENES,
  primary_mta_min_target_coverage = PRIMARY_MTA_MIN_TARGET_COVERAGE,
  stringsAsFactors = FALSE
)
write.table(
  params,
  params_out,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

cat(strrep("=", 120), "\n", sep = "")
cat("02d TCGA WGCNA MODULE-CONSTRUCTION CONTRACT: PASS\n")
cat(strrep("=", 120), "\n", sep = "")
cat("Frozen soft power:          ", frozen_power, "\n", sep = "")
cat("Frozen network type:        ", NETWORK_TYPE, "\n", sep = "")
cat("Frozen TOM type:            ", TOM_TYPE, "\n", sep = "")
cat("Frozen maxBlockSize:        ", MAX_BLOCK_SIZE, " (single block)\n", sep = "")
cat("Frozen minModuleSize:       ", MIN_MODULE_SIZE, "\n", sep = "")
cat("Frozen deepSplit:           ", DEEP_SPLIT, "\n", sep = "")
cat("Frozen mergeCutHeight:      ", MERGE_CUT_HEIGHT, "\n", sep = "")
cat("Frozen PAM respects dendro: ", PAM_RESPECTS_DENDRO, "\n", sep = "")
cat("Frozen kME reassignment:    OFF (reassignThreshold=0)\n")
cat("Frozen minKMEtoStay:        ", MIN_KME_TO_STAY, "\n", sep = "")
cat("\nPrimary MTA assessability:\n")
cat("  minimum evaluable genes:  ", PRIMARY_MTA_MIN_EVALUABLE_GENES, "\n", sep = "")
cat("  minimum target coverage:  ", 100 * PRIMARY_MTA_MIN_TARGET_COVERAGE, "%\n", sep = "")
cat("\nNo modules were calculated by this script.\n")
cat("No target dataset was opened.\n")
cat("\nOutputs:\n")
cat("  ", params_out, "\n", sep = "")
cat("  ", json_out, "\n", sep = "")
cat("  ", md_out, "\n", sep = "")
cat(strrep("=", 120), "\n", sep = "")
