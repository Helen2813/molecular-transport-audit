options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "02d-freeze-tcga-wgcna-module-contract-v2-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"
SFT_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_wgcna_soft_threshold_v1")
POWER_FILE <- file.path(SFT_ROOT, "tcga_wgcna_soft_power_frozen_v1.tsv")
SOURCE_MATRIX <- file.path(SFT_ROOT, "tcga_source_datExpr_1082x10000_v1.rds")

OUT_DIR <- file.path(DATA_ROOT, "paper4_tcbb_wgcna_module_contract_v2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPECTED_POWER <- 8

# Complete frozen WGCNA call contract.
P <- list(
  power = 8,
  networkType = "signed",
  corType = "pearson",
  maxPOutliers = 1,
  quickCor = 0,
  pearsonFallback = "individual",
  cosineCorrelation = FALSE,
  replaceMissingAdjacencies = FALSE,

  TOMType = "signed",
  TOMDenom = "min",
  suppressTOMForZeroAdjacencies = FALSE,
  suppressNegativeTOM = FALSE,

  maxBlockSize = 12000L,
  blockSizePenaltyPower = 5,
  randomSeed = 54321L,
  loadTOM = FALSE,

  deepSplit = 2L,
  detectCutHeight = 0.995,
  minModuleSize = 30L,

  useBranchEigennodeDissim = FALSE,
  pamStage = TRUE,
  pamRespectsDendro = FALSE,

  reassignThreshold = 0,
  minCoreKME = 0.5,
  minCoreKMESize = 10L,
  minKMEtoStay = 0,

  mergeCutHeight = 0.25,
  impute = TRUE,
  trapErrors = FALSE,

  numericLabels = TRUE,
  nThreads = 0L,
  useInternalMatrixAlgebra = FALSE,
  useCorOptionsThroughout = TRUE,
  checkMissingData = TRUE,
  saveTOMs = TRUE
)

PRIMARY_MTA_MIN_EVALUABLE_GENES <- 30L
PRIMARY_MTA_MIN_TARGET_COVERAGE <- 0.80

cat(strrep("=", 122), "\n", sep = "")
cat("Paper 4 / TCBB - freeze COMPLETE TCGA-BRCA WGCNA module-construction contract v2\n")
cat(strrep("=", 122), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("\nReason for v2:\n")
cat("  v1 froze the main module parameters before module detection.\n")
cat("  v2 explicitly freezes remaining scientifically relevant blockwiseModules defaults,\n")
cat("  including minCoreKME/minCoreKMESize and missing/TOM behavior.\n")
cat("  NO module detection occurred between v1 and v2.\n")
cat("\nScientific guard:\n")
cat("  SCAN-B opened:                                      NO\n")
cat("  METABRIC opened:                                    NO\n")
cat("  GSE239948 opened:                                   NO\n")
cat("  Clinical outcomes loaded:                           NO\n")
cat("  Target coverage/preservation inspected:             NO\n")
cat("  WGCNA modules calculated by this script:            NO\n")
cat(strrep("=", 122), "\n\n", sep = "")

for (f in c(POWER_FILE, SOURCE_MATRIX)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

if (!requireNamespace("WGCNA", quietly = TRUE)) stop("WGCNA is not installed.")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed.")

power_tbl <- read.delim(
  POWER_FILE, sep = "\t", header = TRUE, check.names = FALSE,
  quote = "", comment.char = ""
)
if (nrow(power_tbl) != 1) stop("Expected exactly one frozen power record.")
frozen_power <- as.numeric(power_tbl$chosen_power[1])
if (!isTRUE(all.equal(frozen_power, EXPECTED_POWER))) {
  stop("Expected frozen power 8; found ", frozen_power)
}

datExpr <- readRDS(SOURCE_MATRIX)
if (nrow(datExpr) != 1082 || ncol(datExpr) != 10000) {
  stop("Unexpected source matrix dimensions: ", nrow(datExpr), " x ", ncol(datExpr))
}
if (any(!is.finite(as.matrix(datExpr)))) stop("Source matrix contains non-finite values.")
rm(datExpr)
gc(verbose = FALSE)

contract <- list(
  contract_id = "paper4-tcbb-wgcna-module-contract-v2",
  supersedes_for_execution = "paper4-tcbb-wgcna-module-contract-v1",
  script_version = SCRIPT_VERSION,
  status = "FROZEN_BEFORE_FIRST_MODULE_DETECTION",
  amendment_reason = paste(
    "Explicitly records remaining blockwiseModules defaults that can affect",
    "module disbanding/trimming or network construction. No modules had been",
    "calculated before this amendment."
  ),
  scientific_guard = list(
    scanb_opened = FALSE,
    metabric_opened = FALSE,
    gse239948_opened = FALSE,
    clinical_outcomes_loaded = FALSE,
    target_coverage_or_preservation_inspected = FALSE,
    modules_calculated_before_v2_freeze = FALSE
  ),
  source_input = list(
    matrix = SOURCE_MATRIX,
    samples = 1082L,
    genes = 10000L
  ),
  blockwiseModules = P,
  assessability = list(
    minimum_evaluable_genes = PRIMARY_MTA_MIN_EVALUABLE_GENES,
    minimum_target_gene_coverage = PRIMARY_MTA_MIN_TARGET_COVERAGE,
    rule = paste0(
      "Primary-assessable only if evaluable genes >= ",
      PRIMARY_MTA_MIN_EVALUABLE_GENES,
      " AND frozen-gene target coverage >= ",
      100 * PRIMARY_MTA_MIN_TARGET_COVERAGE,
      "%. Otherwise label not assessable."
    )
  ),
  software = list(
    R = R.version.string,
    WGCNA = as.character(packageVersion("WGCNA"))
  )
)

json_out <- file.path(OUT_DIR, "tcga_wgcna_module_contract_v2.json")
jsonlite::write_json(contract, json_out, pretty = TRUE, auto_unbox = TRUE, digits = NA)

# Flatten parameters for exact replay by 02e.
param_df <- data.frame(
  name = names(P),
  value = vapply(P, function(x) paste(x, collapse = ","), character(1)),
  stringsAsFactors = FALSE
)
params_out <- file.path(OUT_DIR, "tcga_wgcna_module_parameters_frozen_v2.tsv")
write.table(param_df, params_out, sep = "\t", row.names = FALSE, quote = FALSE)

md_out <- file.path(OUT_DIR, "tcga_wgcna_module_contract_v2.md")
writeLines(
  c(
    "# Paper 4 / TCBB — Complete TCGA-BRCA WGCNA Module Contract v2",
    "",
    "Status: **FROZEN BEFORE FIRST MODULE DETECTION**",
    "",
    "This v2 contract supersedes v1 for execution only because v1 did not",
    "explicitly record several WGCNA defaults. No module detection occurred",
    "between the v1 and v2 freezes.",
    "",
    "Key additions made explicit in v2:",
    "- `minCoreKME = 0.5`",
    "- `minCoreKMESize = 10`",
    "- `checkMissingData = TRUE`",
    "- `replaceMissingAdjacencies = FALSE`",
    "- `suppressTOMForZeroAdjacencies = FALSE`",
    "- `suppressNegativeTOM = FALSE`",
    "- `impute = TRUE`",
    "- `trapErrors = FALSE`",
    "- `useInternalMatrixAlgebra = FALSE`",
    "- `useCorOptionsThroughout = TRUE`",
    "",
    "All previously frozen scientific choices remain unchanged:",
    "- signed network, Pearson correlation",
    "- soft power 8",
    "- signed TOM with denominator `min`",
    "- one 10,000-gene block",
    "- deepSplit 2",
    "- detectCutHeight 0.995",
    "- minModuleSize 30",
    "- PAM stage on, `pamRespectsDendro = FALSE`",
    "- `reassignThreshold = 0`",
    "- `minKMEtoStay = 0`",
    "- mergeCutHeight 0.25",
    "- numeric module labels",
    "- TOM saving enabled",
    "",
    "Primary MTA assessability remains >=30 evaluable genes and >=80% coverage."
  ),
  md_out
)

cat(strrep("=", 122), "\n", sep = "")
cat("02d v2 COMPLETE WGCNA MODULE CONTRACT: PASS\n")
cat(strrep("=", 122), "\n", sep = "")
cat("Frozen soft power:        ", P$power, "\n", sep = "")
cat("Frozen minCoreKME:        ", P$minCoreKME, "\n", sep = "")
cat("Frozen minCoreKMESize:    ", P$minCoreKMESize, "\n", sep = "")
cat("Frozen minKMEtoStay:      ", P$minKMEtoStay, "\n", sep = "")
cat("Frozen mergeCutHeight:    ", P$mergeCutHeight, "\n", sep = "")
cat("Frozen maxBlockSize:      ", P$maxBlockSize, "\n", sep = "")
cat("\nNo modules were calculated.\n")
cat("No target dataset was opened.\n")
cat("\nOutputs:\n")
cat("  ", params_out, "\n", sep = "")
cat("  ", json_out, "\n", sep = "")
cat("  ", md_out, "\n", sep = "")
cat(strrep("=", 122), "\n", sep = "")
