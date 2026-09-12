options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "02e-build-tcga-source-modules-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"
SFT_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_wgcna_soft_threshold_v1")
CONTRACT_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_wgcna_module_contract_v2")

SOURCE_MATRIX <- file.path(SFT_ROOT, "tcga_source_datExpr_1082x10000_v1.rds")
CONTRACT_JSON <- file.path(CONTRACT_ROOT, "tcga_wgcna_module_contract_v2.json")

OUT_DIR <- file.path(DATA_ROOT, "paper4_tcbb_tcga_source_modules_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(strrep("=", 122), "\n", sep = "")
cat("Paper 4 / TCBB - build frozen TCGA-BRCA source modules\n")
cat(strrep("=", 122), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("\nScientific guard:\n")
cat("  SCAN-B opened:                                      NO\n")
cat("  METABRIC opened:                                    NO\n")
cat("  GSE239948 opened:                                   NO\n")
cat("  Clinical outcomes loaded:                           NO\n")
cat("  Target coverage/preservation used:                  NO\n")
cat("  Module parameters tuned after viewing modules:       NO\n")
cat("  All non-grey modules will be retained automatically: YES\n")
cat(strrep("=", 122), "\n\n", sep = "")

for (f in c(SOURCE_MATRIX, CONTRACT_JSON)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}
if (!requireNamespace("WGCNA", quietly = TRUE)) stop("WGCNA is not installed.")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed.")

suppressPackageStartupMessages(library(WGCNA))

contract <- jsonlite::read_json(CONTRACT_JSON, simplifyVector = TRUE)
P <- contract$blockwiseModules

# Hard replay checks.
if (as.numeric(P$power) != 8) stop("Contract replay failure: power.")
if (P$networkType != "signed") stop("Contract replay failure: networkType.")
if (P$corType != "pearson") stop("Contract replay failure: corType.")
if (P$TOMType != "signed") stop("Contract replay failure: TOMType.")
if (as.integer(P$maxBlockSize) < 10000) stop("Contract does not force a single block.")

datExpr <- readRDS(SOURCE_MATRIX)
if (nrow(datExpr) != 1082 || ncol(datExpr) != 10000) {
  stop("Unexpected source matrix dimensions.")
}
if (any(!is.finite(as.matrix(datExpr)))) stop("Non-finite source expression.")

# One-block invariant.
if (ncol(datExpr) > as.integer(P$maxBlockSize)) {
  stop("Frozen maxBlockSize would split the network; contract violated.")
}

tom_base <- file.path(OUT_DIR, "tcga_source_TOM_v1")

cat("[1/4] Running frozen WGCNA::blockwiseModules ...\n")
cat("  source matrix: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes\n", sep = "")
cat("  signed / Pearson / power 8 / signed TOM\n")
cat("  deepSplit=2, minModuleSize=30, mergeCutHeight=0.25\n")
cat("  This is the first module-detection run under the frozen v2 contract.\n\n")

net <- WGCNA::blockwiseModules(
  datExpr = datExpr,
  checkMissingData = isTRUE(P$checkMissingData),

  maxBlockSize = as.integer(P$maxBlockSize),
  blockSizePenaltyPower = as.numeric(P$blockSizePenaltyPower),
  randomSeed = as.integer(P$randomSeed),
  loadTOM = isTRUE(P$loadTOM),

  corType = P$corType,
  maxPOutliers = as.numeric(P$maxPOutliers),
  quickCor = as.numeric(P$quickCor),
  pearsonFallback = P$pearsonFallback,
  cosineCorrelation = isTRUE(P$cosineCorrelation),

  power = as.numeric(P$power),
  networkType = P$networkType,
  replaceMissingAdjacencies = isTRUE(P$replaceMissingAdjacencies),

  TOMType = P$TOMType,
  TOMDenom = P$TOMDenom,
  suppressTOMForZeroAdjacencies = isTRUE(P$suppressTOMForZeroAdjacencies),
  suppressNegativeTOM = isTRUE(P$suppressNegativeTOM),

  saveTOMs = isTRUE(P$saveTOMs),
  saveTOMFileBase = tom_base,

  deepSplit = as.integer(P$deepSplit),
  detectCutHeight = as.numeric(P$detectCutHeight),
  minModuleSize = as.integer(P$minModuleSize),

  useBranchEigennodeDissim = isTRUE(P$useBranchEigennodeDissim),
  pamStage = isTRUE(P$pamStage),
  pamRespectsDendro = isTRUE(P$pamRespectsDendro),

  reassignThreshold = as.numeric(P$reassignThreshold),
  minCoreKME = as.numeric(P$minCoreKME),
  minCoreKMESize = as.integer(P$minCoreKMESize),
  minKMEtoStay = as.numeric(P$minKMEtoStay),

  mergeCutHeight = as.numeric(P$mergeCutHeight),
  impute = isTRUE(P$impute),
  trapErrors = isTRUE(P$trapErrors),

  numericLabels = isTRUE(P$numericLabels),
  nThreads = as.integer(P$nThreads),
  useInternalMatrixAlgebra = isTRUE(P$useInternalMatrixAlgebra),
  useCorOptionsThroughout = isTRUE(P$useCorOptionsThroughout),

  verbose = 5
)

cat("\n[2/4] Validating frozen module result ...\n")

if (length(net$colors) != ncol(datExpr)) stop("Module-label length mismatch.")
if (length(net$blocks) != ncol(datExpr)) stop("Block-label length mismatch.")
if (length(unique(net$blocks)) != 1) {
  stop("Expected exactly one WGCNA block; found ", length(unique(net$blocks)))
}

labels <- as.integer(net$colors)
colors <- WGCNA::labels2colors(labels)

tab <- sort(table(labels), decreasing = TRUE)
non_grey_labels <- sort(unique(labels[labels != 0]))
n_modules <- length(non_grey_labels)
grey_n <- sum(labels == 0)
grey_fraction <- grey_n / length(labels)

if (n_modules < 1) stop("No non-grey source modules detected.")

cat("  non-grey modules: ", n_modules, "\n", sep = "")
cat("  grey genes:       ", grey_n, " / ", length(labels),
    " (", sprintf("%.2f", 100 * grey_fraction), "%)\n", sep = "")

# Canonical module IDs are label-based and automatic; no module is selected by result.
canonical_id <- rep(NA_character_, length(labels))
for (lab in non_grey_labels) {
  canonical_id[labels == lab] <- sprintf("TCGA_M%03d", lab)
}
canonical_id[labels == 0] <- "GREY"

membership <- data.frame(
  source_gene_order = seq_len(ncol(datExpr)),
  Hugo_Symbol = colnames(datExpr),
  wgcna_label = labels,
  wgcna_color = colors,
  program_id = canonical_id,
  stringsAsFactors = FALSE
)

sizes <- as.data.frame(table(
  wgcna_label = membership$wgcna_label,
  program_id = membership$program_id,
  wgcna_color = membership$wgcna_color
), stringsAsFactors = FALSE)
sizes <- sizes[sizes$Freq > 0, ]
names(sizes)[names(sizes) == "Freq"] <- "n_genes"
sizes$wgcna_label <- as.integer(as.character(sizes$wgcna_label))
sizes <- sizes[order(sizes$wgcna_label), ]

# Eigengene structure: SOURCE ONLY; no target information.
MEs <- as.data.frame(net$MEs, check.names = FALSE)
if (ncol(MEs) < 1) stop("No module eigengenes returned.")

me_cor <- cor(MEs, use = "pairwise.complete.obs", method = "pearson")
diag(me_cor) <- NA_real_
max_positive_me_cor <- if (ncol(MEs) > 1) max(me_cor, na.rm = TRUE) else NA_real_
min_me_cor <- if (ncol(MEs) > 1) min(me_cor, na.rm = TRUE) else NA_real_

cat("  source eigengenes returned: ", ncol(MEs), "\n", sep = "")
cat("  max off-diagonal source ME correlation: ",
    ifelse(is.na(max_positive_me_cor), "NA", signif(max_positive_me_cor, 6)), "\n", sep = "")

cat("\n[3/4] Writing immutable source-module artifacts ...\n")
membership_out <- file.path(OUT_DIR, "tcga_source_module_membership_frozen_v1.tsv")
sizes_out <- file.path(OUT_DIR, "tcga_source_module_summary_v1.tsv")
me_out <- file.path(OUT_DIR, "tcga_source_module_eigengenes_v1.tsv")
me_cor_out <- file.path(OUT_DIR, "tcga_source_module_eigengene_correlations_v1.tsv")
net_out <- file.path(OUT_DIR, "tcga_source_wgcna_network_v1.rds")

write.table(membership, membership_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sizes, sizes_out, sep = "\t", row.names = FALSE, quote = FALSE)

ME_table <- cbind(sample_id = rownames(MEs), MEs)
write.table(ME_table, me_out, sep = "\t", row.names = FALSE, quote = FALSE)

me_cor_df <- cbind(module_eigengene = rownames(me_cor), as.data.frame(me_cor, check.names = FALSE))
write.table(me_cor_df, me_cor_out, sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(net, net_out, compress = "gzip")

# Source-only diagnostic dendrogram; this cannot alter the frozen modules.
png(
  file.path(OUT_DIR, "tcga_source_module_dendrogram_v1.png"),
  width = 2400, height = 1100, res = 150
)
WGCNA::plotDendroAndColors(
  net$dendrograms[[1]],
  colors[net$blockGenes[[1]]],
  "Frozen source modules",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

cat("\n[4/4] Freezing ALL non-grey source modules ...\n")

result <- list(
  result_id = "paper4-tcbb-tcga-source-modules-v1",
  script_version = SCRIPT_VERSION,
  status = "FROZEN_SOURCE_MODULES",
  scientific_guard = list(
    scanb_opened = FALSE,
    metabric_opened = FALSE,
    gse239948_opened = FALSE,
    clinical_outcomes_loaded = FALSE,
    target_coverage_or_preservation_used = FALSE,
    posthoc_module_parameter_tuning = FALSE
  ),
  source = list(samples = nrow(datExpr), genes = ncol(datExpr)),
  contract = CONTRACT_JSON,
  wgcna = list(
    version = as.character(packageVersion("WGCNA")),
    blocks = length(unique(net$blocks)),
    non_grey_modules = n_modules,
    grey_genes = grey_n,
    grey_fraction = grey_fraction,
    module_sizes = setNames(
      as.list(as.integer(sizes$n_genes[sizes$wgcna_label != 0])),
      sizes$program_id[sizes$wgcna_label != 0]
    )
  ),
  source_eigengene_structure = list(
    n_eigengenes = ncol(MEs),
    max_off_diagonal_pearson_correlation = max_positive_me_cor,
    min_off_diagonal_pearson_correlation = min_me_cor
  ),
  freeze_rule = (
    "Every non-grey module returned by the pre-frozen WGCNA contract is retained. "
    "No module is selected or discarded based on size beyond the pre-frozen detection rule, "
    "biology, target coverage, or target preservation."
  )
)

json_out <- file.path(OUT_DIR, "tcga_source_modules_frozen_v1.json")
jsonlite::write_json(result, json_out, pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat(strrep("=", 122), "\n", sep = "")
cat("02e FROZEN TCGA SOURCE MODULES: PASS\n")
cat(strrep("=", 122), "\n", sep = "")
cat("Source samples:       ", nrow(datExpr), "\n", sep = "")
cat("Source genes:         ", ncol(datExpr), "\n", sep = "")
cat("WGCNA blocks:         ", length(unique(net$blocks)), "\n", sep = "")
cat("Non-grey modules:     ", n_modules, "\n", sep = "")
cat("Grey genes:           ", grey_n, " (", sprintf("%.2f", 100 * grey_fraction), "%)\n", sep = "")
cat("\nModule sizes (all frozen, no selection):\n")
for (i in seq_len(nrow(sizes))) {
  if (sizes$wgcna_label[i] != 0) {
    cat("  ", sprintf("%-10s", sizes$program_id[i]),
        " label=", sizes$wgcna_label[i],
        " color=", sizes$wgcna_color[i],
        " n=", sizes$n_genes[i], "\n", sep = "")
  }
}
cat("\nAll non-grey modules are now frozen.\n")
cat("No target dataset was opened.\n")
cat("No target preservation statistic was calculated.\n")
cat("\nOutputs:\n")
cat("  ", membership_out, "\n", sep = "")
cat("  ", sizes_out, "\n", sep = "")
cat("  ", me_out, "\n", sep = "")
cat("  ", me_cor_out, "\n", sep = "")
cat("  ", net_out, "\n", sep = "")
cat("  ", json_out, "\n", sep = "")
cat("  TOM base: ", tom_base, "\n", sep = "")
cat(strrep("=", 122), "\n", sep = "")
