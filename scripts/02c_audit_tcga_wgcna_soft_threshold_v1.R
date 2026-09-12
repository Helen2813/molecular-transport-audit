options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "02c-audit-tcga-wgcna-soft-threshold-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"
TCGA_EXPR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_input_audit_v1",
  "staged_continuous_inputs",
  "TCGA_BRCA_PanCanAtlas2018",
  "data_mrna_seq_v2_rsem.txt"
)
SOURCE_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_tcga_source_universe_v1")
SAMPLE_MANIFEST <- file.path(SOURCE_ROOT, "tcga_source_samples_frozen_v1.tsv")
GENE_MANIFEST <- file.path(SOURCE_ROOT, "tcga_source_gene_universe_frozen_v1.tsv")

OUT_DIR <- file.path(DATA_ROOT, "paper4_tcbb_wgcna_soft_threshold_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# PREDECLARED SOURCE-ONLY WGCNA CONTRACT
# -----------------------------------------------------------------------------
NETWORK_TYPE <- "signed"
COR_FNC <- "cor"
COR_USE <- "p"

POWERS <- c(1:12, seq(14, 30, by = 2))

SFT_R2_TARGET <- 0.80
MEAN_CONNECTIVITY_FLOOR <- 20.0
REQUIRE_NEGATIVE_SLOPE <- TRUE

# Deterministic power-selection rule:
# 1) smallest power with SFT.R.sq >= 0.80, slope < 0, mean.k >= 20.
# 2) if none: among slope < 0 and mean.k >= 20, choose highest SFT.R.sq;
#    ties -> lower power.
# 3) if none: among slope < 0, choose highest SFT.R.sq; ties -> lower power.
# 4) final fallback: highest SFT.R.sq; ties -> lower power.
#
# No target dataset is used anywhere in this script.

cat(paste0(strrep("=", 120), "\n"))
cat("Paper 4 / TCBB - TCGA source-only WGCNA soft-threshold diagnostic\n")
cat(paste0(strrep("=", 120), "\n"))
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("\nScientific guard:\n")
cat("  SCAN-B opened:                                      NO\n")
cat("  METABRIC opened:                                    NO\n")
cat("  GSE239948 opened:                                   NO\n")
cat("  Clinical outcomes loaded:                           NO\n")
cat("  Target coverage/preservation used:                  NO\n")
cat("  Module detection performed:                         NO\n")
cat("  Operation: frozen-source WGCNA power diagnostic only\n")
cat(paste0(strrep("=", 120), "\n\n"))

required <- c(TCGA_EXPR, SAMPLE_MANIFEST, GENE_MANIFEST)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required files:\n", paste(missing, collapse = "\n"))
}

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop(
    "R package 'WGCNA' is not installed.\n",
    "Install it while broadband is available, then rerun:\n",
    "  if (!requireNamespace('BiocManager', quietly=TRUE)) ",
    "install.packages('BiocManager', repos='https://cloud.r-project.org')\n",
    "  BiocManager::install('WGCNA', ask=FALSE, update=FALSE)\n"
  )
}

suppressPackageStartupMessages(library(WGCNA))

# Keep the source analysis deterministic. WGCNA's pickSoftThreshold is
# deterministic for a fixed matrix/power vector, so no random seed is needed.
options(stringsAsFactors = FALSE)

cat("[1/5] Loading frozen source manifests ...\n")
sample_manifest <- read.delim(
  SAMPLE_MANIFEST,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)
gene_manifest <- read.delim(
  GENE_MANIFEST,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

selected_samples <- sample_manifest$sample_id[
  as.integer(sample_manifest$include_in_source_network) == 1
]

if (length(selected_samples) != 1082) {
  stop("Expected 1,082 frozen source samples; found ", length(selected_samples))
}
if (nrow(gene_manifest) != 10000) {
  stop("Expected 10,000 frozen source genes; found ", nrow(gene_manifest))
}
if (!"source_row_index_0based" %in% names(gene_manifest)) {
  stop("Gene manifest lacks source_row_index_0based.")
}

cat("  frozen samples: ", length(selected_samples), "\n", sep = "")
cat("  frozen genes:   ", nrow(gene_manifest), "\n", sep = "")

cat("\n[2/5] Loading TCGA continuous expression and reconstructing frozen source matrix ...\n")
expr <- read.delim(
  TCGA_EXPR,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

if (nrow(expr) != 20531) {
  stop("Expected 20,531 TCGA gene rows; found ", nrow(expr))
}
if (!all(c("Hugo_Symbol", "Entrez_Gene_Id") %in% names(expr))) {
  stop("TCGA expression identifier columns are not as expected.")
}
missing_samples <- setdiff(selected_samples, names(expr))
if (length(missing_samples) > 0) {
  stop("Frozen source samples missing from expression file: ",
       paste(head(missing_samples, 20), collapse = ", "))
}

row_idx <- as.integer(gene_manifest$source_row_index_0based) + 1L
if (any(row_idx < 1L | row_idx > nrow(expr))) {
  stop("Frozen source gene row index is out of bounds.")
}

expr_sub <- expr[row_idx, c("Hugo_Symbol", "Entrez_Gene_Id", selected_samples), drop = FALSE]

manifest_symbols <- as.character(gene_manifest$Hugo_Symbol)
expr_symbols <- as.character(expr_sub$Hugo_Symbol)
if (!identical(manifest_symbols, expr_symbols)) {
  bad <- which(manifest_symbols != expr_symbols)
  stop(
    "Frozen gene manifest does not replay exact TCGA source rows. First mismatch at rank ",
    bad[1], ": manifest=", manifest_symbols[bad[1]],
    " expression=", expr_symbols[bad[1]]
  )
}

x_raw <- as.matrix(expr_sub[, selected_samples, drop = FALSE])
storage.mode(x_raw) <- "double"

if (any(!is.finite(x_raw))) {
  stop("Non-finite value present in frozen TCGA source matrix.")
}
if (any(x_raw < 0)) {
  stop("Negative RSEM value present; frozen expression contract violated.")
}

x_log <- log2(x_raw + 1.0)

# WGCNA requires samples x genes.
datExpr <- t(x_log)
colnames(datExpr) <- manifest_symbols
rownames(datExpr) <- selected_samples

rm(expr, expr_sub, x_raw, x_log)
gc(verbose = FALSE)

cat("  datExpr dimensions: ", nrow(datExpr), " samples x ",
    ncol(datExpr), " genes\n", sep = "")

cat("\n[3/5] WGCNA goodSamplesGenes replay check ...\n")
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
if (!isTRUE(gsg$allOK)) {
  bad_samples <- sum(!gsg$goodSamples)
  bad_genes <- sum(!gsg$goodGenes)
  stop(
    "Frozen source matrix failed goodSamplesGenes: bad samples=",
    bad_samples, ", bad genes=", bad_genes,
    ". Do not silently remove them; review the source-universe contract."
  )
}
cat("  all frozen samples/genes pass goodSamplesGenes: YES\n")

# Save the exact source matrix used for WGCNA so all downstream source analyses
# replay the same object without reinterpreting the raw file.
datExpr_path <- file.path(OUT_DIR, "tcga_source_datExpr_1082x10000_v1.rds")
saveRDS(datExpr, datExpr_path, compress = "xz")
cat("  frozen source matrix RDS: ", datExpr_path, "\n", sep = "")

cat("\n[4/5] Running pickSoftThreshold with predeclared signed/Pearson contract ...\n")
cat("  networkType: signed\n")
cat("  correlation: Pearson\n")
cat("  powers:      ", paste(POWERS, collapse = ", "), "\n", sep = "")
cat("  primary selection threshold: SFT.R.sq >= ", SFT_R2_TARGET,
    ", slope < 0, mean.k >= ", MEAN_CONNECTIVITY_FLOOR, "\n", sep = "")

sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = POWERS,
  networkType = NETWORK_TYPE,
  corFnc = COR_FNC,
  corOptions = list(use = COR_USE),
  verbose = 5,
  blockSize = 2000
)

fit <- as.data.frame(sft$fitIndices, check.names = FALSE)

required_fit_cols <- c("Power", "SFT.R.sq", "slope", "mean.k.")
missing_fit_cols <- setdiff(required_fit_cols, names(fit))
if (length(missing_fit_cols) > 0) {
  stop(
    "Unexpected WGCNA fitIndices columns; missing: ",
    paste(missing_fit_cols, collapse = ", "),
    "\nAvailable: ", paste(names(fit), collapse = ", ")
  )
}

fit$qualifies_primary_rule <- (
  is.finite(fit$SFT.R.sq) &
  fit$SFT.R.sq >= SFT_R2_TARGET &
  is.finite(fit$slope) &
  fit$slope < 0 &
  is.finite(fit$mean.k.) &
  fit$mean.k. >= MEAN_CONNECTIVITY_FLOOR
)

choose_lowest_power <- function(df) {
  df <- df[order(df$Power), , drop = FALSE]
  df[1, , drop = FALSE]
}

choose_best_r2 <- function(df) {
  df <- df[order(-df$SFT.R.sq, df$Power), , drop = FALSE]
  df[1, , drop = FALSE]
}

qualified <- fit[fit$qualifies_primary_rule, , drop = FALSE]

if (nrow(qualified) > 0) {
  chosen <- choose_lowest_power(qualified)
  selection_reason <- paste0(
    "primary rule: lowest power with SFT.R.sq >= ",
    SFT_R2_TARGET, ", negative slope, mean.k >= ",
    MEAN_CONNECTIVITY_FLOOR
  )
  fallback_level <- 0L
} else {
  pool1 <- fit[
    is.finite(fit$slope) & fit$slope < 0 &
      is.finite(fit$mean.k.) & fit$mean.k. >= MEAN_CONNECTIVITY_FLOOR &
      is.finite(fit$SFT.R.sq),
    , drop = FALSE
  ]

  if (nrow(pool1) > 0) {
    chosen <- choose_best_r2(pool1)
    selection_reason <- paste0(
      "fallback 1: no power met R2 target; highest SFT.R.sq among ",
      "negative-slope powers with mean.k >= ", MEAN_CONNECTIVITY_FLOOR
    )
    fallback_level <- 1L
  } else {
    pool2 <- fit[
      is.finite(fit$slope) & fit$slope < 0 & is.finite(fit$SFT.R.sq),
      , drop = FALSE
    ]

    if (nrow(pool2) > 0) {
      chosen <- choose_best_r2(pool2)
      selection_reason <- (
        "fallback 2: highest SFT.R.sq among negative-slope powers"
      )
      fallback_level <- 2L
    } else {
      pool3 <- fit[is.finite(fit$SFT.R.sq), , drop = FALSE]
      if (nrow(pool3) == 0) {
        stop("No finite WGCNA SFT.R.sq values were produced.")
      }
      chosen <- choose_best_r2(pool3)
      selection_reason <- "fallback 3: highest finite SFT.R.sq"
      fallback_level <- 3L
    }
  }
}

chosen_power <- as.numeric(chosen$Power[1])

fit_out <- file.path(OUT_DIR, "tcga_wgcna_soft_threshold_diagnostics_v1.tsv")
write.table(
  fit,
  fit_out,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

selection_out <- file.path(OUT_DIR, "tcga_wgcna_soft_power_frozen_v1.tsv")
selection_table <- data.frame(
  script_version = SCRIPT_VERSION,
  network_type = NETWORK_TYPE,
  correlation = "Pearson",
  SFT_R2_target = SFT_R2_TARGET,
  mean_connectivity_floor = MEAN_CONNECTIVITY_FLOOR,
  chosen_power = chosen_power,
  chosen_SFT_R_sq = as.numeric(chosen$SFT.R.sq[1]),
  chosen_slope = as.numeric(chosen$slope[1]),
  chosen_mean_k = as.numeric(chosen$mean.k.[1]),
  fallback_level = fallback_level,
  selection_reason = selection_reason,
  stringsAsFactors = FALSE
)
write.table(
  selection_table,
  selection_out,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Diagnostic figure is source-only and does not alter the frozen selection rule.
png(
  filename = file.path(OUT_DIR, "tcga_wgcna_soft_threshold_diagnostic_v1.png"),
  width = 1800,
  height = 800,
  res = 150
)
par(mfrow = c(1, 2))

plot(
  fit$Power,
  fit$SFT.R.sq,
  type = "b",
  xlab = "Soft-threshold power",
  ylab = "Scale-free topology fit (R^2)",
  main = "Source-only scale-free topology fit"
)
abline(h = SFT_R2_TARGET, lty = 2)
abline(v = chosen_power, lty = 3)

plot(
  fit$Power,
  fit$mean.k.,
  type = "b",
  xlab = "Soft-threshold power",
  ylab = "Mean connectivity",
  main = "Source-only mean connectivity"
)
abline(h = MEAN_CONNECTIVITY_FLOOR, lty = 2)
abline(v = chosen_power, lty = 3)

dev.off()

cat("\n[5/5] Freezing source-only soft power ...\n")
cat("  chosen power:     ", chosen_power, "\n", sep = "")
cat("  SFT.R.sq:         ", signif(chosen$SFT.R.sq[1], 6), "\n", sep = "")
cat("  slope:            ", signif(chosen$slope[1], 6), "\n", sep = "")
cat("  mean connectivity:", signif(chosen$mean.k.[1], 6), "\n", sep = "")
cat("  fallback level:   ", fallback_level, "\n", sep = "")
cat("  reason:           ", selection_reason, "\n", sep = "")

contract_md <- file.path(OUT_DIR, "tcga_wgcna_soft_threshold_contract_v1.md")
writeLines(
  c(
    "# Paper 4 / TCBB — TCGA-BRCA WGCNA Soft-Threshold Contract v1",
    "",
    "Status: **FROZEN**",
    "",
    "No SCAN-B, METABRIC, GSE239948, clinical outcome, target coverage, or",
    "target preservation result was used.",
    "",
    "## Frozen source-network choices",
    "",
    "- Network type: **signed**",
    "- Source correlation: **Pearson**",
    paste0("- Candidate powers: `", paste(POWERS, collapse = ", "), "`"),
    paste0("- Scale-free fit target: **R^2 >= ", SFT_R2_TARGET, "**"),
    paste0("- Required slope under primary rule: **negative**"),
    paste0("- Mean-connectivity floor under primary rule: **",
           MEAN_CONNECTIVITY_FLOOR, "**"),
    paste0("- Frozen soft power: **", chosen_power, "**"),
    paste0("- Selection path: **fallback level ", fallback_level, "**"),
    "",
    "The power was selected by the predeclared deterministic rule in",
    "`02c_audit_tcga_wgcna_soft_threshold_v1.R`; no module-preservation",
    "result was available or inspected."
  ),
  con = contract_md
)

cat("\n", strrep("=", 120), "\n", sep = "")
cat("02c TCGA SOURCE-ONLY WGCNA SOFT-THRESHOLD: PASS\n")
cat(strrep("=", 120), "\n", sep = "")
cat("Frozen network type: signed\n")
cat("Frozen source correlation: Pearson\n")
cat("Frozen soft power: ", chosen_power, "\n", sep = "")
cat("Fallback level:    ", fallback_level, "\n", sep = "")
cat("\nNo target dataset was opened.\n")
cat("No modules or preservation statistics were calculated.\n")
cat("\nOutputs:\n")
cat("  ", fit_out, "\n", sep = "")
cat("  ", selection_out, "\n", sep = "")
cat("  ", datExpr_path, "\n", sep = "")
cat("  ", contract_md, "\n", sep = "")
cat("  ", file.path(OUT_DIR, "tcga_wgcna_soft_threshold_diagnostic_v1.png"), "\n", sep = "")
cat(strrep("=", 120), "\n", sep = "")
