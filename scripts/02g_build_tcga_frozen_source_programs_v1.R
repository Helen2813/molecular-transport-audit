options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "02g-build-tcga-frozen-source-programs-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"
MODULE_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_tcga_source_modules_v1")
CONTRACT_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_source_program_contract_v1")
SFT_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_wgcna_soft_threshold_v1")

SOURCE_MATRIX <- file.path(SFT_ROOT, "tcga_source_datExpr_1082x10000_v1.rds")
NETWORK_RDS <- file.path(MODULE_ROOT, "tcga_source_wgcna_network_v1.rds")
MEMBERSHIP_TSV <- file.path(MODULE_ROOT, "tcga_source_module_membership_frozen_v1.tsv")
CONTRACT_JSON <- file.path(CONTRACT_ROOT, "tcga_source_program_representation_contract_v1.json")

OUT_DIR <- file.path(DATA_ROOT, "paper4_tcbb_frozen_source_programs_v1")
EDGE_DIR <- file.path(OUT_DIR, "source_edge_vectors")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(EDGE_DIR, recursive = TRUE, showWarnings = FALSE)

ORIENTATION_TOL <- 1e-12

cat(strrep("=", 124), "\n", sep = "")
cat("Paper 4 / TCBB - build frozen TCGA source-program weights and source structure\n")
cat(strrep("=", 124), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("\nScientific guard:\n")
cat("  SCAN-B opened:                                      NO\n")
cat("  METABRIC opened:                                    NO\n")
cat("  GSE239948 opened:                                   NO\n")
cat("  Clinical outcomes loaded:                           NO\n")
cat("  Target preservation statistics calculated:          NO\n")
cat("  Source program representation contract frozen:      YES\n")
cat(strrep("=", 124), "\n\n", sep = "")

for (f in c(SOURCE_MATRIX, NETWORK_RDS, MEMBERSHIP_TSV, CONTRACT_JSON)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed.")
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is not installed.")

contract <- jsonlite::read_json(CONTRACT_JSON, simplifyVector = TRUE)
if (contract$status != "FROZEN_BEFORE_SOURCE_PROGRAM_STATISTICS") {
  stop("Source-program representation contract is not in expected frozen state.")
}

datExpr <- readRDS(SOURCE_MATRIX)
net <- readRDS(NETWORK_RDS)
membership <- read.delim(
  MEMBERSHIP_TSV, sep = "\t", header = TRUE, check.names = FALSE,
  quote = "", comment.char = ""
)

if (nrow(datExpr) != 1082 || ncol(datExpr) != 10000) stop("Unexpected source matrix dimensions.")
if (nrow(membership) != 10000) stop("Unexpected membership row count.")
if (!identical(as.character(membership$Hugo_Symbol), colnames(datExpr))) {
  stop("Membership gene order does not exactly match frozen source matrix.")
}

program_ids <- sort(unique(membership$program_id[membership$program_id != "GREY"]))
if (length(program_ids) != 12L) stop("Expected exactly 12 frozen non-grey programs.")

all_weights <- list()
all_summaries <- list()
score_table <- data.frame(sample_id = rownames(datExpr), stringsAsFactors = FALSE)

cat("[1/3] Computing frozen source representations for all 12 programs ...\n")

for (idx in seq_along(program_ids)) {
  pid <- program_ids[idx]
  rows <- which(membership$program_id == pid)
  genes <- membership$Hugo_Symbol[rows]
  label <- unique(membership$wgcna_label[rows])
  color <- unique(membership$wgcna_color[rows])

  if (length(label) != 1 || length(color) != 1) {
    stop("Program ", pid, " has inconsistent WGCNA label/color.")
  }

  X <- as.matrix(datExpr[, genes, drop = FALSE])
  storage.mode(X) <- "double"

  # Exact frozen within-source z-standardization.
  Xz <- scale(X, center = TRUE, scale = TRUE)
  if (any(!is.finite(Xz))) stop("Non-finite standardized value in ", pid)

  # First right singular vector = PC1 loading vector.
  pc <- svd(Xz, nu = 1, nv = 1)
  loading <- as.numeric(pc$v[, 1])
  names(loading) <- genes
  score <- as.numeric(pc$u[, 1] * pc$d[1])

  unweighted_score <- rowMeans(Xz)
  orient_cor <- suppressWarnings(cor(score, unweighted_score, method = "pearson"))

  flip <- FALSE
  orientation_path <- "correlation_with_unweighted_mean"

  if (is.finite(orient_cor) && orient_cor < -ORIENTATION_TOL) {
    flip <- TRUE
  } else if (!is.finite(orient_cor) || abs(orient_cor) <= ORIENTATION_TOL) {
    orientation_path <- "maximum_absolute_loading_fallback"
    max_abs <- max(abs(loading))
    candidates <- which(abs(abs(loading) - max_abs) <= .Machine$double.eps^0.5)
    candidate_genes <- genes[candidates]
    chosen_gene <- sort(candidate_genes)[1]
    chosen_idx <- which(genes == chosen_gene)[1]
    if (loading[chosen_idx] < 0) flip <- TRUE
  }

  if (flip) {
    loading <- -loading
    score <- -score
    orient_cor <- suppressWarnings(cor(score, unweighted_score, method = "pearson"))
  }

  # Existing MTA sign-score orientation object.
  loading_sign <- sign(loading)
  signed_score <- as.numeric(Xz %*% loading_sign / length(loading_sign))
  signed_score_cor_pc1 <- suppressWarnings(cor(score, signed_score, method = "pearson"))

  # Source edge structure in the frozen gene order.
  R <- crossprod(Xz) / (nrow(Xz) - 1)
  diag(R) <- 1
  edges <- R[upper.tri(R, diag = FALSE)]

  if (length(edges) != choose(length(genes), 2)) {
    stop("Edge-vector length mismatch in ", pid)
  }
  if (any(!is.finite(edges))) stop("Non-finite source edge in ", pid)

  abs_edges <- abs(edges)
  qs <- as.numeric(quantile(
    edges, probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
    names = FALSE, type = 7
  ))

  pc1_var_explained <- (pc$d[1]^2) / sum(Xz^2)

  # Cross-check against WGCNA merged-module eigengene.
  me_name <- paste0("ME", label)
  if (!me_name %in% colnames(net$MEs)) {
    stop("Expected WGCNA eigengene column ", me_name, " not found for ", pid)
  }
  wgcna_me <- as.numeric(net$MEs[, me_name])
  abs_cor_wgcna_me <- abs(cor(score, wgcna_me, method = "pearson"))

  weight_df <- data.frame(
    program_id = pid,
    wgcna_label = as.integer(label),
    wgcna_color = as.character(color),
    source_gene_order_within_program = seq_along(genes),
    Hugo_Symbol = genes,
    source_pc1_loading = loading,
    source_pc1_loading_sign = loading_sign,
    source_abs_loading_rank = rank(-abs(loading), ties.method = "first"),
    stringsAsFactors = FALSE
  )
  all_weights[[pid]] <- weight_df

  edge_file <- file.path(EDGE_DIR, paste0(pid, "_source_edges_v1.rds"))
  saveRDS(
    list(
      program_id = pid,
      genes_in_frozen_order = genes,
      correlation = "Pearson",
      edge_vector_order = "upper triangle in frozen source-gene order",
      source_edges = edges
    ),
    edge_file,
    compress = "gzip"
  )

  summary_row <- data.frame(
    program_id = pid,
    wgcna_label = as.integer(label),
    wgcna_color = as.character(color),
    n_genes = length(genes),
    n_edges = length(edges),
    orientation_path = orientation_path,
    orientation_flip_applied = as.integer(flip),
    source_pc1_vs_unweighted_mean_cor = orient_cor,
    source_signed_score_vs_pc1_cor = signed_score_cor_pc1,
    abs_cor_with_wgcna_module_eigengene = abs_cor_wgcna_me,
    source_pc1_variance_explained = pc1_var_explained,
    source_coherence_mean_abs_r = mean(abs_edges),
    source_coherence_median_abs_r = median(abs_edges),
    source_edge_mean_signed_r = mean(edges),
    source_edge_sd_signed_r = sd(edges),
    source_edge_q05 = qs[1],
    source_edge_q25 = qs[2],
    source_edge_median = qs[3],
    source_edge_q75 = qs[4],
    source_edge_q95 = qs[5],
    source_edge_fraction_positive = mean(edges > 0),
    source_edge_file = edge_file,
    stringsAsFactors = FALSE
  )
  all_summaries[[pid]] <- summary_row

  score_table[[paste0(pid, "_source_pc1")]] <- score
  score_table[[paste0(pid, "_source_signed_score")]] <- signed_score

  cat(
    "  ", sprintf("%-10s", pid),
    " n=", sprintf("%4d", length(genes)),
    " PC1=", sprintf("%.3f", pc1_var_explained),
    " mean|r|=", sprintf("%.3f", mean(abs_edges)),
    " signedScore~PC1=", sprintf("%.3f", signed_score_cor_pc1),
    " |PC1~WGCNA_ME|=", sprintf("%.6f", abs_cor_wgcna_me),
    "\n",
    sep = ""
  )

  rm(X, Xz, pc, R, edges, abs_edges, weight_df)
  gc(verbose = FALSE)
}

weights <- do.call(rbind, all_weights)
summary <- do.call(rbind, all_summaries)
rownames(weights) <- NULL
rownames(summary) <- NULL

cat("\n[2/3] Validating frozen source-program bundle ...\n")

if (nrow(summary) != 12L) stop("Expected 12 program summaries.")
if (nrow(weights) != sum(summary$n_genes)) stop("Weight-table size mismatch.")
if (any(summary$abs_cor_with_wgcna_module_eigengene < 0.999999)) {
  bad <- summary$program_id[summary$abs_cor_with_wgcna_module_eigengene < 0.999999]
  stop("Source PC1 does not replay WGCNA module eigengene for: ", paste(bad, collapse = ", "))
}
if (any(summary$source_pc1_vs_unweighted_mean_cor < -ORIENTATION_TOL)) {
  stop("At least one source PC1 orientation remains negative.")
}

weights_out <- file.path(OUT_DIR, "tcga_frozen_source_program_weights_v1.tsv")
summary_out <- file.path(OUT_DIR, "tcga_frozen_source_program_summary_v1.tsv")
scores_out <- file.path(OUT_DIR, "tcga_frozen_source_program_scores_v1.tsv")

write.table(weights, weights_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(summary, summary_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(score_table, scores_out, sep = "\t", row.names = FALSE, quote = FALSE)

result <- list(
  result_id = "paper4-tcbb-frozen-source-programs-v1",
  script_version = SCRIPT_VERSION,
  status = "FROZEN_SOURCE_PROGRAMS",
  scientific_guard = list(
    scanb_opened = FALSE,
    metabric_opened = FALSE,
    gse239948_opened = FALSE,
    clinical_outcomes_loaded = FALSE,
    target_preservation_statistics_calculated = FALSE
  ),
  representation_contract = CONTRACT_JSON,
  programs = lapply(seq_len(nrow(summary)), function(i) {
    as.list(summary[i, , drop = FALSE])
  }),
  primary_source_coherence_for_later_matching = "source_coherence_mean_abs_r"
)

json_out <- file.path(OUT_DIR, "tcga_frozen_source_programs_v1.json")
jsonlite::write_json(result, json_out, pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\n[3/3] Creating SHA-256 manifest ...\n")
manifest_files <- c(
  CONTRACT_JSON,
  SOURCE_MATRIX,
  NETWORK_RDS,
  MEMBERSHIP_TSV,
  weights_out,
  summary_out,
  scores_out,
  json_out,
  list.files(EDGE_DIR, full.names = TRUE, pattern = "\\.rds$")
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
manifest_out <- file.path(OUT_DIR, "tcga_frozen_source_programs_sha256_manifest_v1.tsv")
write.table(manifest, manifest_out, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n", strrep("=", 124), "\n", sep = "")
cat("02g FROZEN TCGA SOURCE PROGRAMS: PASS\n")
cat(strrep("=", 124), "\n", sep = "")
cat("Frozen programs: ", nrow(summary), "\n", sep = "")
cat("Frozen weighted genes across programs: ", nrow(weights), "\n", sep = "")
cat("\nProgram summary:\n")
for (i in seq_len(nrow(summary))) {
  cat(
    "  ", sprintf("%-10s", summary$program_id[i]),
    " n=", sprintf("%4d", summary$n_genes[i]),
    " PC1var=", sprintf("%.3f", summary$source_pc1_variance_explained[i]),
    " mean|r|=", sprintf("%.3f", summary$source_coherence_mean_abs_r[i]),
    " mean_r=", sprintf("%.3f", summary$source_edge_mean_signed_r[i]),
    " signScore~PC1=", sprintf("%.3f", summary$source_signed_score_vs_pc1_cor[i]),
    "\n",
    sep = ""
  )
}
cat("\nPrimary coherence metric for later null matching: mean absolute source Pearson edge correlation\n")
cat("No target dataset was opened.\n")
cat("No target preservation statistic was calculated.\n")
cat("\nOutputs:\n")
cat("  ", weights_out, "\n", sep = "")
cat("  ", summary_out, "\n", sep = "")
cat("  ", scores_out, "\n", sep = "")
cat("  ", json_out, "\n", sep = "")
cat("  ", manifest_out, "\n", sep = "")
cat("  ", EDGE_DIR, "\n", sep = "")
cat(strrep("=", 124), "\n", sep = "")
