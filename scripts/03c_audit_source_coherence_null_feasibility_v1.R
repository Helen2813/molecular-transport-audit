options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "03c-audit-source-coherence-null-feasibility-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"

SOURCE_MATRIX <- file.path(
  DATA_ROOT,
  "paper4_tcbb_wgcna_soft_threshold_v1",
  "tcga_source_datExpr_1082x10000_v1.rds"
)
SOURCE_UNIVERSE <- file.path(
  DATA_ROOT,
  "paper4_tcbb_tcga_source_universe_v1",
  "tcga_source_gene_universe_frozen_v1.tsv"
)
SOURCE_WEIGHTS <- file.path(
  DATA_ROOT,
  "paper4_tcbb_frozen_source_programs_v1",
  "tcga_frozen_source_program_weights_v1.tsv"
)
MAPPING_ROOT <- file.path(DATA_ROOT, "paper4_tcbb_target_mapping_v1")
SCANB_MAP <- file.path(MAPPING_ROOT, "scanb_gene_mapping_frozen_v1.tsv")
METABRIC_MAP <- file.path(MAPPING_ROOT, "metabric_gene_mapping_frozen_v1.tsv")
COVERAGE <- file.path(MAPPING_ROOT, "target_program_gene_coverage_v1.tsv")

CONTRACT_JSON <- file.path(
  DATA_ROOT,
  "paper4_tcbb_coherence_null_feasibility_contract_v1",
  "coherence_null_feasibility_contract_v1.json"
)

OUT_DIR <- file.path(DATA_ROOT, "paper4_tcbb_coherence_null_feasibility_audit_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(strrep("=", 124), "\n", sep = "")
cat("Paper 4 / TCBB - audit feasibility of naive source-coherence-matched random panels\n")
cat(strrep("=", 124), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("\nScientific guard:\n")
cat("  Target outcomes loaded:                              NO\n")
cat("  Target clinical characteristics loaded:              NO\n")
cat("  Target expression values loaded:                     NO\n")
cat("  Target correlation matrices calculated:              NO\n")
cat("  Target preservation statistics calculated:           NO\n")
cat("  Source Pearson correlations calculated:              YES\n")
cat("  Target identifier availability used:                 YES\n")
cat(strrep("=", 124), "\n\n", sep = "")

required <- c(
  SOURCE_MATRIX, SOURCE_UNIVERSE, SOURCE_WEIGHTS,
  SCANB_MAP, METABRIC_MAP, COVERAGE, CONTRACT_JSON
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) stop("Missing required file(s):\n", paste(missing, collapse = "\n"))

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed.")

contract <- jsonlite::read_json(CONTRACT_JSON, simplifyVector = TRUE)
N_PANELS <- as.integer(contract$feasibility_design$random_panels_per_target_program)
BASE_SEED <- as.integer(contract$feasibility_design$base_seed)
PRIMARY_ABS_TOL <- as.numeric(contract$diagnostic_tolerances$primary_absolute)
PRIMARY_REL_TOL <- as.numeric(contract$diagnostic_tolerances$primary_relative)
MIN_MATCHES <- as.integer(contract$naive_rejection_feasibility_rule$minimum_matches_out_of_100)

cat("[1/4] Loading frozen source data and identifier-only target maps ...\n")
datExpr <- readRDS(SOURCE_MATRIX)
if (nrow(datExpr) != 1082 || ncol(datExpr) != 10000) stop("Unexpected source matrix dimensions.")

universe <- read.delim(SOURCE_UNIVERSE, sep = "\t", header = TRUE, check.names = FALSE,
                       quote = "", comment.char = "")
weights <- read.delim(SOURCE_WEIGHTS, sep = "\t", header = TRUE, check.names = FALSE,
                      quote = "", comment.char = "")
scanb_map <- read.delim(SCANB_MAP, sep = "\t", header = TRUE, check.names = FALSE,
                        quote = "", comment.char = "")
metabric_map <- read.delim(METABRIC_MAP, sep = "\t", header = TRUE, check.names = FALSE,
                           quote = "", comment.char = "")
coverage <- read.delim(COVERAGE, sep = "\t", header = TRUE, check.names = FALSE,
                       quote = "", comment.char = "")

norm <- function(x) toupper(trimws(as.character(x)))

source_genes <- norm(colnames(datExpr))
if (length(unique(source_genes)) != 10000) stop("Source genes not unique after normalization.")
gene_index <- setNames(seq_along(source_genes), source_genes)

program_genes <- split(norm(weights$Hugo_Symbol), weights$program_id)

target_symbols <- list(
  SCANB_GSE96058 = unique(norm(scanb_map$target_normalized_symbol)),
  METABRIC = unique(norm(metabric_map$target_normalized_symbol))
)
target_symbols <- lapply(target_symbols, function(x) x[nzchar(x)])

cat("  source genes: 10,000\n")
cat("  target mapped symbols: SCAN-B=", length(target_symbols$SCANB_GSE96058),
    ", METABRIC=", length(target_symbols$METABRIC), "\n", sep = "")

cat("\n[2/4] Calculating full SOURCE Pearson correlation matrix ...\n")
cat("  (10,000 x 10,000; target expression is not involved)\n")

# cor() centers/scales internally. Matrix is retained only for this source-side design audit.
source_cor <- cor(datExpr, use = "everything", method = "pearson")
if (any(!is.finite(source_cor))) stop("Non-finite source correlation.")
abs_cor <- abs(source_cor)
rm(source_cor)
diag(abs_cor) <- 0
gc(verbose = FALSE)

mean_abs_coherence <- function(idx) {
  m <- length(idx)
  if (m < 2) return(NA_real_)
  # Sum symmetric matrix, diagonal is zero; divide by m*(m-1).
  sum(abs_cor[idx, idx, drop = FALSE]) / (m * (m - 1))
}

cat("\n[3/4] Running predeclared unconditional-panel feasibility audit ...\n")

detail_rows <- list()
summary_rows <- list()
row_counter <- 1L
detail_counter <- 1L

targets <- c("SCANB_GSE96058", "METABRIC")

for (t_i in seq_along(targets)) {
  target <- targets[t_i]
  avail <- intersect(source_genes, target_symbols[[target]])
  avail_set <- unique(avail)

  cov_t <- coverage[coverage$target == target, , drop = FALSE]
  if (nrow(cov_t) != 12L) stop("Coverage table does not have 12 rows for ", target)

  cat("\n  ", target, ": target-available frozen source universe = ",
      length(avail_set), "\n", sep = "")

  for (p_i in seq_along(program_genes)) {
    pid <- names(program_genes)[p_i]
    pg <- program_genes[[pid]]
    eval_genes <- intersect(pg, avail_set)
    m <- length(eval_genes)

    cov_row <- cov_t[cov_t$program_id == pid, , drop = FALSE]
    if (nrow(cov_row) != 1L) stop("Coverage row missing/duplicated for ", target, " ", pid)
    if (m != as.integer(cov_row$evaluable_genes)) {
      stop("Evaluable-gene replay mismatch for ", target, " ", pid)
    }

    eval_idx <- unname(gene_index[eval_genes])
    target_coh <- mean_abs_coherence(eval_idx)

    # Exclude ALL frozen genes of the tested program, not only evaluable members.
    candidate_genes <- setdiff(avail_set, pg)
    candidate_idx <- unname(gene_index[candidate_genes])

    if (length(candidate_idx) < m) {
      stop("Candidate pool too small for ", target, " ", pid)
    }

    # Deterministic target/program seed schedule.
    seed <- BASE_SEED + (t_i * 1000L) + p_i
    set.seed(seed)

    vals <- numeric(N_PANELS)
    for (b in seq_len(N_PANELS)) {
      idx <- sample(candidate_idx, size = m, replace = FALSE)
      vals[b] <- mean_abs_coherence(idx)

      detail_rows[[detail_counter]] <- data.frame(
        target = target,
        program_id = pid,
        panel_id = b,
        panel_size = m,
        source_coherence_mean_abs_r = vals[b],
        stringsAsFactors = FALSE
      )
      detail_counter <- detail_counter + 1L
    }

    abs_gap <- abs(vals - target_coh)
    rel_gap <- abs_gap / target_coh

    primary_match <- (abs_gap <= PRIMARY_ABS_TOL) & (rel_gap <= PRIMARY_REL_TOL)
    n_match <- sum(primary_match)

    q <- as.numeric(quantile(vals, probs = c(.01, .05, .50, .95, .99),
                             names = FALSE, type = 7))

    percentile <- mean(vals <= target_coh)

    summary_rows[[row_counter]] <- data.frame(
      target = target,
      program_id = pid,
      panel_size = m,
      candidate_pool_size = length(candidate_idx),
      program_source_coherence = target_coh,
      random_q01 = q[1],
      random_q05 = q[2],
      random_median = q[3],
      random_q95 = q[4],
      random_q99 = q[5],
      random_max = max(vals),
      program_percentile_among_unconditional_random = percentile,
      nearest_abs_gap = min(abs_gap),
      nearest_rel_gap = min(rel_gap),
      n_match_abs_0p01 = sum(abs_gap <= 0.01),
      n_match_abs_0p02 = sum(abs_gap <= 0.02),
      n_match_rel_5pct = sum(rel_gap <= 0.05),
      n_match_rel_10pct = sum(rel_gap <= 0.10),
      n_primary_match = n_match,
      naive_rejection_feasible = as.integer(n_match >= MIN_MATCHES),
      primary_assessable = as.integer(cov_row$primary_assessable),
      seed = seed,
      stringsAsFactors = FALSE
    )
    row_counter <- row_counter + 1L

    cat(
      "    ", sprintf("%-10s", pid),
      " m=", sprintf("%4d", m),
      " moduleC=", sprintf("%.3f", target_coh),
      " randomMed=", sprintf("%.3f", median(vals)),
      " randomMax=", sprintf("%.3f", max(vals)),
      " matches=", sprintf("%3d", n_match), "/", N_PANELS,
      " feasible=", ifelse(n_match >= MIN_MATCHES, "YES", "NO"),
      "\n",
      sep = ""
    )
  }
}

summary <- do.call(rbind, summary_rows)
detail <- do.call(rbind, detail_rows)

cat("\n[4/4] Writing feasibility audit ...\n")
summary_out <- file.path(OUT_DIR, "coherence_null_feasibility_summary_v1.tsv")
detail_out <- file.path(OUT_DIR, "coherence_null_unconditional_panel_draws_v1.tsv")
write.table(summary, summary_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(detail, detail_out, sep = "\t", row.names = FALSE, quote = FALSE)

result <- list(
  result_id = "paper4-tcbb-coherence-null-feasibility-audit-v1",
  script_version = SCRIPT_VERSION,
  status = "SOURCE_SIDE_DESIGN_AUDIT_COMPLETE",
  scientific_guard = list(
    target_outcomes_loaded = FALSE,
    target_clinical_characteristics_loaded = FALSE,
    target_expression_values_loaded = FALSE,
    target_preservation_statistics_calculated = FALSE
  ),
  design_contract = CONTRACT_JSON,
  summary = lapply(seq_len(nrow(summary)), function(i) as.list(summary[i, , drop = FALSE])),
  interpretation_rule = paste0(
    "naive rejection sampling is practically feasible only if at least ",
    MIN_MATCHES, "/", N_PANELS,
    " panels meet abs gap <= ", PRIMARY_ABS_TOL,
    " and relative gap <= ", 100 * PRIMARY_REL_TOL, "%"
  )
)

json_out <- file.path(OUT_DIR, "coherence_null_feasibility_audit_v1.json")
jsonlite::write_json(result, json_out, pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\n", strrep("=", 124), "\n", sep = "")
cat("03c SOURCE-COHERENCE NULL FEASIBILITY AUDIT: PASS\n")
cat(strrep("=", 124), "\n", sep = "")
for (target in targets) {
  s <- summary[summary$target == target, , drop = FALSE]
  cat(target, ":\n", sep = "")
  cat("  naive rejection feasible: ",
      sum(s$naive_rejection_feasible == 1), "/12 programs\n", sep = "")
  cat("  primary-assessable programs: ",
      sum(s$primary_assessable == 1), "/12\n", sep = "")
}
cat("\nNo target expression value was loaded.\n")
cat("No target preservation statistic was calculated.\n")
cat("This audit does NOT yet define the final matched-null generator.\n")
cat("\nOutputs:\n")
cat("  ", summary_out, "\n", sep = "")
cat("  ", detail_out, "\n", sep = "")
cat("  ", json_out, "\n", sep = "")
cat(strrep("=", 124), "\n", sep = "")
