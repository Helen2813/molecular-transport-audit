options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "05i-run-postprimary-comparator-benchmark-v4-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"

CONTRACT_JSON <- file.path(
  DATA_ROOT,
  "paper4_tcbb_postprimary_comparator_benchmark_contract_v1",
  "postprimary_comparator_benchmark_contract_v1.json"
)

INTERPRETATION_CONTRACT_JSON <- file.path(
  DATA_ROOT,
  "paper4_tcbb_comparator_interpretation_headtohead_contract_v1",
  "comparator_interpretation_headtohead_contract_v1.json"
)

NETREP_PERM_AMENDMENT_JSON <- file.path(
  DATA_ROOT,
  "paper4_tcbb_netrep_permutation_amendment_v1",
  "netrep_permutation_amendment_v1.json"
)

SOURCE_MATRIX <- file.path(
  DATA_ROOT,
  "paper4_tcbb_wgcna_soft_threshold_v1",
  "tcga_source_datExpr_1082x10000_v1.rds"
)
SOURCE_MEMBERSHIP <- file.path(
  DATA_ROOT,
  "paper4_tcbb_tcga_source_modules_v1",
  "tcga_source_module_membership_frozen_v1.tsv"
)

SCANB_EXPR <- file.path(
  DATA_ROOT,
  "SCANB_GSE96058",
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
SCANB_PRIMARY <- file.path(
  DATA_ROOT,
  "paper4_tcbb_scanb_pairing_contract_v1",
  "scanb_primary_profiles_frozen_v1.tsv"
)
SCANB_EVAL <- file.path(
  DATA_ROOT,
  "paper4_tcbb_exact_variation_evaluability_correction_v1",
  "scanb_target_exact_variation_evaluability_v1.tsv"
)

METABRIC_EXPR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_input_audit_v1",
  "staged_continuous_inputs",
  "METABRIC",
  "data_mrna_illumina_microarray.txt"
)
METABRIC_EVAL <- file.path(
  DATA_ROOT,
  "paper4_tcbb_exact_variation_evaluability_correction_v1",
  "metabric_target_exact_variation_evaluability_v1.tsv"
)

CLASSIFICATION <- file.path(
  DATA_ROOT,
  "paper4_tcbb_final_mapping_specificity_null_v2",
  "primary_pooled_final_classification_v2.tsv"
)

OUT_DIR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_postprimary_comparator_benchmark_v4"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

NETREP_PERM <- 10000L
WGCNA_PERM <- 200L
BASE_SEED <- 20260919L

TARGETS <- list(
  SCANB_GSE96058 = list(
    code = 1L,
    n = 3273L,
    expected_genes = 9220L
  ),
  METABRIC = list(
    code = 2L,
    n = 1980L,
    expected_genes = 8485L
  )
)

canon_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[[:space:]]+", " ", x)
  x
}

safe_scale <- function(x, label) {
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(label, ": non-finite values reached scaling.")
  }
  mu <- colMeans(x)
  s <- apply(x, 2, sd)
  if (any(!is.finite(s)) || any(s <= 0)) {
    bad <- colnames(x)[!is.finite(s) | s <= 0]
    stop(
      label,
      ": exact/non-finite constant gene(s) reached scaling: ",
      paste(head(bad, 20), collapse = ", ")
    )
  }
  out <- sweep(x, 2, mu, "-")
  out <- sweep(out, 2, s, "/")
  colnames(out) <- colnames(x)
  rownames(out) <- rownames(x)
  out
}

aggregate_duplicate_rows <- function(values, symbols) {
  if (nrow(values) != length(symbols)) {
    stop("Duplicate-aggregation input length mismatch.")
  }
  if (any(symbols == "")) {
    stop("Blank gene symbol reached duplicate aggregation.")
  }
  summed <- rowsum(
    values,
    group = symbols,
    reorder = FALSE
  )
  first_groups <- rownames(summed)
  counts <- tabulate(
    match(symbols, first_groups),
    nbins = length(first_groups)
  )
  averaged <- summed / counts
  rownames(averaged) <- first_groups
  averaged
}

load_scanb_target <- function(common_genes) {
  primary <- read.delim(
    SCANB_PRIMARY,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  titles <- as.character(primary$primary_title)
  if (length(titles) != 3273L) {
    stop("SCAN-B primary manifest is not 3,273 profiles.")
  }

  # v4 computational bugfix only:
  # The leading gene-ID field can have a blank header. Work by POSITION,
  # not by column name, and separately audit row names as an identifier source.
  raw <- read.csv(
    gzfile(SCANB_EXPR),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  sample_pos <- match(titles, colnames(raw))
  if (anyNA(sample_pos)) {
    missing_titles <- titles[is.na(sample_pos)]
    stop(
      "SCAN-B expression missing primary title(s): ",
      paste(head(missing_titles, 20), collapse = ", ")
    )
  }

  candidate_pos <- setdiff(seq_len(ncol(raw)), unique(sample_pos))
  candidate_labels <- character(0)
  candidate_vectors <- list()
  candidate_unique_overlap <- integer(0)
  candidate_row_hits <- integer(0)

  if (length(candidate_pos) > 0) {
    for (j in candidate_pos) {
      vec <- canon_symbol(raw[[j]])
      nm <- colnames(raw)[j]
      shown_name <- if (is.na(nm) || nm == "") "<blank>" else nm
      label <- paste0("COLUMN[", j, "] name=", shown_name)

      candidate_labels <- c(candidate_labels, label)
      candidate_vectors[[length(candidate_vectors) + 1L]] <- vec
      candidate_unique_overlap <- c(
        candidate_unique_overlap,
        length(intersect(unique(vec), common_genes))
      )
      candidate_row_hits <- c(
        candidate_row_hits,
        sum(vec %in% common_genes)
      )
    }
  }

  rn <- canon_symbol(rownames(raw))
  candidate_labels <- c(candidate_labels, "ROW_NAMES")
  candidate_vectors[[length(candidate_vectors) + 1L]] <- rn
  candidate_unique_overlap <- c(
    candidate_unique_overlap,
    length(intersect(unique(rn), common_genes))
  )
  candidate_row_hits <- c(
    candidate_row_hits,
    sum(rn %in% common_genes)
  )

  cat("  SCAN-B identifier candidates:\n")
  for (k in seq_along(candidate_vectors)) {
    cat(
      "    ", candidate_labels[k],
      ": unique overlap=", candidate_unique_overlap[k],
      "/", length(common_genes),
      "; matching rows=", candidate_row_hits[k],
      "\n",
      sep = ""
    )
  }

  best <- which.max(candidate_unique_overlap)
  symbols <- candidate_vectors[[best]]

  cat(
    "  SCAN-B identifier source selected: ",
    candidate_labels[best],
    "; unique overlap=", candidate_unique_overlap[best],
    "/", length(common_genes),
    "\n",
    sep = ""
  )

  if (candidate_unique_overlap[best] != length(common_genes)) {
    stop(
      "SCAN-B identifier replay mismatch: best candidate '",
      candidate_labels[best],
      "' covers ", candidate_unique_overlap[best],
      " of ", length(common_genes),
      " frozen corrected common genes."
    )
  }

  keep <- symbols %in% common_genes
  values <- as.matrix(raw[keep, sample_pos, drop = FALSE])
  colnames(values) <- titles
  storage.mode(values) <- "double"

  fpkm <- pmax(2^values - 0.1, 0)
  transformed <- log2(fpkm + 1)
  rm(values, fpkm)
  gc(verbose = FALSE)

  agg <- aggregate_duplicate_rows(
    transformed,
    symbols[keep]
  )
  rm(transformed, raw)
  gc(verbose = FALSE)

  missing <- setdiff(common_genes, rownames(agg))
  if (length(missing) > 0) {
    stop(
      "SCAN-B common-gene reconstruction missing after positional parse: ",
      paste(head(missing, 20), collapse = ", ")
    )
  }

  x <- t(agg[common_genes, , drop = FALSE])
  colnames(x) <- common_genes
  rownames(x) <- titles

  if (nrow(x) != 3273L || ncol(x) != length(common_genes)) {
    stop(
      "SCAN-B reconstructed matrix shape mismatch: ",
      nrow(x), " x ", ncol(x)
    )
  }

  x
}

load_metabric_target <- function(common_genes) {
  raw <- read.delim(
    METABRIC_EXPR,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (!("Hugo_Symbol" %in% colnames(raw))) {
    stop("METABRIC Hugo_Symbol column not found.")
  }

  symbols <- canon_symbol(raw$Hugo_Symbol)
  sample_cols <- setdiff(
    colnames(raw),
    c("Hugo_Symbol", "Entrez_Gene_Id")
  )

  if (length(sample_cols) != 1980L) {
    stop(
      "Expected 1,980 METABRIC expression samples; found ",
      length(sample_cols)
    )
  }

  keep <- symbols %in% common_genes
  values <- as.matrix(raw[keep, sample_cols, drop = FALSE])
  storage.mode(values) <- "double"

  agg <- aggregate_duplicate_rows(
    values,
    symbols[keep]
  )
  rm(values, raw)
  gc(verbose = FALSE)

  missing <- setdiff(common_genes, rownames(agg))
  if (length(missing) > 0) {
    stop(
      "METABRIC common-gene reconstruction missing: ",
      paste(head(missing, 20), collapse = ", ")
    )
  }

  x <- t(agg[common_genes, , drop = FALSE])
  colnames(x) <- common_genes
  rownames(x) <- sample_cols
  x
}

load_target_evaluable_symbols <- function(target, source_order) {
  path <- if (target == "SCANB_GSE96058") SCANB_EVAL else METABRIC_EVAL

  ev <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required <- c(
    "Hugo_Symbol",
    "statistically_evaluable_corrected"
  )
  if (!all(required %in% colnames(ev))) {
    stop(
      target,
      ": corrected evaluability manifest missing required columns."
    )
  }

  good <- as.character(
    ev$statistically_evaluable_corrected
  ) %in% c("1", "TRUE", "True", "true")

  symbols <- canon_symbol(ev$Hugo_Symbol[good])
  common <- source_order[source_order %in% symbols]

  expected <- TARGETS[[target]]$expected_genes
  if (length(common) != expected) {
    stop(
      target,
      ": expected ",
      expected,
      " corrected common genes, found ",
      length(common)
    )
  }
  common
}

signed_adjacency_power8 <- function(cor_mat) {
  cor_mat <- pmax(pmin(cor_mat, 1), -1)
  adj <- ((1 + cor_mat) / 2)^8
  diag(adj) <- 1
  dimnames(adj) <- dimnames(cor_mat)
  adj
}

netrep_result_pair <- function(res) {
  if (!is.list(res)) {
    stop("NetRep output is not a list.")
  }

  if (!is.null(res$observed) && !is.null(res$p.values)) {
    return(res)
  }

  # simplify=FALSE: discovery -> test -> result.
  for (a in res) {
    if (!is.list(a)) next
    for (b in a) {
      if (
        is.list(b) &&
        !is.null(b$observed) &&
        !is.null(b$p.values)
      ) {
        return(b)
      }
    }
  }

  stop("Could not locate observed/p.values pair in NetRep result.")
}

as_netrep_tables <- function(pair, target) {
  obs <- as.data.frame(pair$observed, check.names = FALSE)
  p <- as.data.frame(pair$p.values, check.names = FALSE)

  if (is.null(rownames(obs))) {
    stop("NetRep observed matrix has no module row names.")
  }
  if (!identical(rownames(obs), rownames(p))) {
    stop("NetRep observed and p-value module rows differ.")
  }

  modules <- rownames(obs)

  long_obs <- do.call(
    rbind,
    lapply(
      colnames(obs),
      function(stat) {
        data.frame(
          target = target,
          program_id = modules,
          statistic = stat,
          observed = as.numeric(obs[[stat]]),
          raw_p = if (stat %in% colnames(p)) {
            as.numeric(p[[stat]])
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      }
    )
  )

  long_obs$bh_q <- NA_real_
  for (stat in unique(long_obs$statistic)) {
    idx <- which(
      long_obs$statistic == stat &
      is.finite(long_obs$raw_p)
    )
    if (length(idx) > 0) {
      long_obs$bh_q[idx] <- p.adjust(
        long_obs$raw_p[idx],
        method = "BH"
      )
    }
  }

  long_obs
}

find_table_with_col <- function(x, wanted) {
  if (is.matrix(x) || is.data.frame(x)) {
    if (wanted %in% colnames(x)) {
      return(as.data.frame(x, check.names = FALSE))
    }
  }

  if (is.list(x)) {
    for (el in x) {
      hit <- find_table_with_col(el, wanted)
      if (!is.null(hit)) return(hit)
    }
  }

  NULL
}

wgcna_category <- function(z) {
  out <- rep(NA_character_, length(z))
  out[is.finite(z) & z < 2] <- "NO_EVIDENCE_ZSUMMARY_LT2"
  out[is.finite(z) & z >= 2 & z < 10] <- "MODERATE_ZSUMMARY_2_TO_LT10"
  out[is.finite(z) & z >= 10] <- "STRONG_ZSUMMARY_GE10"
  out
}

cat(strrep("=", 152), "\n", sep = "")
cat("Paper 4 / TCBB - post-primary NetRep + WGCNA comparator benchmark\n")
cat(strrep("=", 152), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n\n", sep = "")

cat("Frozen comparator role:\n")
cat("  Prespecified primary MTA analysis:                 NO — comparator only\n")
cat("  Comparator contract frozen before comparator stats:YES\n")
cat("  NetRep permutations:                               ", NETREP_PERM, "\n", sep = "")
cat("  WGCNA modulePreservation permutations:             ", WGCNA_PERM, "\n", sep = "")
cat("  Primary MTA classifications changed:               NO\n")
cat("  v4 SCAN-B parser:                                  positional + row-name identifier audit\n")
cat("  05i3 interpretation matrix/head-to-head frozen:    YES\n")
cat("  05i4 NetRep permutation amendment frozen:          10,000 permutations\n")
cat(strrep("=", 152), "\n\n", sep = "")

required <- c(
  CONTRACT_JSON,
  INTERPRETATION_CONTRACT_JSON,
  NETREP_PERM_AMENDMENT_JSON,
  SOURCE_MATRIX,
  SOURCE_MEMBERSHIP,
  SCANB_EXPR,
  SCANB_PRIMARY,
  SCANB_EVAL,
  METABRIC_EXPR,
  METABRIC_EVAL,
  CLASSIFICATION
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required file(s):\n", paste(missing, collapse = "\n"))
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required.")
}
if (!requireNamespace("NetRep", quietly = TRUE)) {
  stop("R package 'NetRep' is required.")
}
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("R package 'WGCNA' is required.")
}

contract <- jsonlite::read_json(
  CONTRACT_JSON,
  simplifyVector = TRUE
)
if (as.character(contract$status) !=
    "FROZEN_POST_PRIMARY_BEFORE_FIRST_COMPARATOR_STATISTIC") {
  stop("05i2 comparator contract has unexpected status.")
}

interpretation_contract <- jsonlite::read_json(
  INTERPRETATION_CONTRACT_JSON,
  simplifyVector = TRUE
)
if (as.character(interpretation_contract$status) !=
    "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT") {
  stop("05i3 comparator-interpretation contract has unexpected status.")
}

perm_amendment <- jsonlite::read_json(
  NETREP_PERM_AMENDMENT_JSON,
  simplifyVector = TRUE
)
if (as.character(perm_amendment$status) !=
    "FROZEN_BEFORE_FIRST_COMPARATOR_STATISTIC") {
  stop("05i4 NetRep permutation amendment has unexpected status.")
}
if (as.integer(perm_amendment$NetRep$amended_nPerm) != NETREP_PERM) {
  stop("05i4 amendment does not match runner NetRep nPerm.")
}

if (as.character(packageVersion("NetRep")) != "1.2.10") {
  stop("NetRep version drift.")
}
if (as.character(packageVersion("WGCNA")) != "1.74") {
  stop("WGCNA version drift.")
}

source_full <- readRDS(SOURCE_MATRIX)
if (!is.matrix(source_full) && !is.data.frame(source_full)) {
  stop("Frozen TCGA source matrix has unexpected type.")
}
source_full <- as.matrix(source_full)
storage.mode(source_full) <- "double"

if (!identical(dim(source_full), c(1082L, 10000L))) {
  stop(
    "Unexpected frozen TCGA source dimensions: ",
    paste(dim(source_full), collapse = " x ")
  )
}

source_order <- canon_symbol(colnames(source_full))
colnames(source_full) <- source_order

membership <- read.delim(
  SOURCE_MEMBERSHIP,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
membership$Hugo_Symbol <- canon_symbol(membership$Hugo_Symbol)

if (nrow(membership) != 10000L) {
  stop("Frozen source membership is not 10,000 rows.")
}

midx <- match(source_order, membership$Hugo_Symbol)
if (anyNA(midx)) {
  stop("Source membership does not cover frozen source matrix.")
}
membership <- membership[midx, , drop = FALSE]

classification <- read.delim(
  CLASSIFICATION,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

all_programs <- sprintf("TCGA_M%03d", 1:12)

netrep_all <- list()
wgcna_all <- list()
integrated_all <- list()

for (target in names(TARGETS)) {
  cfg <- TARGETS[[target]]

  cat("\n", strrep("=", 152), "\n", sep = "")
  cat("TARGET: ", target, "\n", sep = "")
  cat(strrep("=", 152), "\n", sep = "")

  target_dir <- file.path(OUT_DIR, target)
  dir.create(
    target_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  common_genes <- load_target_evaluable_symbols(
    target,
    source_order
  )
  cat("  corrected common genes: ", length(common_genes), "\n", sep = "")

  sidx <- match(common_genes, source_order)
  source_x <- source_full[, sidx, drop = FALSE]
  colnames(source_x) <- common_genes

  target_x <- if (target == "SCANB_GSE96058") {
    load_scanb_target(common_genes)
  } else {
    load_metabric_target(common_genes)
  }

  if (nrow(target_x) != cfg$n) {
    stop(
      target,
      ": expected ",
      cfg$n,
      " target samples, found ",
      nrow(target_x)
    )
  }
  if (!identical(colnames(source_x), colnames(target_x))) {
    stop(target, ": source/target common-gene order mismatch.")
  }

  source_z <- safe_scale(
    source_x,
    paste0(target, " source")
  )
  target_z <- safe_scale(
    target_x,
    paste0(target, " target")
  )

  assign_idx <- match(
    common_genes,
    membership$Hugo_Symbol
  )
  module_assignment <- as.character(
    membership$program_id[assign_idx]
  )
  names(module_assignment) <- common_genes

  module_color <- as.character(
    membership$wgcna_color[assign_idx]
  )
  names(module_color) <- common_genes

  if (anyNA(module_assignment) || anyNA(module_color)) {
    stop(target, ": missing frozen source module labels/colors.")
  }

  # ---------------------------------------------------------------------------
  # NetRep
  # ---------------------------------------------------------------------------
  netrep_rds <- file.path(
    target_dir,
    "netrep_modulePreservation_v4.rds"
  )
  netrep_tsv <- file.path(
    target_dir,
    "netrep_modulePreservation_long_v4.tsv"
  )

  if (!file.exists(netrep_rds)) {
    cat("\n[NetRep] Building Pearson correlation matrices ...\n")
    source_cor <- cor(
      source_z,
      use = "pairwise.complete.obs",
      method = "pearson"
    )
    target_cor <- cor(
      target_z,
      use = "pairwise.complete.obs",
      method = "pearson"
    )
    dimnames(source_cor) <- list(common_genes, common_genes)
    dimnames(target_cor) <- list(common_genes, common_genes)

    cat("[NetRep] Building signed power-8 adjacency matrices ...\n")
    source_net <- signed_adjacency_power8(source_cor)
    target_net <- signed_adjacency_power8(target_cor)

    set.seed(BASE_SEED + cfg$code * 1000L + 11L)

    cat(
      "[NetRep] Running ",
      NETREP_PERM,
      " permutations on all 12 frozen source modules ...\n",
      sep = ""
    )

    netrep_res <- NetRep::modulePreservation(
      network = list(
        source = source_net,
        target = target_net
      ),
      data = list(
        source = source_z,
        target = target_z
      ),
      correlation = list(
        source = source_cor,
        target = target_cor
      ),
      moduleAssignments = list(
        source = module_assignment
      ),
      modules = list(
        source = all_programs
      ),
      backgroundLabel = "GREY",
      discovery = "source",
      test = list(source = "target"),
      selfPreservation = FALSE,
      nThreads = 1L,
      nPerm = NETREP_PERM,
      null = "overlap",
      alternative = "greater",
      simplify = FALSE,
      verbose = TRUE
    )

    saveRDS(
      netrep_res,
      netrep_rds,
      compress = "xz"
    )

    pair <- netrep_result_pair(netrep_res)
    netrep_long <- as_netrep_tables(
      pair,
      target
    )
    write.table(
      netrep_long,
      netrep_tsv,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    rm(
      source_cor,
      target_cor,
      source_net,
      target_net,
      netrep_res,
      pair
    )
    gc(verbose = FALSE)
  } else {
    cat("\n[NetRep] Completed RDS found — reusing.\n")
    netrep_res <- readRDS(netrep_rds)
    pair <- netrep_result_pair(netrep_res)
    netrep_long <- as_netrep_tables(pair, target)
    if (!file.exists(netrep_tsv)) {
      write.table(
        netrep_long,
        netrep_tsv,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
    rm(netrep_res, pair)
    gc(verbose = FALSE)
  }

  netrep_all[[target]] <- netrep_long

  # ---------------------------------------------------------------------------
  # WGCNA modulePreservation
  # ---------------------------------------------------------------------------
  wgcna_rds <- file.path(
    target_dir,
    "wgcna_modulePreservation_v4.rds"
  )
  wgcna_tsv <- file.path(
    target_dir,
    "wgcna_modulePreservation_summary_v4.tsv"
  )

  if (!file.exists(wgcna_rds)) {
    cat(
      "\n[WGCNA] Running modulePreservation with ",
      WGCNA_PERM,
      " permutations ...\n",
      sep = ""
    )

    multiData <- list(
      source = list(data = source_z),
      target = list(data = target_z)
    )
    multiColor <- list(
      source = module_color
    )

    wgcna_res <- WGCNA::modulePreservation(
      multiData = multiData,
      multiColor = multiColor,
      multiWeights = NULL,
      dataIsExpr = TRUE,
      networkType = "signed",
      corFnc = "cor",
      corOptions = "use = 'p'",
      referenceNetworks = 1L,
      testNetworks = list(2L),
      nPermutations = WGCNA_PERM,
      includekMEallInSummary = FALSE,
      restrictSummaryForGeneralNetworks = TRUE,
      calculateQvalue = FALSE,
      randomSeed = BASE_SEED + cfg$code * 1000L + 21L,
      maxGoldModuleSize = 5000L,
      maxModuleSize = 5000L,
      quickCor = 0L,
      ccTupletSize = 2L,
      calculateCor.kIMall = FALSE,
      calculateClusterCoeff = FALSE,
      useInterpolation = FALSE,
      checkData = TRUE,
      greyName = "grey",
      goldName = "gold",
      savePermutedStatistics = FALSE,
      loadPermutedStatistics = FALSE,
      plotInterpolation = FALSE,
      discardInvalidOutput = TRUE,
      parallelCalculation = FALSE,
      verbose = 2L,
      indent = 0L
    )

    saveRDS(
      wgcna_res,
      wgcna_rds,
      compress = "xz"
    )
  } else {
    cat("\n[WGCNA] Completed RDS found — reusing.\n")
    wgcna_res <- readRDS(wgcna_rds)
  }

  ztab <- find_table_with_col(
    wgcna_res,
    "Zsummary.pres"
  )
  rtab <- find_table_with_col(
    wgcna_res,
    "medianRank.pres"
  )

  if (is.null(ztab) || is.null(rtab)) {
    stop(
      target,
      ": could not locate WGCNA Zsummary/medianRank output tables."
    )
  }

  ztab$wgcna_color <- rownames(ztab)
  rtab$wgcna_color <- rownames(rtab)

  zkeep <- intersect(
    c(
      "wgcna_color",
      "moduleSize",
      "Zsummary.pres",
      "Zdensity.pres",
      "Zconnectivity.pres"
    ),
    colnames(ztab)
  )
  rkeep <- intersect(
    c(
      "wgcna_color",
      "moduleSize",
      "medianRank.pres",
      "medianRankDensity.pres",
      "medianRankConnectivity.pres"
    ),
    colnames(rtab)
  )

  wgcna_summary <- merge(
    ztab[, zkeep, drop = FALSE],
    rtab[, rkeep, drop = FALSE],
    by = "wgcna_color",
    all = TRUE,
    suffixes = c(".Z", ".rank")
  )

  color_map <- unique(
    membership[
      membership$program_id != "GREY",
      c("program_id", "wgcna_color"),
      drop = FALSE
    ]
  )

  wgcna_summary <- merge(
    color_map,
    wgcna_summary,
    by = "wgcna_color",
    all.x = TRUE
  )
  wgcna_summary$target <- target

  if ("Zsummary.pres" %in% colnames(wgcna_summary)) {
    wgcna_summary$Zsummary_category <- wgcna_category(
      as.numeric(wgcna_summary$Zsummary.pres)
    )
  }

  write.table(
    wgcna_summary,
    wgcna_tsv,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  wgcna_all[[target]] <- wgcna_summary

  # ---------------------------------------------------------------------------
  # Integrated audit table
  # ---------------------------------------------------------------------------
  cls <- classification[
    classification$target == target,
    ,
    drop = FALSE
  ]

  net_cor <- netrep_long[
    netrep_long$statistic %in% c("cor.cor", "cor.contrib"),
    ,
    drop = FALSE
  ]
  if (nrow(net_cor) > 0) {
    net_wide <- reshape(
      net_cor[
        ,
        c(
          "program_id",
          "statistic",
          "observed",
          "raw_p",
          "bh_q"
        ),
        drop = FALSE
      ],
      idvar = "program_id",
      timevar = "statistic",
      direction = "wide"
    )
  } else {
    net_wide <- data.frame(program_id = all_programs)
  }

  keep_cls <- intersect(
    c(
      "program_id",
      "primary_assessable",
      "primary_classification",
      "rho_edge",
      "rho_load",
      "coverage",
      "n_evaluable_genes"
    ),
    colnames(cls)
  )

  integrated <- merge(
    cls[, keep_cls, drop = FALSE],
    net_wide,
    by = "program_id",
    all.x = TRUE
  )
  integrated <- merge(
    integrated,
    wgcna_summary,
    by = "program_id",
    all.x = TRUE,
    suffixes = c("", ".wgcna")
  )

  integrated$target <- target
  integrated$comparator_role <-
    "POST_PRIMARY_NO_RECLASSIFICATION"

  integrated_path <- file.path(
    target_dir,
    "integrated_mta_netrep_wgcna_comparison_v4.tsv"
  )
  write.table(
    integrated,
    integrated_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  integrated_all[[target]] <- integrated

  rm(
    source_x,
    target_x,
    source_z,
    target_z,
    wgcna_res,
    ztab,
    rtab
  )
  gc(verbose = FALSE)

  cat(
    "\nTARGET COMPLETE: ",
    target,
    "\n",
    sep = ""
  )
}

cat("\n[FINAL] Writing cross-target comparator summaries ...\n")

netrep_long_all <- do.call(
  rbind,
  netrep_all
)
wgcna_summary_all <- do.call(
  rbind,
  wgcna_all
)
integrated_all_df <- do.call(
  rbind,
  integrated_all
)

netrep_all_path <- file.path(
  OUT_DIR,
  "netrep_all_targets_long_v4.tsv"
)
wgcna_all_path <- file.path(
  OUT_DIR,
  "wgcna_modulePreservation_all_targets_v4.tsv"
)
integrated_all_path <- file.path(
  OUT_DIR,
  "mta_netrep_wgcna_integrated_all_targets_v4.tsv"
)

write.table(
  netrep_long_all,
  netrep_all_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  wgcna_summary_all,
  wgcna_all_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  integrated_all_df,
  integrated_all_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

master <- list(
  script_version = SCRIPT_VERSION,
  status = "POSTPRIMARY_COMPARATOR_BENCHMARK_COMPLETE",
  provenance = list(
    comparator_contract =
      "Comparator benchmark frozen post-primary before first comparator statistic.",
    v2_bugfix =
      paste(
        "Runner v1 stopped before any comparator statistic because base R",
        "read.csv() auto-consumed the leading SCAN-B gene-identifier field",
        "as row names. Runner v2 forces row.names=NULL and validates the",
        "identifier column against the frozen corrected common-gene universe.",
        "Scientific comparator contract is unchanged."
      ),
    v1_comparator_statistics_calculated = FALSE
  ),
  targets = names(TARGETS),
  NetRep = list(
    version = as.character(packageVersion("NetRep")),
    permutations = NETREP_PERM,
    all_seven_statistics_reported = TRUE,
    no_custom_composite_classification = TRUE
  ),
  WGCNA = list(
    version = as.character(packageVersion("WGCNA")),
    permutations = WGCNA_PERM,
    headline_outputs = c(
      "Zsummary.pres",
      "medianRank.pres"
    )
  ),
  primary_mta_classification_changed = FALSE,
  interpretation_contract = INTERPRETATION_CONTRACT_JSON,
  netrep_permutation_amendment = NETREP_PERM_AMENDMENT_JSON,
  runner_bugfix_provenance = list(
    v1 = "stopped before comparator statistics: SCAN-B gene identifier interpreted incorrectly",
    v2 = "explicit-column parser prepared but superseded before comparator execution",
    v3 = "stopped before comparator statistics: blank-named gene-ID column accessed unsafely by name",
    v4 = "positional-column plus row-name identifier audit; scientific comparator contract unchanged"
  ),
  interpretation_rules = list(
    netrep_omnibus_binary_preserved_call = FALSE,
    netrep_structural_axis = "cor.cor",
    netrep_loading_axis =
      "cor.contrib interpreted only with avg.contrib BH q < 0.05",
    wgcna_Zsummary_categories =
      "<2 none; 2-<10 moderate; >=10 strong",
    cross_method_pvalues_directly_equated = FALSE,
    mandatory_headtohead_corruption_ladder = TRUE
  ),
  output_files = list(
    NetRep_long = netrep_all_path,
    WGCNA_summary = wgcna_all_path,
    integrated = integrated_all_path
  )
)

master_path <- file.path(
  OUT_DIR,
  "postprimary_comparator_benchmark_v4.json"
)
jsonlite::write_json(
  master,
  master_path,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

cat("\n", strrep("=", 152), "\n", sep = "")
cat("05i POST-PRIMARY COMPARATOR BENCHMARK: COMPLETE\n")
cat(strrep("=", 152), "\n", sep = "")
cat("NetRep:                             COMPLETE for both targets\n")
cat("WGCNA modulePreservation:           COMPLETE for both targets\n")
cat("NetRep statistics retained:         ALL 7\n")
cat("WGCNA headline metrics:             Zsummary + medianRank\n")
cat("Primary MTA classifications changed:NO\n")
cat("\nNetRep:    ", netrep_all_path, "\n", sep = "")
cat("WGCNA:     ", wgcna_all_path, "\n", sep = "")
cat("Integrated:", integrated_all_path, "\n", sep = "")
cat("Master:    ", master_path, "\n", sep = "")
cat(strrep("=", 152), "\n", sep = "")
