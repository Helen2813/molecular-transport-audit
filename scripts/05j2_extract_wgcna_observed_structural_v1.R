options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "05j2-extract-wgcna-observed-structural-v1-no-cli"
DATA_ROOT <- "D:/paper4_tcbb_data"

MEMBERSHIP <- file.path(
  DATA_ROOT,
  "paper4_tcbb_tcga_source_modules_v1",
  "tcga_source_module_membership_frozen_v1.tsv"
)

V7_ROOT <- file.path(
  DATA_ROOT,
  "paper4_tcbb_postprimary_comparator_benchmark_v7"
)

RDS_FILES <- c(
  SCANB_GSE96058 = file.path(
    V7_ROOT,
    "SCANB_GSE96058",
    "wgcna_modulePreservation_v7.rds"
  ),
  METABRIC = file.path(
    V7_ROOT,
    "METABRIC",
    "wgcna_modulePreservation_v7.rds"
  )
)

CONTRACT <- file.path(
  DATA_ROOT,
  "paper4_tcbb_structural_corruption_headtohead_execution_contract_v1",
  "structural_corruption_headtohead_execution_contract_v1.json"
)

OUT_DIR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_structural_corruption_headtohead_execution_contract_v1"
)
OUT_TSV <- file.path(
  OUT_DIR,
  "wgcna_observed_structural_baseline_v1.tsv"
)

cat(strrep("=", 154), "\n", sep = "")
cat("Paper 4 / TCBB - extract already-computed WGCNA raw structural observations\n")
cat(strrep("=", 154), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n\n", sep = "")
cat("Scientific guard:\n")
cat("  WGCNA modulePreservation recalculated:             NO\n")
cat("  Permutations recalculated:                         NO\n")
cat("  Operation: read completed v7 RDS and export raw observed cor.cor/cor.kIM only\n")
cat(strrep("=", 154), "\n\n", sep = "")

for (f in c(MEMBERSHIP, CONTRACT, unname(RDS_FILES))) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required.")
}

contract <- jsonlite::read_json(CONTRACT, simplifyVector = TRUE)
if (as.character(contract$status) !=
    "FROZEN_BEFORE_FIRST_CORRUPTION_HEADTOHEAD_COMPARATOR_STATISTIC") {
  stop("05j2 execution contract has unexpected status.")
}

membership <- read.delim(
  MEMBERSHIP,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

color_map <- unique(
  membership[
    membership$program_id != "GREY",
    c("program_id", "wgcna_color"),
    drop = FALSE
  ]
)

find_structural_table <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    cn <- colnames(x)
    if (!is.null(cn) && all(c("cor.cor", "cor.kIM") %in% cn)) {
      return(as.data.frame(x, check.names = FALSE))
    }
  }
  if (is.list(x)) {
    for (el in x) {
      hit <- find_structural_table(el)
      if (!is.null(hit)) return(hit)
    }
  }
  NULL
}

all_rows <- list()

for (target in names(RDS_FILES)) {
  cat("Reading completed WGCNA RDS for ", target, " ...\n", sep = "")
  res <- readRDS(RDS_FILES[[target]])

  # Prefer the documented preservation$observed object; recursive search is a
  # defensive fallback for minor nesting differences.
  tab <- NULL
  if (!is.null(res$preservation) && !is.null(res$preservation$observed)) {
    tab <- find_structural_table(res$preservation$observed)
  }
  if (is.null(tab)) tab <- find_structural_table(res)
  if (is.null(tab)) {
    stop(target, ": could not find a raw observed WGCNA table containing cor.cor and cor.kIM.")
  }

  tab$wgcna_color <- rownames(tab)
  tab <- merge(
    color_map,
    tab,
    by = "wgcna_color",
    all.x = TRUE
  )
  tab$target <- target

  keep <- c(
    "target",
    "program_id",
    "wgcna_color",
    "moduleSize",
    "cor.cor",
    "cor.kIM"
  )
  keep <- keep[keep %in% colnames(tab)]

  out <- tab[, keep, drop = FALSE]
  out <- out[order(out$program_id), , drop = FALSE]
  all_rows[[target]] <- out
}

final <- do.call(rbind, all_rows)
rownames(final) <- NULL

if (!all(c("target", "program_id", "cor.cor", "cor.kIM") %in% colnames(final))) {
  stop("Extracted WGCNA table is missing required structural columns.")
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.table(
  final,
  OUT_TSV,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\nExtracted rows: ", nrow(final), "\n", sep = "")
print(final, row.names = FALSE)

cat("\n", strrep("=", 154), "\n", sep = "")
cat("05j2 WGCNA OBSERVED STRUCTURAL EXTRACTION: COMPLETE\n")
cat(strrep("=", 154), "\n", sep = "")
cat("No WGCNA preservation statistic was recomputed.\n")
cat("Output: ", OUT_TSV, "\n", sep = "")
cat(strrep("=", 154), "\n", sep = "")
