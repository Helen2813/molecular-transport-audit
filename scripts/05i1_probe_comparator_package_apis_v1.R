options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "05i1-probe-comparator-package-apis-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"

SOURCE_MODULE_JSON <- file.path(
  DATA_ROOT,
  "paper4_tcbb_tcga_source_modules_v1",
  "tcga_source_modules_frozen_v1.json"
)
SOURCE_MEMBERSHIP <- file.path(
  DATA_ROOT,
  "paper4_tcbb_tcga_source_modules_v1",
  "tcga_source_module_membership_frozen_v1.tsv"
)
CLASSIFICATION <- file.path(
  DATA_ROOT,
  "paper4_tcbb_final_mapping_specificity_null_v2",
  "primary_pooled_final_classification_v2.tsv"
)

cat(strrep("=", 148), "\n", sep = "")
cat("Paper 4 / TCBB - probe exact NetRep and WGCNA comparator APIs before comparator execution\n")
cat(strrep("=", 148), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n\n", sep = "")

cat("Scientific guard:\n")
cat("  Target expression values read:                    NO\n")
cat("  NetRep statistics calculated:                     NO\n")
cat("  WGCNA modulePreservation calculated:              NO\n")
cat("  Comparator settings frozen here:                  NO\n")
cat("  Operation: inspect installed package APIs/help + frozen source metadata only\n")
cat(strrep("=", 148), "\n\n", sep = "")

required_files <- c(SOURCE_MODULE_JSON, SOURCE_MEMBERSHIP, CLASSIFICATION)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
}

if (!requireNamespace("NetRep", quietly = TRUE)) {
  stop("R package 'NetRep' is not installed.")
}
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("R package 'WGCNA' is not installed.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is not installed.")
}

cat("[1/5] Installed package versions\n")
cat("  R:      ", R.version.string, "\n", sep = "")
cat("  NetRep: ", as.character(packageVersion("NetRep")), "\n", sep = "")
cat("  WGCNA:  ", as.character(packageVersion("WGCNA")), "\n", sep = "")

cat("\n[2/5] NetRep exported objects and likely comparator entry points\n")
exports <- sort(getNamespaceExports("NetRep"))
cat(paste(exports, collapse = "\n"), "\n")

candidate_names <- c(
  "netrep",
  "networkProperties",
  "modulePreservation",
  "plotModule",
  "plotModules",
  "preservation"
)

for (nm in candidate_names) {
  if (exists(nm, envir = asNamespace("NetRep"), inherits = FALSE)) {
    obj <- get(nm, envir = asNamespace("NetRep"))
    if (is.function(obj)) {
      cat("\n--- NetRep::", nm, " formals ---\n", sep = "")
      print(formals(obj))
    }
  }
}

cat("\n[3/5] WGCNA::modulePreservation formals\n")
print(formals(WGCNA::modulePreservation))

cat("\n[4/5] Local package help excerpts\n")

print_help <- function(topic, package, max_lines = 220L) {
  cat("\n", strrep("-", 120), "\n", sep = "")
  cat("HELP: ", package, "::", topic, "\n", sep = "")
  cat(strrep("-", 120), "\n", sep = "")

  h <- tryCatch(
    utils::help(topic, package = package),
    error = function(e) NULL
  )
  if (is.null(h) || length(h) == 0) {
    cat("<help topic not found>\n")
    return(invisible(NULL))
  }

  rd <- tryCatch(
    utils:::.getHelpFile(h),
    error = function(e) NULL
  )
  if (is.null(rd)) {
    cat("<could not read help file>\n")
    return(invisible(NULL))
  }

  txt <- capture.output(tools::Rd2txt(rd))
  if (length(txt) > max_lines) {
    txt <- c(
      txt[seq_len(max_lines)],
      paste0("... <truncated after ", max_lines, " lines>")
    )
  }
  cat(paste(txt, collapse = "\n"), "\n")
}

for (topic in c("netrep", "networkProperties")) {
  if (exists(topic, envir = asNamespace("NetRep"), inherits = FALSE)) {
    print_help(topic, "NetRep")
  }
}
print_help("modulePreservation", "WGCNA")

cat("\n[5/5] Frozen source-module and target-assessability metadata\n")

module_meta <- jsonlite::read_json(
  SOURCE_MODULE_JSON,
  simplifyVector = TRUE
)
cat("  source-module status: ",
    as.character(module_meta$status), "\n", sep = "")
cat("  non-grey source modules: ",
    as.integer(module_meta$wgcna$non_grey_modules), "\n", sep = "")

membership <- read.delim(
  SOURCE_MEMBERSHIP,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cat("  membership rows: ", nrow(membership), "\n", sep = "")
cat("  membership columns: ",
    paste(colnames(membership), collapse = ", "), "\n", sep = "")

non_grey <- membership[membership$program_id != "GREY", , drop = FALSE]
sizes <- aggregate(
  Hugo_Symbol ~ program_id + wgcna_label + wgcna_color,
  data = non_grey,
  FUN = length
)
colnames(sizes)[colnames(sizes) == "Hugo_Symbol"] <- "n_source_genes"
sizes <- sizes[order(sizes$program_id), ]
print(sizes, row.names = FALSE)

cls <- read.delim(
  CLASSIFICATION,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

wanted_cols <- intersect(
  c(
    "target",
    "program_id",
    "n_evaluable_genes",
    "coverage",
    "primary_assessable",
    "primary_classification"
  ),
  colnames(cls)
)

cat("\n  corrected primary target assessability/classification columns:\n")
print(
  cls[, wanted_cols, drop = FALSE],
  row.names = FALSE
)

cat("\n", strrep("=", 148), "\n", sep = "")
cat("05i1 COMPARATOR PACKAGE-API PROBE: COMPLETE\n")
cat(strrep("=", 148), "\n", sep = "")
cat("No comparator statistic was calculated.\n")
cat("Important provenance: 05i0 found NO pre-primary frozen NetRep/modulePreservation analysis contract.\n")
cat("The next step must therefore be labeled a post-primary comparator benchmark, with its settings\n")
cat("frozen before any NetRep or modulePreservation statistic is inspected.\n")
cat(strrep("=", 148), "\n", sep = "")
