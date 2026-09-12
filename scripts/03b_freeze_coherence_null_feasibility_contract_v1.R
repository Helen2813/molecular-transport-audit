options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "03b-freeze-coherence-null-feasibility-contract-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"
OUT_DIR <- file.path(DATA_ROOT, "paper4_tcbb_coherence_null_feasibility_contract_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Diagnostic-only feasibility audit parameters.
N_RANDOM_PANELS_PER_TARGET_PROGRAM <- 100L
BASE_SEED <- 20260912L

ABS_TOLERANCES <- c(0.01, 0.02)
REL_TOLERANCES <- c(0.05, 0.10)

# Predeclared practical rule:
# naive rejection sampling is considered practically feasible only if at least
# 5 of 100 unconditional same-size random panels satisfy BOTH:
#   abs coherence gap <= 0.02
#   relative coherence gap <= 10%
# This does not define the final null; it only determines whether naive rejection
# is a viable generator or whether a source-conditioned sampler is required.
PRIMARY_ABS_TOL <- 0.02
PRIMARY_REL_TOL <- 0.10
MIN_MATCHES_FOR_NAIVE_FEASIBILITY <- 5L

cat(strrep("=", 122), "\n", sep = "")
cat("Paper 4 / TCBB - freeze source-coherence null FEASIBILITY audit contract\n")
cat(strrep("=", 122), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("\nScientific guard:\n")
cat("  Target clinical outcomes loaded:                     NO\n")
cat("  Target clinical characteristics loaded:              NO\n")
cat("  Target expression values loaded:                     NO\n")
cat("  Target correlations/preservation calculated:          NO\n")
cat("  Source correlations permitted:                       YES\n")
cat("  Target identifier availability permitted:            YES\n")
cat("  Final matched-null algorithm selected here:          NO\n")
cat(strrep("=", 122), "\n\n", sep = "")

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed.")

contract <- list(
  contract_id = "paper4-tcbb-coherence-null-feasibility-v1",
  script_version = SCRIPT_VERSION,
  status = "FROZEN_BEFORE_FEASIBILITY_AUDIT",
  scientific_guard = list(
    target_outcomes_loaded = FALSE,
    target_clinical_characteristics_loaded = FALSE,
    target_expression_values_loaded = FALSE,
    target_preservation_statistics_calculated = FALSE
  ),
  feasibility_design = list(
    random_panels_per_target_program = N_RANDOM_PANELS_PER_TARGET_PROGRAM,
    base_seed = BASE_SEED,
    panel_size = "exactly the target-evaluable frozen-program gene count",
    candidate_universe = paste0(
      "frozen 10,000-gene TCGA source universe intersected with the target's ",
      "frozen mapped-symbol universe, excluding ALL genes of the tested frozen program"
    ),
    sampling = "uniform without replacement within each panel",
    coherence = "mean absolute source Pearson correlation over all nonredundant panel edges",
    target_coherence = paste0(
      "same source-coherence statistic recomputed on the exact target-evaluable subset ",
      "of the frozen program"
    )
  ),
  diagnostic_tolerances = list(
    absolute = ABS_TOLERANCES,
    relative = REL_TOLERANCES,
    primary_absolute = PRIMARY_ABS_TOL,
    primary_relative = PRIMARY_REL_TOL
  ),
  naive_rejection_feasibility_rule = list(
    minimum_matches_out_of_100 = MIN_MATCHES_FOR_NAIVE_FEASIBILITY,
    match_definition = paste0(
      "absolute coherence gap <= ", PRIMARY_ABS_TOL,
      " AND relative coherence gap <= ", 100 * PRIMARY_REL_TOL, "%"
    ),
    interpretation = paste0(
      "If fewer than ", MIN_MATCHES_FOR_NAIVE_FEASIBILITY,
      " of 100 unconditional same-size panels match, naive rejection sampling ",
      "is declared impractical for that target/program and the final null must use ",
      "a separately frozen source-conditioned generator. No tolerance will be widened ",
      "after seeing target preservation."
    )
  ),
  final_null_not_yet_frozen = TRUE
)

json_out <- file.path(OUT_DIR, "coherence_null_feasibility_contract_v1.json")
jsonlite::write_json(contract, json_out, pretty = TRUE, auto_unbox = TRUE, digits = NA)

md_out <- file.path(OUT_DIR, "coherence_null_feasibility_contract_v1.md")
writeLines(
  c(
    "# Paper 4 / TCBB — Source-Coherence Null Feasibility Contract v1",
    "",
    "Status: **FROZEN BEFORE FEASIBILITY AUDIT**",
    "",
    "This is a design audit, not the final specificity null.",
    "",
    "- 100 unconditional same-size random panels per target/program.",
    "- Panels draw from the frozen target-available TCGA 10k universe.",
    "- All genes of the tested program are excluded from its candidate pool.",
    "- Panel coherence is mean absolute source Pearson edge correlation.",
    "- Target-program coherence is recomputed on its exact target-evaluable source subset.",
    "- Primary practical match: absolute gap <= 0.02 AND relative gap <= 10%.",
    "- Naive rejection is called practically feasible only if >=5/100 random panels match.",
    "",
    "If the criterion fails, the final matched null will use a separately frozen",
    "source-conditioned panel generator. The tolerance will not be relaxed using",
    "SCAN-B or METABRIC preservation results."
  ),
  md_out
)

cat(strrep("=", 122), "\n", sep = "")
cat("03b COHERENCE-NULL FEASIBILITY CONTRACT: PASS\n")
cat(strrep("=", 122), "\n", sep = "")
cat("Random panels per target/program: ", N_RANDOM_PANELS_PER_TARGET_PROGRAM, "\n", sep = "")
cat("Primary coherence match: abs gap <= ", PRIMARY_ABS_TOL,
    " AND relative gap <= ", 100 * PRIMARY_REL_TOL, "%\n", sep = "")
cat("Naive rejection feasibility: >=", MIN_MATCHES_FOR_NAIVE_FEASIBILITY,
    "/100 matching panels\n", sep = "")
cat("\nNo source correlation matrix was calculated by this script.\n")
cat("No target expression/preservation statistic was calculated.\n")
cat("\nOutputs:\n")
cat("  ", json_out, "\n", sep = "")
cat("  ", md_out, "\n", sep = "")
cat(strrep("=", 122), "\n", sep = "")
