options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "05i2-freeze-postprimary-comparator-benchmark-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"

SOURCE_MODULE_JSON <- file.path(
  DATA_ROOT,
  "paper4_tcbb_tcga_source_modules_v1",
  "tcga_source_modules_frozen_v1.json"
)
WGCNA_CONTRACT <- file.path(
  DATA_ROOT,
  "paper4_tcbb_wgcna_module_contract_v2",
  "tcga_wgcna_module_contract_v2.json"
)
CLASSIFICATION <- file.path(
  DATA_ROOT,
  "paper4_tcbb_final_mapping_specificity_null_v2",
  "primary_pooled_final_classification_v2.tsv"
)
EXACT_VARIATION_JSON <- file.path(
  DATA_ROOT,
  "paper4_tcbb_exact_variation_evaluability_correction_v1",
  "exact_variation_evaluability_correction_v1.json"
)
PLATFORM_MASTER <- file.path(
  DATA_ROOT,
  "paper4_tcbb_scanb_platform_sensitivity_v1",
  "scanb_platform_sensitivity_v1.json"
)

OUT_DIR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_postprimary_comparator_benchmark_contract_v1"
)
OUT_JSON <- file.path(
  OUT_DIR,
  "postprimary_comparator_benchmark_contract_v1.json"
)

RESULT_DIR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_postprimary_comparator_benchmark_v1"
)

NETREP_PERM <- 1000L
WGCNA_PERM <- 200L
BASE_SEED <- 20260919L

cat(strrep("=", 150), "\n", sep = "")
cat("Paper 4 / TCBB - freeze post-primary NetRep + WGCNA comparator benchmark\n")
cat(strrep("=", 150), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n\n", sep = "")

cat("Provenance:\n")
cat("  Primary MTA target results already observed:       YES\n")
cat("  Pre-primary comparator contract existed:           NO\n")
cat("  Comparator statistics inspected before this freeze:NO\n")
cat("  Role:                                              POST-PRIMARY BENCHMARK\n")
cat("  Primary MTA classification can change:             NO\n")
cat(strrep("=", 150), "\n\n", sep = "")

required <- c(
  SOURCE_MODULE_JSON,
  WGCNA_CONTRACT,
  CLASSIFICATION,
  EXACT_VARIATION_JSON,
  PLATFORM_MASTER
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

if (dir.exists(RESULT_DIR)) {
  existing <- list.files(
    RESULT_DIR,
    recursive = TRUE,
    full.names = TRUE
  )
  existing <- existing[file.info(existing)$isdir %in% FALSE]
  if (length(existing) > 0) {
    stop(
      "Comparator result directory already contains files BEFORE contract freeze:\n",
      paste(existing, collapse = "\n"),
      "\nDo not overwrite or continue until provenance is reviewed."
    )
  }
}

source_meta <- jsonlite::read_json(
  SOURCE_MODULE_JSON,
  simplifyVector = TRUE
)
if (as.character(source_meta$status) !=
    "FROZEN_SOURCE_MODULES_FINALIZED_FROM_FIRST_RUN") {
  stop("Unexpected frozen source-module status.")
}
if (as.integer(source_meta$wgcna$non_grey_modules) != 12L) {
  stop("Expected exactly 12 frozen source modules.")
}

wgcna_meta <- jsonlite::read_json(
  WGCNA_CONTRACT,
  simplifyVector = TRUE
)
if (as.character(wgcna_meta$status) !=
    "FROZEN_BEFORE_FIRST_MODULE_DETECTION") {
  stop("Unexpected frozen source WGCNA contract status.")
}

if (as.character(packageVersion("NetRep")) != "1.2.10") {
  stop(
    "Expected NetRep 1.2.10, found ",
    as.character(packageVersion("NetRep"))
  )
}
if (as.character(packageVersion("WGCNA")) != "1.74") {
  stop(
    "Expected WGCNA 1.74, found ",
    as.character(packageVersion("WGCNA"))
  )
}

contract <- list(
  contract_id =
    "paper4-tcbb-postprimary-comparator-benchmark-v1",
  script_version = SCRIPT_VERSION,
  status =
    "FROZEN_POST_PRIMARY_BEFORE_FIRST_COMPARATOR_STATISTIC",

  provenance = list(
    primary_mta_results_already_observed = TRUE,
    preprimary_comparator_contract_existed = FALSE,
    comparator_statistics_inspected_before_freeze = FALSE,
    interpretation =
      paste(
        "NetRep and WGCNA module-preservation analyses are fixed",
        "post-primary external comparators. They must not be described",
        "as prespecified primary analyses and cannot alter the frozen",
        "MTA primary classifications."
      )
  ),

  source_modules = list(
    source = "TCGA_BRCA",
    samples = 1082L,
    frozen_non_grey_modules = 12L,
    module_membership =
      "exact frozen 02e/02e1 WGCNA source membership",
    source_network_definition =
      "signed Pearson WGCNA, soft power 8"
  ),

  targets = list(
    SCANB_GSE96058 = list(
      biological_samples = 3273L,
      corrected_target_evaluable_genes = 9220L
    ),
    METABRIC = list(
      biological_samples = 1980L,
      corrected_target_evaluable_genes = 8485L
    )
  ),

  common_gene_rule = list(
    definition =
      paste(
        "For each target use the intersection of the frozen TCGA",
        "10,000-gene source universe with the corrected exact-variation",
        "target-evaluable universe from 04i3."
      ),
    source_order_retained = TRUE,
    alias_or_fuzzy_rescue = FALSE,
    target_duplicate_symbol_rule =
      "arithmetic mean on the already-frozen transformed scale",
    target_specific_gene_dropping_after_this_freeze = FALSE
  ),

  benchmark_scope = list(
    compute_modules =
      paste(
        "Compute all 12 frozen non-grey source modules when supported",
        "by the comparator package."
      ),
    headline_comparable_modules =
      paste(
        "Primary MTA-assessable target/program pairs only:",
        "12 SCAN-B and 10 METABRIC."
      ),
    metabric_outside_primary_scope =
      c("TCGA_M006", "TCGA_M011"),
    outside_scope_rule =
      paste(
        "If conventional comparators return values for M006/M011 in",
        "METABRIC, retain them in audit outputs but label them outside",
        "the primary comparable set; do not use them to revise MTA",
        "assessability."
      )
  ),

  expression_preparation = list(
    source =
      "frozen log2(RSEM+1) TCGA source matrix",
    scanb =
      paste(
        "published log2(FPKM+0.1) -> back-transform -> log2(FPKM+1),",
        "then duplicate symbols averaged"
      ),
    metabric =
      "as-published continuous log2 microarray, duplicate symbols averaged",
    comparator_data_matrix =
      paste(
        "Within each cohort and target-specific common-gene universe,",
        "gene-wise center and scale before comparator calculation."
      ),
    correlation = "Pearson"
  ),

  NetRep = list(
    package = "NetRep",
    version = "1.2.10",
    function_name = "NetRep::modulePreservation",
    data =
      "cohort-wise gene-z-standardized transformed expression",
    correlation =
      "Pearson gene-gene correlation matrix",
    network =
      "signed WGCNA-style adjacency ((1 + Pearson r)/2)^8",
    moduleAssignments =
      "frozen source program_id labels; GREY is background",
    modules = "all 12 frozen non-grey programs",
    backgroundLabel = "GREY",
    discovery = "TCGA source",
    test = "external target",
    selfPreservation = FALSE,
    nThreads = 1L,
    nPerm = NETREP_PERM,
    null = "overlap",
    alternative = "greater",
    simplify = FALSE,
    statistics_reported = c(
      "coherence",
      "avg.contrib",
      "cor.contrib",
      "cor.cor",
      "avg.cor",
      "avg.weight",
      "cor.degree"
    ),
    conceptual_anchor_to_mta = list(
      structural_axis = "cor.cor",
      loading_order_axis = "cor.contrib"
    ),
    multiplicity =
      paste(
        "Raw NetRep permutation p-values retained. In addition,",
        "Benjamini-Hochberg q-values are calculated separately within",
        "each target and each NetRep statistic across the 12 frozen",
        "modules. These q-values are comparator descriptors only."
      ),
    no_custom_composite_classification = TRUE,
    seed_rule =
      "set.seed(20260919 + target_code*1000 + 11), target_code SCANB=1 METABRIC=2"
  ),

  WGCNA_modulePreservation = list(
    package = "WGCNA",
    version = "1.74",
    function_name = "WGCNA::modulePreservation",
    dataIsExpr = TRUE,
    networkType = "signed",
    corFnc = "cor",
    corOptions = "use = 'p'",
    referenceNetworks = 1L,
    testNetworks = "2",
    nPermutations = WGCNA_PERM,
    includekMEallInSummary = FALSE,
    restrictSummaryForGeneralNetworks = TRUE,
    calculateQvalue = FALSE,
    maxGoldModuleSize = 5000L,
    maxModuleSize = 5000L,
    quickCor = 0L,
    ccTupletSize = 2L,
    calculateCor_kIMall = FALSE,
    calculateClusterCoeff = FALSE,
    useInterpolation = FALSE,
    checkData = TRUE,
    savePermutedStatistics = FALSE,
    loadPermutedStatistics = FALSE,
    plotInterpolation = FALSE,
    parallelCalculation = FALSE,
    primary_outputs = c(
      "Zsummary.pres",
      "medianRank.pres",
      "Zdensity.pres",
      "Zconnectivity.pres"
    ),
    conventional_Zsummary_interpretation = list(
      no_evidence = "Zsummary.pres < 2",
      moderate = "2 <= Zsummary.pres < 10",
      strong = "Zsummary.pres >= 10"
    ),
    medianRank_interpretation =
      "lower medianRank.pres indicates stronger relative preservation",
    seed_rule =
      "randomSeed = 20260919 + target_code*1000 + 21, target_code SCANB=1 METABRIC=2"
  ),

  cross_method_reporting = list(
    join_with_frozen_mta = c(
      "rho_edge",
      "rho_load",
      "primary_classification",
      "primary_assessable"
    ),
    no_reclassification = TRUE,
    no_threshold_tuning_after_comparator_results = TRUE,
    no_module_exclusion_based_on_comparator_results = TRUE
  ),

  computation = list(
    process_targets_sequentially = TRUE,
    NetRep_single_thread = TRUE,
    WGCNA_parallelCalculation = FALSE,
    checkpoint_each_target_and_method = TRUE,
    base_seed = BASE_SEED
  )
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  contract,
  OUT_JSON,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

cat("\n", strrep("=", 150), "\n", sep = "")
cat("05i2 POST-PRIMARY COMPARATOR BENCHMARK CONTRACT: PASS\n")
cat(strrep("=", 150), "\n", sep = "")
cat("Role:                                POST-PRIMARY comparator only\n")
cat("Targets:                             SCAN-B + METABRIC\n")
cat("NetRep permutations:                 ", NETREP_PERM, "\n", sep = "")
cat("WGCNA modulePreservation perms:      ", WGCNA_PERM, "\n", sep = "")
cat("NetRep statistics:                   ALL 7 package statistics\n")
cat("WGCNA headline outputs:              Zsummary + medianRank\n")
cat("Primary MTA classifications changed: NO\n")
cat("Comparator thresholds tuned later:   NO\n")
cat("\nOutput: ", OUT_JSON, "\n", sep = "")
cat(strrep("=", 150), "\n", sep = "")
