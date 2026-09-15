options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "05i3-freeze-comparator-interpretation-and-headtohead-v1-no-cli"

DATA_ROOT <- "D:/paper4_tcbb_data"

BASELINE_CONTRACT <- file.path(
  DATA_ROOT,
  "paper4_tcbb_postprimary_comparator_benchmark_contract_v1",
  "postprimary_comparator_benchmark_contract_v1.json"
)

OUT_DIR <- file.path(
  DATA_ROOT,
  "paper4_tcbb_comparator_interpretation_headtohead_contract_v1"
)
OUT_JSON <- file.path(
  OUT_DIR,
  "comparator_interpretation_headtohead_contract_v1.json"
)

# Comparator result locations that must still contain NO scientific comparator
# result before this interpretation/head-to-head contract is frozen.
POSSIBLE_RESULT_DIRS <- file.path(
  DATA_ROOT,
  c(
    "paper4_tcbb_postprimary_comparator_benchmark_v1",
    "paper4_tcbb_postprimary_comparator_benchmark_v2",
    "paper4_tcbb_postprimary_comparator_benchmark_v3"
  )
)

HEADTOHEAD_FRACTIONS <- c(0.00, 0.10, 0.25, 0.50, 0.75, 1.00)
HEADTOHEAD_NONZERO_REPEATS <- 100L
HEADTOHEAD_BOOTSTRAP_REPEATS <- 2000L
HEADTOHEAD_BASE_SEED <- 20260920L

cat(strrep("=", 154), "\n", sep = "")
cat("Paper 4 / TCBB - freeze comparator interpretation + mandatory known-answer head-to-head plan\n")
cat(strrep("=", 154), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n\n", sep = "")

cat("Provenance guard:\n")
cat("  Primary MTA results already observed:                 YES\n")
cat("  05i2 comparator settings already frozen:              YES\n")
cat("  NetRep/WGCNA scientific comparator results observed:   NO\n")
cat("  Interpretation rules frozen before comparator results: YES\n")
cat("  Head-to-head corruption benchmark required regardless\n")
cat("    of baseline comparator outcome:                      YES\n")
cat(strrep("=", 154), "\n\n", sep = "")

if (!file.exists(BASELINE_CONTRACT)) {
  stop("Missing 05i2 baseline comparator contract: ", BASELINE_CONTRACT)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required.")
}

baseline <- jsonlite::read_json(
  BASELINE_CONTRACT,
  simplifyVector = TRUE
)
if (as.character(baseline$status) !=
    "FROZEN_POST_PRIMARY_BEFORE_FIRST_COMPARATOR_STATISTIC") {
  stop("05i2 comparator contract has unexpected status.")
}

# Allow empty directories/checkpoint folders from failed input parsing, but do
# not allow any actual scientific comparator artifact to exist before this freeze.
scientific_patterns <- c(
  "netrep_modulePreservation.*\\.rds$",
  "netrep_modulePreservation.*\\.tsv$",
  "wgcna_modulePreservation.*\\.rds$",
  "wgcna_modulePreservation.*\\.tsv$",
  "integrated_mta_netrep_wgcna.*\\.tsv$",
  "netrep_all_targets.*\\.tsv$",
  "wgcna_modulePreservation_all_targets.*\\.tsv$",
  "mta_netrep_wgcna_integrated_all_targets.*\\.tsv$"
)

found_scientific <- character(0)
for (d in POSSIBLE_RESULT_DIRS) {
  if (!dir.exists(d)) next
  ff <- list.files(d, recursive = TRUE, full.names = TRUE)
  ff <- ff[file.info(ff)$isdir %in% FALSE]
  if (length(ff) == 0) next
  bn <- basename(ff)
  keep <- rep(FALSE, length(ff))
  for (pat in scientific_patterns) {
    keep <- keep | grepl(pat, bn, ignore.case = TRUE)
  }
  found_scientific <- c(found_scientific, ff[keep])
}
if (length(found_scientific) > 0) {
  stop(
    "Scientific comparator artifacts already exist before 05i3 freeze:\n",
    paste(found_scientific, collapse = "\n"),
    "\nDo not continue until provenance is reviewed."
  )
}

contract <- list(
  contract_id =
    "paper4-tcbb-comparator-interpretation-headtohead-v1",
  script_version = SCRIPT_VERSION,
  status =
    "FROZEN_BEFORE_FIRST_NETREP_OR_WGCNA_COMPARATOR_RESULT",

  provenance = list(
    primary_mta_results_already_observed = TRUE,
    baseline_05i2_settings_frozen = TRUE,
    comparator_results_observed_before_this_freeze = FALSE,
    rationale = paste(
      "This contract closes interpretive degrees of freedom before",
      "the first NetRep or WGCNA comparator result is inspected."
    )
  ),

  netrep_interpretation = list(
    all_seven_statistics_retained = TRUE,
    no_single_omnibus_preserved_call = TRUE,
    reason = paste(
      "NetRep defines seven distinct preservation statistics with",
      "different interpretations and caveats; the package does not",
      "define a canonical all-seven binary module-preserved rule.",
      "No post-result composite will be invented."
    ),
    statistic_specific_testing = list(
      p_values =
        "NetRep one-sided permutation p-values from the frozen 05i2 run",
      multiplicity =
        paste(
          "BH separately within target and statistic across all 12",
          "frozen source modules, exactly as frozen in 05i2."
        )
    ),
    closest_predeclared_analogs = list(
      structural =
        "NetRep cor.cor versus MTA rho_edge",
      loading_order =
        "NetRep cor.contrib versus MTA rho_load",
      loading_order_guard =
        paste(
          "cor.contrib is interpreted as loading-order support only",
          "when avg.contrib is also BH q < 0.05, following the",
          "package documentation caution that cor.contrib is not",
          "meaningful when average contribution is not supported."
        )
    ),
    descriptive_axis_support = list(
      structural =
        "cor.cor BH q < 0.05",
      loading =
        "avg.contrib BH q < 0.05 AND cor.contrib BH q < 0.05"
    ),
    remaining_statistics =
      paste(
        "coherence, avg.contrib, avg.cor, avg.weight, cor.degree",
        "remain fully reported; none may be discarded because its",
        "result is inconvenient."
      )
  ),

  wgcna_interpretation = list(
    headline_continuous = c(
      "Zsummary.pres",
      "medianRank.pres",
      "Zdensity.pres",
      "Zconnectivity.pres"
    ),
    Zsummary_categories = list(
      no_evidence = "Zsummary.pres < 2",
      moderate = "2 <= Zsummary.pres < 10",
      strong = "Zsummary.pres >= 10"
    ),
    medianRank =
      paste(
        "Continuous relative-preservation rank only; lower is stronger.",
        "No new absolute medianRank cut-off is introduced."
      ),
    no_cross_method_pvalue_equating = TRUE,
    reason =
      paste(
        "MTA, NetRep, and WGCNA use different nulls and multiplicity",
        "families. Cross-method inference therefore emphasizes",
        "continuous statistics and explicitly method-specific calls."
      )
  ),

  cross_method_descriptive_comparisons = list(
    primary_assessable_pairs_only_for_headline = TRUE,
    pairs = list(
      mta_edge_vs_netrep_corcor =
        "Spearman across assessable programs within each target",
      mta_loading_vs_netrep_corcontrib =
        "Spearman across assessable programs within each target",
      mta_edge_vs_wgcna_zsummary =
        "Spearman across assessable programs within each target",
      mta_edge_vs_negative_wgcna_medianrank =
        "Spearman across assessable programs within each target"
    ),
    inferential_p_values_for_cross_method_rank_correlations = FALSE,
    reason =
      "Only 12 SCAN-B and 10 METABRIC primary-assessable programs."
  ),

  frozen_interpretation_matrix = list(
    broad_saturation_in_conventional_comparators = paste(
      "If conventional comparators also support preservation for nearly",
      "all primary-assessable programs and do not materially separate",
      "the same programs, conclude that MTA has NOT demonstrated",
      "superior binary discrimination on these real cohorts.",
      "This may support convergent preservation evidence, but not",
      "incremental validity."
    ),
    axis_specific_agreement = paste(
      "If an MTA axis-specific distinction is mirrored by the closest",
      "NetRep/WGCNA axis-specific statistic, interpret as convergent",
      "axis validity, not superiority."
    ),
    mta_only_axis_specific_distinction = paste(
      "If MTA separates an axis while conventional comparators remain",
      "broadly preserved, report additional decomposition observed in",
      "that specific case. Do NOT call it validated incremental",
      "discrimination unless the frozen known-answer head-to-head",
      "benchmark also supports greater corruption sensitivity."
    ),
    comparator_only_failure = paste(
      "If established comparators show clear preservation failure where",
      "MTA remains Strong, treat this as a challenge to MTA",
      "interpretation and discuss it explicitly."
    ),
    mixed_unstructured_discordance = paste(
      "If discrepancies occur in both directions without a coherent",
      "pattern, conclude that incremental validity is not demonstrated;",
      "do not relabel discordance as additional resolution."
    ),
    continuous_rank_concordance = paste(
      "Positive cross-program rank concordance is evidence of convergent",
      "validity only. It is not evidence that MTA is superior."
    )
  ),

  mandatory_headtohead_known_answer = list(
    run_regardless_of_baseline_comparator_outcome = TRUE,
    scientific_question =
      paste(
        "How sensitively do continuous preservation statistics respond",
        "when the known source-to-target gene correspondence is",
        "progressively corrupted?"
      ),
    same_target_program_scope =
      "the same 22 primary-assessable target/program pairs used in 05d",
    exact_corruption_fractions = HEADTOHEAD_FRACTIONS,
    exact_nonzero_repeats_per_fraction = HEADTOHEAD_NONZERO_REPEATS,
    corruption_realizations = paste(
      "Reuse the exact 05d corrupted mappings if they were serialized.",
      "If not serialized, replay the frozen 05d generator/seeds exactly",
      "and require reproduction of the saved 05d MTA edge-rho results",
      "within numerical tolerance before accepting comparator results."
    ),
    mta_statistic =
      "existing 05d rho_edge on the frozen corrected edge subset",
    netrep_primary_statistic = list(
      name = "cor.cor",
      mode =
        "observed continuous statistic only; NetRep nPerm=0 in ladder",
      reason =
        "closest predeclared structural analogue to MTA rho_edge"
    ),
    netrep_secondary_statistics =
      "retain all other observed NetRep statistics, but do not select among them post hoc",
    wgcna_point_statistics = list(
      mode =
        "WGCNA modulePreservation nPermutations=0 for the corruption ladder",
      retain =
        paste(
          "all observed preservation statistics returned by WGCNA;",
          "pre-identify cor.cor and cor.kIM as structural/connectivity",
          "descriptors if present in the installed 1.74 output"
        ),
      reason =
        paste(
          "The ladder compares continuous response to known corruption;",
          "permutation significance at every corrupted replicate would",
          "mix statistical power with the operating-response question."
        )
    ),
    no_binary_reclassification_in_ladder = TRUE,
    within_method_normalization = list(
      formula =
        "(S(f) - median[S(1.0)]) / (S(0) - median[S(1.0)])",
      orientation =
        "statistics are oriented so larger means stronger preservation",
      use =
        paste(
          "compare curve shape across methods without equating their",
          "raw numerical scales"
        )
    ),
    operating_summaries = c(
      "median normalized response at each corruption fraction",
      "AUC of normalized response versus corruption fraction",
      "linear-interpolated f50 where normalized response first reaches 0.5",
      "Spearman correlation of corruption fraction with median raw statistic"
    ),
    uncertainty = list(
      paired_bootstrap_replicates = HEADTOHEAD_BOOTSTRAP_REPEATS,
      base_seed = HEADTOHEAD_BASE_SEED,
      pairing_rule =
        paste(
          "Within each fraction, resample the same corruption replicate",
          "indices across methods so AUC/f50 contrasts are paired."
        ),
      interval = "95% percentile interval"
    ),
    falsifiable_interpretation = list(
      mta_more_responsive = paste(
        "If MTA has lower normalized AUC / lower f50 than a comparator",
        "with a paired contrast interval excluding zero in the expected",
        "direction, this supports greater sensitivity to mapping",
        "corruption for this known-answer task only."
      ),
      comparator_more_responsive = paste(
        "If a comparator has lower normalized AUC / lower f50 than MTA",
        "with the paired interval excluding zero, MTA's claim of",
        "greater corruption discrimination is weakened/challenged."
      ),
      similar_response = paste(
        "If paired contrasts overlap zero materially, no incremental",
        "operating advantage is demonstrated."
      ),
      flatter_curve_note = paste(
        "In this known-answer experiment, a flatter curve is interpreted",
        "as lower sensitivity to deliberately wrong correspondence, not",
        "automatically as desirable robustness. No universal method",
        "superiority is inferred from this benchmark alone."
      )
    )
  ),

  manuscript_reporting = list(
    primary_table_rule =
      paste(
        "Lead with continuous effect sizes and uncertainty; categorical",
        "MTA class appears only as a final summary column."
      ),
    netrep_binary_omnibus_column = FALSE,
    technical_ceiling_rule =
      "compare only matched-n biological transport and technical repeatability",
    degradation_rule =
      "show continuous corruption curves; do not invent a post-hoc Strong-to-not-Strong breakpoint"
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

cat("\n", strrep("=", 154), "\n", sep = "")
cat("05i3 COMPARATOR INTERPRETATION + HEAD-TO-HEAD CONTRACT: PASS\n")
cat(strrep("=", 154), "\n", sep = "")
cat("NetRep omnibus binary preserved call:              NO\n")
cat("NetRep structural axis:                            cor.cor\n")
cat("NetRep loading axis guard:                         avg.contrib + cor.contrib\n")
cat("WGCNA binary threshold:                            conventional Zsummary only\n")
cat("Cross-method p-value equivalence assumed:           NO\n")
cat("Mixed discordance => incremental validity:          NOT DEMONSTRATED\n")
cat("Head-to-head corruption ladder required:            YES, regardless of 05i outcome\n")
cat("Head-to-head fractions:                             ",
    paste(HEADTOHEAD_FRACTIONS, collapse = ", "), "\n", sep = "")
cat("Head-to-head nonzero repeats/fraction:              ",
    HEADTOHEAD_NONZERO_REPEATS, "\n", sep = "")
cat("Comparator ladder uses point statistics:            YES\n")
cat("Primary MTA classification changed:                 NO\n")
cat("\nOutput: ", OUT_JSON, "\n", sep = "")
cat(strrep("=", 154), "\n", sep = "")
