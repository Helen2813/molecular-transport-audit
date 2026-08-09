export type RunStatus =
  | "queued"
  | "running"
  | "completed"
  | "failed";

export type RunStage =
  | "queued"
  | "starting"
  | "direct_preservation"
  | "permutation_inference"
  | "reliability"
  | "random_controls"
  | "summary_assembly"
  | "verification"
  | "completed"
  | "failed";

export type ClassificationShort =
  | "strong"
  | "partial"
  | "limited"
  | "no_clear";

export interface RunDetail {
  run_id: string;
  analysis_id: string;

  status: RunStatus;
  stage: RunStage;

  validated: boolean;
  validation_status: "PASS" | null;

  pipeline_version: string | null;
  schema_version: string | null;

  created_at_utc: string;
  started_at_utc: string | null;
  finished_at_utc: string | null;
  generated_at_utc: string | null;

  reference_cohort: string | null;
  external_cohort: string | null;
  outcome_loaded: boolean | null;

  module_count: number | null;
  error: string | null;
}

export interface ModuleResult {
  module: string;

  coverage: {
    n_frozen_genes: number | null;
    n_common_genes: number | null;
    fraction: number | null;
  };

  direct_preservation: {
    edge_spearman: number | null;
    loading_spearman: number | null;
    external_pc1_variance_explained:
      | number
      | null;
    pc1_orientation_correlation:
      | number
      | null;
  };

  inference: {
    edge_permutation_p: number | null;
    edge_q_bh_8: number | null;
    edge_extreme_count: number | null;

    loading_permutation_p: number | null;
    loading_q_bh_8: number | null;
    loading_extreme_count: number | null;
  };

  reliability: {
    split_half_median: number | null;
    split_half_q05: number | null;
    split_half_q95: number | null;
    valid_repeats: number | null;

    minimum_gene_loo_correlation:
      | number
      | null;

    median_gene_loo_correlation:
      | number
      | null;
  };

  random_control: {
    matching_mode: string;
    n_random_panels: number | null;
    null_median: number | null;
    null_q05: number | null;
    null_q95: number | null;
    empirical_p: number | null;
    fallback_draws: number | null;
  };

  classification: {
    short: ClassificationShort;
    full: string;
  };
}

export interface AuditSummary {
  schema_version: "0.1.0";
  assembler_version: string;
  generated_at_utc: string;

  run: {
    analysis_id: string;
    status: "completed";
    reference_cohort: string;
    external_cohort: string;
    outcome_loaded: boolean;
  };

  parameters: {
    permutations: number;
    split_half_repeats: number;
    random_panels: number;
    random_panel_matching_mode: string;
  };

  diagnostics: {
    gene_leave_one_out_rows: number;
    module_count: number;
  };

  modules: ModuleResult[];
}
