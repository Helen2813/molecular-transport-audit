include {
    SCORE_PROGRAMS
} from './modules/score_programs'

include {
    VERIFY_SCORING
} from './modules/verify_scoring'


params.expression = "${launchDir}/reference_results/osteosarcoma_locked/scoring_fixture/TARGET_OS_primary_strict_expression.csv"

params.weights = "${launchDir}/reference_results/osteosarcoma_locked/scoring_fixture/primary_strict_gene_weights.csv"

params.expected_scores = "${launchDir}/reference_results/osteosarcoma_locked/scoring_fixture/TARGET_OS_expected_primary_scores.csv"

params.expected_coverage = "${launchDir}/reference_results/osteosarcoma_locked/scoring_fixture/TARGET_OS_expected_strict_coverage.csv"

params.scoring_config = "${launchDir}/configs/scoring_target_regression.json"

params.scoring_script = "${launchDir}/workflow/scripts/score_programs.py"

params.verifier_script = "${launchDir}/workflow/scripts/verify_scoring.py"


workflow {

    main:

    expression = file(
        params.expression,
        checkIfExists: true
    )

    weights = file(
        params.weights,
        checkIfExists: true
    )

    scoring_config = file(
        params.scoring_config,
        checkIfExists: true
    )

    scoring_script = file(
        params.scoring_script,
        checkIfExists: true
    )

    expected_scores = file(
        params.expected_scores,
        checkIfExists: true
    )

    expected_coverage = file(
        params.expected_coverage,
        checkIfExists: true
    )

    verifier_script = file(
        params.verifier_script,
        checkIfExists: true
    )

    SCORE_PROGRAMS(
        expression,
        weights,
        scoring_config,
        scoring_script
    )

    VERIFY_SCORING(
        SCORE_PROGRAMS.out.scores,
        SCORE_PROGRAMS.out.coverage,
        expected_scores,
        expected_coverage,
        verifier_script
    )

    publish:

    observed_scores = SCORE_PROGRAMS.out.scores

    observed_coverage = SCORE_PROGRAMS.out.coverage

    scoring_manifest = SCORE_PROGRAMS.out.manifest

    verification = VERIFY_SCORING.out.verification
}


output {

    observed_scores {
        path "score_programs"
        mode "copy"
    }

    observed_coverage {
        path "score_programs"
        mode "copy"
    }

    scoring_manifest {
        path "score_programs"
        mode "copy"
    }

    verification {
        path "verification"
        mode "copy"
    }
}
