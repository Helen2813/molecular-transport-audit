process VERIFY_SCORING {

    tag "verify-frozen-scoring"

    input:
    path observed_scores
    path observed_coverage
    path expected_scores
    path expected_coverage
    path verifier_script

    output:
    path "scoring_nextflow_verification.json",
        emit: verification

    script:
    """
    python3 ${verifier_script} \
        --observed-scores ${observed_scores} \
        --observed-coverage ${observed_coverage} \
        --expected-scores ${expected_scores} \
        --expected-coverage ${expected_coverage} \
        --output scoring_nextflow_verification.json
    """
}
