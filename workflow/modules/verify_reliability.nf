process VERIFY_RELIABILITY {

    tag "verify-frozen-program-reliability"

    input:
    path observed_reliability
    path observed_loo
    path expected_reliability
    path expected_loo
    path verifier_script

    output:
    path "reliability_nextflow_verification.json",
        emit: verification

    script:
    """
    python3 ${verifier_script} \
        --observed-reliability ${observed_reliability} \
        --observed-loo ${observed_loo} \
        --expected-reliability ${expected_reliability} \
        --expected-loo ${expected_loo} \
        --output reliability_nextflow_verification.json
    """
}
