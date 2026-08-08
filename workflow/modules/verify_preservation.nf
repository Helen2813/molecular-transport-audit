process VERIFY_PRESERVATION {

    tag "verify-direct-preservation"

    input:
    path observed_structure
    path observed_coverage
    path expected_structure
    path expected_coverage
    path verifier_script

    output:
    path "preservation_nextflow_verification.json",
        emit: verification

    script:
    """
    python3 ${verifier_script} \
        --observed-structure ${observed_structure} \
        --observed-coverage ${observed_coverage} \
        --expected-structure ${expected_structure} \
        --expected-coverage ${expected_coverage} \
        --output preservation_nextflow_verification.json
    """
}
