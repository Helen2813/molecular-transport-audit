process VERIFY_UNIFIED_AUDIT {

    tag "verify-unified-audit"

    input:
    path module_summary
    path summary
    path expected_direct
    path expected_coverage
    path expected_reliability
    path audit_config
    path verifier_script

    output:
    path "unified_audit_verification.json",
        emit: verification

    script:
    """
    python3 ${verifier_script} \
        --module-summary ${module_summary} \
        --summary ${summary} \
        --expected-direct ${expected_direct} \
        --expected-coverage ${expected_coverage} \
        --expected-reliability ${expected_reliability} \
        --config ${audit_config} \
        --output unified_audit_verification.json
    """
}
