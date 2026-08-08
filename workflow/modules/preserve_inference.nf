process PRESERVE_INFERENCE {

    tag "preservation-permutation-inference"

    input:
    path reference_expression
    path external_expression
    path weights
    path audit_config
    path inference_script

    output:
    path "observed_permutation_preservation.csv",
        emit: inference

    path "permutation_manifest.json",
        emit: manifest

    script:
    """
    python3 ${inference_script} \
        --reference ${reference_expression} \
        --external ${external_expression} \
        --weights ${weights} \
        --config ${audit_config} \
        --output observed_permutation_preservation.csv \
        --manifest permutation_manifest.json
    """
}
