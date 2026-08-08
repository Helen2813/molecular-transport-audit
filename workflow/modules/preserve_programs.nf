process PRESERVE_PROGRAMS {

    tag "direct-molecular-preservation"

    input:
    path reference_expression
    path external_expression
    path weights
    path preservation_config
    path preservation_script

    output:
    path "observed_direct_preservation.csv",
        emit: structure

    path "observed_preservation_coverage.csv",
        emit: coverage

    path "preservation_manifest.json",
        emit: manifest

    script:
    """
    python3 ${preservation_script} \
        --reference ${reference_expression} \
        --external ${external_expression} \
        --weights ${weights} \
        --config ${preservation_config} \
        --structure observed_direct_preservation.csv \
        --coverage observed_preservation_coverage.csv \
        --manifest preservation_manifest.json
    """
}
