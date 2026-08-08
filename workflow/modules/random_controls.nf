process RANDOM_CONTROLS {

    tag "corrected-random-panel-controls"

    input:
    path reference_array
    path external_array
    path background_genes
    path weights
    path observed_structure
    path audit_config
    path random_control_script

    output:
    path "observed_corrected_random_controls.csv",
        emit: controls

    path "random_controls_manifest.json",
        emit: manifest

    script:
    """
    python3 ${random_control_script} \
        --reference-array ${reference_array} \
        --external-array ${external_array} \
        --genes ${background_genes} \
        --weights ${weights} \
        --observed-structure ${observed_structure} \
        --config ${audit_config} \
        --output observed_corrected_random_controls.csv \
        --manifest random_controls_manifest.json
    """
}
