process RUN_MULTIPLICITY {

    tag "${input_csv.name}"

    input:
    path input_csv
    path core_dir
    path runner
    val p_column
    val q_column
    val support_column
    val alpha

    output:
    path "multiplicity_results.csv", emit: table
    path "multiplicity_run_manifest.json", emit: manifest

    script:
    """
    export PYTHONPATH="${core_dir}"

    python3 "${runner}" \
        --input "${input_csv}" \
        --output multiplicity_results.csv \
        --manifest multiplicity_run_manifest.json \
        --p-column "${p_column}" \
        --q-column "${q_column}" \
        --support-column "${support_column}" \
        --alpha "${alpha}"
    """
}
