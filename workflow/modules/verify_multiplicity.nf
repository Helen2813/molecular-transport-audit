process VERIFY_MULTIPLICITY {

    tag "locked-reference"

    input:
    path observed_csv
    path reference_csv
    path tolerance_file
    path verifier
    val key_columns
    val q_column
    val support_column

    output:
    path "multiplicity_verification.json", emit: verification

    script:
    """
    python3 "${verifier}" \
        --observed "${observed_csv}" \
        --reference "${reference_csv}" \
        --tolerances "${tolerance_file}" \
        --output multiplicity_verification.json \
        --keys "${key_columns}" \
        --q-column "${q_column}" \
        --support-column "${support_column}"
    """
}
