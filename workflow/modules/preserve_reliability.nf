process PRESERVE_RELIABILITY {

    tag "frozen-program-reliability"

    input:
    path external_expression
    path weights
    path reliability_config
    path reliability_script

    output:
    path "observed_reliability.csv",
        emit: reliability

    path "observed_gene_leave_one_out.csv",
        emit: loo

    path "reliability_manifest.json",
        emit: manifest

    script:
    """
    python3 ${reliability_script} \
        --external ${external_expression} \
        --weights ${weights} \
        --config ${reliability_config} \
        --reliability observed_reliability.csv \
        --loo observed_gene_leave_one_out.csv \
        --manifest reliability_manifest.json
    """
}
