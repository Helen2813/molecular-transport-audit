process ASSEMBLE_AUDIT_SUMMARY {

    tag "assemble-audit-summary"

    input:
    path direct_structure
    path coverage
    path inference
    path reliability
    path gene_loo
    path random_controls
    path audit_config
    path assembler_script

    output:
    path "module_summary.csv",
        emit: module_summary

    path "summary.json",
        emit: summary

    script:
    """
    python3 ${assembler_script} \
        --direct ${direct_structure} \
        --coverage ${coverage} \
        --inference ${inference} \
        --reliability ${reliability} \
        --loo ${gene_loo} \
        --random-controls ${random_controls} \
        --config ${audit_config} \
        --module-summary module_summary.csv \
        --summary summary.json
    """
}
