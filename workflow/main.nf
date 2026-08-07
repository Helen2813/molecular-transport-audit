#!/usr/bin/env nextflow

include { RUN_MULTIPLICITY } from './modules/multiplicity'
include { VERIFY_MULTIPLICITY } from './modules/verify_multiplicity'


params.input = (
    "${projectDir}/../reference_results/" +
    "osteosarcoma_locked/artifacts/" +
    "projectwide_primary_multiplicity__" +
    "paper4_projectwide_primary_multiplicity.csv"
)

params.reference = params.input

params.tolerances = (
    "${projectDir}/../reference_results/" +
    "osteosarcoma_locked/tolerances.yaml"
)

params.core_dir = (
    "${projectDir}/../core"
)

params.p_column = "primary_p"

params.q_column = (
    "projectwide_q_12"
)

params.support_column = (
    "projectwide_fdr_supported"
)

params.key_columns = (
    "cohort,endpoint,module_label"
)

params.alpha = 0.05


workflow {

    input_ch = Channel.fromPath(
        params.input,
        checkIfExists: true
    )

    reference_ch = Channel.fromPath(
        params.reference,
        checkIfExists: true
    )

    tolerance_ch = Channel.fromPath(
        params.tolerances,
        checkIfExists: true
    )

    core_ch = Channel.fromPath(
        params.core_dir,
        checkIfExists: true
    )

    runner_ch = Channel.fromPath(
        "${projectDir}/bin/run_multiplicity.py",
        checkIfExists: true
    )

    verifier_ch = Channel.fromPath(
        "${projectDir}/bin/verify_multiplicity.py",
        checkIfExists: true
    )


    RUN_MULTIPLICITY(
        input_ch,
        core_ch,
        runner_ch,
        params.p_column,
        params.q_column,
        params.support_column,
        params.alpha
    )


    VERIFY_MULTIPLICITY(
        RUN_MULTIPLICITY.out.table,
        reference_ch,
        tolerance_ch,
        verifier_ch,
        params.key_columns,
        params.q_column,
        params.support_column
    )


    RUN_MULTIPLICITY.out.table.view {
        value ->
        "Multiplicity output: ${value}"
    }


    RUN_MULTIPLICITY.out.manifest.view {
        value ->
        "Run manifest: ${value}"
    }


    VERIFY_MULTIPLICITY.out.verification.view {
        value ->
        "Verification report: ${value}"
    }
}
