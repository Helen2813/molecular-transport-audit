include {
    PRESERVE_RELIABILITY
} from './modules/preserve_reliability'

include {
    VERIFY_RELIABILITY
} from './modules/verify_reliability'


params.external = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/GSE239948_external_module_expression.csv"

params.weights = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/primary_canine_program_weights.csv"

params.expected_reliability = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/expected_split_half_reliability.csv"

params.expected_loo = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/expected_gene_leave_one_out.csv"

params.reliability_config = "${launchDir}/configs/preservation_reliability_gse239948_regression.json"

params.reliability_script = "${launchDir}/workflow/scripts/preserve_reliability.py"

params.verifier_script = "${launchDir}/workflow/scripts/verify_reliability.py"


workflow {

    main:

    external_expression = file(
        params.external,
        checkIfExists: true
    )

    weights = file(
        params.weights,
        checkIfExists: true
    )

    expected_reliability = file(
        params.expected_reliability,
        checkIfExists: true
    )

    expected_loo = file(
        params.expected_loo,
        checkIfExists: true
    )

    reliability_config = file(
        params.reliability_config,
        checkIfExists: true
    )

    reliability_script = file(
        params.reliability_script,
        checkIfExists: true
    )

    verifier_script = file(
        params.verifier_script,
        checkIfExists: true
    )

    PRESERVE_RELIABILITY(
        external_expression,
        weights,
        reliability_config,
        reliability_script
    )

    VERIFY_RELIABILITY(
        PRESERVE_RELIABILITY.out.reliability,
        PRESERVE_RELIABILITY.out.loo,
        expected_reliability,
        expected_loo,
        verifier_script
    )

    publish:

    observed_reliability = PRESERVE_RELIABILITY.out.reliability
    observed_loo = PRESERVE_RELIABILITY.out.loo
    reliability_manifest = PRESERVE_RELIABILITY.out.manifest
    verification = VERIFY_RELIABILITY.out.verification
}


output {

    observed_reliability {
        path "preserve_reliability"
        mode "copy"
    }

    observed_loo {
        path "preserve_reliability"
        mode "copy"
    }

    reliability_manifest {
        path "preserve_reliability"
        mode "copy"
    }

    verification {
        path "verification"
        mode "copy"
    }
}
