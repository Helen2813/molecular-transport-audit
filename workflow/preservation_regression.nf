include {
    PRESERVE_PROGRAMS
} from './modules/preserve_programs'

include {
    VERIFY_PRESERVATION
} from './modules/verify_preservation'


params.reference = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/DOG2_reference_module_expression.csv"

params.external = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/GSE239948_external_module_expression.csv"

params.weights = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/primary_canine_program_weights.csv"

params.expected_structure = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/expected_direct_preservation.csv"

params.expected_coverage = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/expected_preservation_coverage.csv"

params.preservation_config = "${launchDir}/configs/preservation_gse239948_regression.json"

params.preservation_script = "${launchDir}/workflow/scripts/preserve_programs.py"

params.verifier_script = "${launchDir}/workflow/scripts/verify_preservation.py"


workflow {

    main:

    reference_expression = file(
        params.reference,
        checkIfExists: true
    )

    external_expression = file(
        params.external,
        checkIfExists: true
    )

    weights = file(
        params.weights,
        checkIfExists: true
    )

    expected_structure = file(
        params.expected_structure,
        checkIfExists: true
    )

    expected_coverage = file(
        params.expected_coverage,
        checkIfExists: true
    )

    preservation_config = file(
        params.preservation_config,
        checkIfExists: true
    )

    preservation_script = file(
        params.preservation_script,
        checkIfExists: true
    )

    verifier_script = file(
        params.verifier_script,
        checkIfExists: true
    )

    PRESERVE_PROGRAMS(
        reference_expression,
        external_expression,
        weights,
        preservation_config,
        preservation_script
    )

    VERIFY_PRESERVATION(
        PRESERVE_PROGRAMS.out.structure,
        PRESERVE_PROGRAMS.out.coverage,
        expected_structure,
        expected_coverage,
        verifier_script
    )

    publish:

    observed_structure = PRESERVE_PROGRAMS.out.structure
    observed_coverage = PRESERVE_PROGRAMS.out.coverage
    preservation_manifest = PRESERVE_PROGRAMS.out.manifest
    verification = VERIFY_PRESERVATION.out.verification

}


output {

    observed_structure {
        path "preserve_programs"
        mode "copy"
    }

    observed_coverage {
        path "preserve_programs"
        mode "copy"
    }

    preservation_manifest {
        path "preserve_programs"
        mode "copy"
    }

    verification {
        path "verification"
        mode "copy"
    }
}
