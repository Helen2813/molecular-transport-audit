include {
    PRESERVE_PROGRAMS
} from './modules/preserve_programs'

include {
    PRESERVE_INFERENCE
} from './modules/preserve_inference'

include {
    PRESERVE_RELIABILITY
} from './modules/preserve_reliability'

include {
    RANDOM_CONTROLS
} from './modules/random_controls'

include {
    ASSEMBLE_AUDIT_SUMMARY
} from './modules/assemble_audit_summary'

include {
    VERIFY_UNIFIED_AUDIT
} from './modules/verify_unified_audit'


params.reference = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/DOG2_reference_module_expression.csv"

params.external = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/GSE239948_external_module_expression.csv"

params.weights = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/primary_canine_program_weights.csv"

params.background_reference = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/random_panel_DOG2_prestandardization.npy"

params.background_external = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/random_panel_GSE239948_prestandardization.npy"

params.background_genes = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/random_panel_background_genes.csv"

params.expected_direct = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/expected_direct_preservation.csv"

params.expected_coverage = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/expected_preservation_coverage.csv"

params.expected_reliability = "${launchDir}/reference_results/osteosarcoma_locked/preservation_fixture/expected_split_half_reliability.csv"

params.audit_config = "${launchDir}/configs/gse239948_unified_audit.json"

params.direct_script = "${launchDir}/workflow/scripts/preserve_programs.py"

params.inference_script = "${launchDir}/workflow/scripts/preserve_inference.py"

params.reliability_script = "${launchDir}/workflow/scripts/preserve_reliability.py"

params.random_control_script = "${launchDir}/workflow/scripts/run_random_controls.py"

params.assembler_script = "${launchDir}/workflow/scripts/assemble_audit_summary.py"

params.verifier_script = "${launchDir}/workflow/scripts/verify_unified_audit.py"


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

    background_reference = file(
        params.background_reference,
        checkIfExists: true
    )

    background_external = file(
        params.background_external,
        checkIfExists: true
    )

    background_genes = file(
        params.background_genes,
        checkIfExists: true
    )

    expected_direct = file(
        params.expected_direct,
        checkIfExists: true
    )

    expected_coverage = file(
        params.expected_coverage,
        checkIfExists: true
    )

    expected_reliability = file(
        params.expected_reliability,
        checkIfExists: true
    )

    audit_config = file(
        params.audit_config,
        checkIfExists: true
    )

    direct_script = file(
        params.direct_script,
        checkIfExists: true
    )

    inference_script = file(
        params.inference_script,
        checkIfExists: true
    )

    reliability_script = file(
        params.reliability_script,
        checkIfExists: true
    )

    random_control_script = file(
        params.random_control_script,
        checkIfExists: true
    )

    assembler_script = file(
        params.assembler_script,
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
        audit_config,
        direct_script
    )

    PRESERVE_INFERENCE(
        reference_expression,
        external_expression,
        weights,
        audit_config,
        inference_script
    )

    PRESERVE_RELIABILITY(
        external_expression,
        weights,
        audit_config,
        reliability_script
    )

    RANDOM_CONTROLS(
        background_reference,
        background_external,
        background_genes,
        weights,
        PRESERVE_PROGRAMS.out.structure,
        audit_config,
        random_control_script
    )

    ASSEMBLE_AUDIT_SUMMARY(
        PRESERVE_PROGRAMS.out.structure,
        PRESERVE_PROGRAMS.out.coverage,
        PRESERVE_INFERENCE.out.inference,
        PRESERVE_RELIABILITY.out.reliability,
        PRESERVE_RELIABILITY.out.loo,
        RANDOM_CONTROLS.out.controls,
        audit_config,
        assembler_script
    )

    VERIFY_UNIFIED_AUDIT(
        ASSEMBLE_AUDIT_SUMMARY.out.module_summary,
        ASSEMBLE_AUDIT_SUMMARY.out.summary,
        expected_direct,
        expected_coverage,
        expected_reliability,
        audit_config,
        verifier_script
    )

    publish:

    direct_structure = PRESERVE_PROGRAMS.out.structure
    direct_coverage = PRESERVE_PROGRAMS.out.coverage
    direct_manifest = PRESERVE_PROGRAMS.out.manifest

    permutation_inference = PRESERVE_INFERENCE.out.inference
    permutation_manifest = PRESERVE_INFERENCE.out.manifest

    reliability = PRESERVE_RELIABILITY.out.reliability
    gene_loo = PRESERVE_RELIABILITY.out.loo
    reliability_manifest = PRESERVE_RELIABILITY.out.manifest

    random_controls = RANDOM_CONTROLS.out.controls
    random_controls_manifest = RANDOM_CONTROLS.out.manifest

    module_summary = ASSEMBLE_AUDIT_SUMMARY.out.module_summary
    summary = ASSEMBLE_AUDIT_SUMMARY.out.summary

    verification = VERIFY_UNIFIED_AUDIT.out.verification
}


output {

    direct_structure {
        path "direct"
        mode "copy"
    }

    direct_coverage {
        path "direct"
        mode "copy"
    }

    direct_manifest {
        path "direct"
        mode "copy"
    }

    permutation_inference {
        path "inference"
        mode "copy"
    }

    permutation_manifest {
        path "inference"
        mode "copy"
    }

    reliability {
        path "reliability"
        mode "copy"
    }

    gene_loo {
        path "reliability"
        mode "copy"
    }

    reliability_manifest {
        path "reliability"
        mode "copy"
    }

    random_controls {
        path "random_controls"
        mode "copy"
    }

    random_controls_manifest {
        path "random_controls"
        mode "copy"
    }

    module_summary {
        path "summary"
        mode "copy"
    }

    summary {
        path "summary"
        mode "copy"
    }

    verification {
        path "verification"
        mode "copy"
    }
}
