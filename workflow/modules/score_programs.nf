process SCORE_PROGRAMS {

    tag "frozen-program-scoring"

    input:
    path expression
    path weights
    path scoring_config
    path scoring_script

    output:
    path "observed_scores.csv",
        emit: scores

    path "observed_coverage.csv",
        emit: coverage

    path "scoring_manifest.json",
        emit: manifest

    script:
    """
    python3 ${scoring_script} \
        --expression ${expression} \
        --weights ${weights} \
        --config ${scoring_config} \
        --scores observed_scores.csv \
        --coverage observed_coverage.csv \
        --manifest scoring_manifest.json
    """
}
