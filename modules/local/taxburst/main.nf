process TAXBURST {
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxburst:0.3.2--pyhdfd78af_0' :
        'biocontainers/taxburst:0.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(kraken_report)

    output:
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml", emit: versions

    script:
    """
    taxburst -F krona ${kraken_report} -o ${meta.id}_taxburst.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxburst: \$(taxburst --version 2>&1 | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
