process CHECK_BROKEN_LINKS {
    label 'process_long'

    conda "conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' :
        'biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0' }"

    input:
    tuple val(accno), val(genome_fna), val(genome_gff)

    output:
    tuple val(accno), val(genome_fna), val(genome_gff), emit: valid_genomes_ch
    path "versions.yml"                               , emit: versions

    script:
    """
    echo "Checking URL: ${genome_fna}"
    # Use curl to check if the URL returns 404
    if curl -Is "${genome_fna}" | grep -q "404 Not Found"; then
        echo "Broken link: ${genome_fna}"
        exit 0  # Exit successfully but don't emit anything
    else
        # Correctly echo the variables
        echo "{accno:'${accno}', genome_fna:'${genome_fna}', genome_gff:'${genome_gff}'}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """
}
