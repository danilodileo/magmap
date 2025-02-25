process CHECK_DUPLICATES {
    label 'process_low'
    tag "${meta.id}"

    conda "conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'biocontainers/ubuntu:24.04' }"

    input:
    tuple val(meta), path(fna)

    output:
    stdout emit: result
    path "versions.yml"              , emit: versions

    script:
    prefix = task.ext.prefix ?: meta.id

    """
    zgrep -H '>' *.fna.gz | sed 's/^[^:]*://' | sort | uniq -d > temp_dupes.txt
    zgrep -l -F -f temp_dupes.txt *.fna.gz | sort -u > duplicate_contig_names.txt || touch duplicate_contig_names.txt
    rm temp_dupes.txt
    cat duplicate_contig_names.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zgrep: \$( zgrep --version | sed 's/.*/1.5/' | head -n 1 )
    END_VERSIONS
    """
}
