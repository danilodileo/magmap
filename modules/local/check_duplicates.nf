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
    path "duplicates.txt"            , emit: duplicates_file, optional: true
    tuple val(meta), path(fna)       , emit: fna, optional: true
    path "versions.yml"              , emit: versions

    script:
    prefix = task.ext.prefix ?: meta.id

    """
    # Find duplicate contig names across all files
    zgrep -H '>' *.fna.gz | sed 's/^[^:]*://' | sort | uniq -d > duplicate_contig_names.txt

    # If duplicates are found, identify which files contain them
    if [ -s duplicate_contig_names.txt ]; then
        while read contig; do
            for file in *.fna.gz; do
                if zgrep -q "\$contig" \$file; then
                    echo "\$file" >> duplicates.txt
                fi
            done
        done < duplicate_contig_names.txt

        # Remove duplicates and sort the file list
        sort -u duplicates.txt -o duplicates.txt
    fi

        # Remove duplicates.txt if it is empty
    if [ ! -s duplicates.txt ]; then
        rm -f duplicates.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zgrep: \$( zgrep --version | sed 's/.*/1.5/' | head -n 1 )
    END_VERSIONS
    """
}
