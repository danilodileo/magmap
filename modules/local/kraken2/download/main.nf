process KRAKEN2_DOWNLOAD_DB {
    tag "$db_name"
    label 'process_low'
    storeDir "${params.kraken2_db_store_dir ?: params.outdir}/kraken2_databases"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.21.3' :
        'biocontainers/wget:1.21.3' }"
    
    input:
    val db_name
    val db_url
    
    output:
    path "${db_name}/", emit: db_dir
    path "versions.yml", emit: versions
    
    script:
    """
    # Download the database
    wget -O ${db_name}.tar.gz ${db_url}
    
    # Create directory and extract
    mkdir -p ${db_name}
    tar -xzf ${db_name}.tar.gz -C ${db_name} --strip-components=1
    
    # Clean up archive
    rm ${db_name}.tar.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>&1 | head -n1 | cut -d' ' -f3)
        tar: \$(tar --version 2>&1 | head -n1 | sed 's/^.*(GNU tar) //; s/ .*\$//')
    END_VERSIONS
    """
}