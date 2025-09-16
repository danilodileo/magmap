process KRAKEN2_DOWNLOAD_DB {
    tag "$db_name"
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data':
        'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e' }"
    
    input:
    val db_name
    val db_url
    
    output:
    path "${db_name}/" , emit: db_dir
    path "versions.yml", emit: versions
    
    script:
    """
    # Download the database
    wget -O ${db_name}.tar.gz ${db_url}
    
    # Create directory and extract
    mkdir -p ${db_name}
    tar -xzf ${db_name}.tar.gz -C ${db_name}
    
    # Clean up archive
    rm ${db_name}.tar.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>&1 | head -n1 | cut -d' ' -f3)
        tar: \$(tar --version 2>&1 | head -n1 | sed 's/^.*(GNU tar) //; s/ .*\$//')
    END_VERSIONS
    """
}