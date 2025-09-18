process FORMAT_KRONA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/07c6b82f970ba651b6273df1378adb205050c131c40bac7f80d9354e088e0865/data':
        'community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186' }"

    input:
    tuple val(meta), path(krona_file)

    output:
    tuple val(meta), path("*.tsv"), emit: format_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    # Read krona file (count + tab-separated taxonomy)
    krona_data <- read_tsv("${krona_file}", 
        col_names = c("count", "taxonomy"),
        col_types = "dc",
        show_col_types = FALSE
    )

    # Calculate total reads for fraction calculation
    total_reads <- sum(krona_data\$count, na.rm = TRUE)

    # Parse taxonomy and create tabular format
    formatted_data <- krona_data %>%
        filter(count > 0) %>%
        mutate(
            fraction = count / total_reads,
            # Split taxonomy by tabs and extract levels
            tax_levels = str_split(taxonomy, "\\t")
        ) %>%
        # Create separate columns for each taxonomic rank
        rowwise() %>%
        mutate(
            superkingdom = ifelse(length(tax_levels) >= 1, tax_levels[[1]], 'unclassified'),
            phylum = ifelse(length(tax_levels) >= 2, tax_levels[[2]], 'unclassified'),
            class = ifelse(length(tax_levels) >= 3, tax_levels[[3]], 'unclassified'),
            order = ifelse(length(tax_levels) >= 4, tax_levels[[4]], 'unclassified'),
            family = ifelse(length(tax_levels) >= 5, tax_levels[[5]], 'unclassified'),
            genus = ifelse(length(tax_levels) >= 6, tax_levels[[6]], 'unclassified'),
            species = ifelse(length(tax_levels) >= 7, tax_levels[[7]], 'unclassified')
        ) %>%
        ungroup() %>%
        # Clean up taxonomy names (remove prefixes like k__, p__, etc.)
        mutate(
            superkingdom = str_remove(superkingdom, "^[kd]__") %>% str_replace_all("_", " "),
            phylum = str_remove(phylum, "^p__") %>% str_replace_all("_", " "),
            class = str_remove(class, "^c__") %>% str_replace_all("_", " "),
            order = str_remove(order, "^o__") %>% str_replace_all("_", " "),
            family = str_remove(family, "^f__") %>% str_replace_all("_", " "),
            genus = str_remove(genus, "^g__") %>% str_replace_all("_", " "),
            species = str_remove(species, "^s__") %>% str_replace_all("_", " ")
        ) %>%
        # Select final columns
        select(fraction, superkingdom, phylum, class, order, family, genus, species)

    # Write output
    write_tsv(formatted_data, "${prefix}.tsv")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    tidyverse: ", packageVersion('tidyverse'))
        ),
        "versions.yml"
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat > ${prefix}.tsv << 'EOF'
    fraction	superkingdom	phylum	class	order	family	genus	species
    0.06582267833109018	Eukaryota	Chordata	Mammalia	Artiodactyla	Suidae	Sus	Sus scrofa
    0.0259084791386272	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia	Escherichia coli
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.3.1
        tidyverse: 2.0.0
    END_VERSIONS
    """
}