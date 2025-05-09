process COLLECT_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-dtplyr=1.3.1 conda-forge::r-data.table=1.14.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), val(samples), path(trimlogs), path(bbduklogs), path(idxstats), path(fcs)

    output:
    path "${meta.id}.overall_stats.tsv.gz", emit: overall_stats
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def read_trimlogs = "tibble(sample = character(), m = character(), v = numeric()) %>%"
    if ( trimlogs ) {
        read_trimlogs = """
    read_tsv(
        pipe("grep -H 'Reads written (passing filters)' *trimming_report.txt"),
        col_names = c('s'), col_types = 'c'
    ) %>%
        transmute(
            sample    = str_replace(s, '^(.+)_\\\\d\\\\.fastq.*', '\\\\1'),
            n_trimmed = str_replace(s, '.* (\\\\d+) .*', '\\\\1') %>% as.integer() * 2
        ) %>%
        pivot_longer(2:ncol(.), names_to = 'm', values_to = 'v') %>%
        """
    }

    def read_bbduklogs = "tibble(sample = character(), m = character(), v = numeric())"
    if ( bbduklogs ) {
        read_bbduklogs = """
        read_fwf(pipe("grep -H 'Result:' *.bbduk.log"), fwf_widths(c(1e5), c('c')), col_types = 'c') %>%
            transmute(
                sample = str_remove(c, '.bbduk.log.*'),
                m = 'n_non_contaminated',
                v = str_replace(c, '.*Result:\\\\D*(\\\\d+).*', '\\\\1') %>% as.numeric()
            )
        """
    }

    """
    #!/usr/bin/env Rscript

    #library(data.table)
    #library(dtplyr)
    library(dplyr)
    library(readr)
    library(purrr)
    library(tidyr)
    library(stringr)

    TYPE_ORDER = c('n_trimmed', 'n_non_contaminated', 'idxs_n_mapped', 'idxs_n_unmapped', 'CDS', 'rRNA', 'tRNA', 'tmRNA')

    t <- ${read_trimlogs}
        # add samtools idxstats output
        union(
            read_tsv(
                pipe("grep -vH '*' *.idxstats"),
                col_names = c('c', 'length', 'idxs_n_mapped', 'idxs_n_unmapped'),
                col_types = 'ciii'
            ) %>%
                separate(c, c('s', 'chr'), sep = ':') %>%
                mutate(sample = str_replace(s, '^(.*)\\\\.idxstats', '\\\\1')) %>%
                group_by(sample) %>%
                summarise(
                    idxs_n_mapped   = sum(idxs_n_mapped),
                    idxs_n_unmapped = sum(idxs_n_unmapped)
                ) %>%
                pivot_longer(2:ncol(.), names_to = 'm', values_to = 'v')
        ) %>%
        union(
            # Total observation after featureCounts
            read_tsv(Sys.glob('*.counts.tsv.gz'), col_types = 'cciicicid', id = 'fname') %>%
                mutate(m = str_replace(fname, '.*\\\\.(.+)\\\\.counts.tsv.gz', '\\\\1')) %>%
                group_by(sample, m) %>% summarise(v = sum(count), .groups = 'drop')
        ) %>%
        union(
            ${read_bbduklogs}
        )

    # Write the table in wide format
    t %>%
        mutate(m = parse_factor(m, levels = unique(c(TYPE_ORDER, t\$m)), ordered = TRUE)) %>%
        arrange(sample, m) %>%
        pivot_wider(names_from = m, values_from = v, values_fill = 0) %>%
        write_tsv('${prefix}.overall_stats.tsv.gz')

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    dtplyr: ", packageVersion('dtplyr')),
            paste0("    data.table: ", packageVersion('data.table')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    purrr: ", packageVersion('purrr')),
            paste0("    tidyr: ", packageVersion('tidyr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml"
    )
    """
}
