//
// Create a metadata table that contains information from checkM/checkM2, GTDB and GTDB-tk
//

workflow METADATA {
    take:

    main:

    //
    // INPUT: if user provides, populate ch_metadata
    //
    ch_gtdb_metadata = Channel.empty()

    if ( params.gtdb_metadata) {
        Channel
            .of(params.gtdb_metadata.split(','))
            .map { file(it) }
            .splitCsv( sep: '\t', header: true)
            .map {
                [
                    accno: it.accession - ~/^[A-Z][A-Z]_/,
                    checkm_completeness: it.checkm_completeness,
                    checkm_contamination: it.checkm_contamination,
                    checkm_strain_heterogeneity: it.checkm_strain_heterogeneity,
                    contig_count: it.contig_count,
                    genome_size: it.genome_size,
                    gtdb_genome_representative: it.gtdb_genome_representative,
                    gtdb_representative: it.gtdb_representative,
                    gtdb_taxonomy: it.gtdb_taxonomy
                ]
            }
            .set { ch_gtdb_metadata }
    }

    //
    // INPUT: CheckM metadata
    //
    ch_checkm_metadata = Channel.empty()
    if ( params.checkm_metadata) {
        Channel
            .of(params.checkm_metadata.split(','))
            .map { file(it) }
            .splitCsv( sep: '\t', header: true)
            .map { [ [ it["Bin Id"] ],
                [
                    checkm_completeness: it.Completeness,
                    checkm_contamination: it.Contamination,
                    contig_count: it["# contigs"],
                    checkm_strain_heterogeneity: it["Strain heterogeneity"],
                    genome_size: it["Genome size (bp)"]
                ] ]
            }
            .set { ch_checkm_metadata }
    }

    //
    // INPUT: GTDB-Tk metadata
    //
    ch_gtdbtk_metadata = Channel.empty()
    if ( params.gtdbtk_metadata) {
        Channel
            .of(params.gtdbtk_metadata.split(','))
            .map { file(it) }
            .splitCsv( sep: '\t', header: true)
            .map { [ [it.user_genome],
                [
                    gtdb_genome_representative: it.user_genome,
                    gtdb_representative: "f",
                    gtdb_taxonomy: it.classification
                ] ]
            }
            .set { ch_gtdbtk_metadata }
    }

    //
    // gtdbtk_metadata and checkm_metadata need to be joined
    //
    if ( params.gtdbtk_metadata && params.checkm_metadata && params.gtdb_metadata) {
        ch_gtdbtk_metadata
        .map{ accno, gtdbtk -> [ accno, gtdbtk ]}
        .join( ch_checkm_metadata
            .map{ accno, checkm -> [ accno, checkm ] }
        )
        .map { accno, gtdbtk, checkm ->
                [
                    accno,
                    [
                    accno: accno[0],
                    checkm_completeness: checkm.checkm_completeness,
                    checkm_contamination: checkm.checkm_contamination,
                    checkm_strain_heterogeneity: checkm.checkm_strain_heterogeneity,
                    contig_count: checkm.contig_count,
                    genome_size: checkm.genome_size,
                    gtdb_genome_representative: gtdbtk.gtdb_genome_representative,
                    gtdb_representative: gtdbtk.gtdb_representative,
                    gtdb_taxonomy: gtdbtk.gtdb_taxonomy
                    ]
                ]
            }
    .set { ch_checkm_gtdb_metadata }

    ch_checkm_gtdb_metadata
            .map {
                accno, gtdbtk_checkm -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 !in it }
            .map{ it[0] }
            .join(
                ch_checkm_gtdb_metadata.map { accno, gtdbtk_checkm -> [ gtdbtk_checkm.accno, gtdbtk_checkm ] }
            )
            .map{ it[1]}
        .set{ ch_gtdbtk_checkm_filtered }

        ch_checkm_gtdb_metadata
            .map {
                accno, gtdbtk_checkm -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 in it }
            .map{ it[0] }
            .join(
                ch_gtdb_metadata.map { gtdb -> [ gtdb.accno, gtdb ] }
            )
            .map{ it[1]}
            .mix(ch_gtdbtk_checkm_filtered)
            .set{ ch_metadata }

    } else if ( !params.gtdbtk_metadata && params.checkm_metadata && params.gtdb_metadata) {
        ch_checkm_metadata
            .map{ accno, checkm -> [
                    accno,
                    [
                    accno: accno[0],
                    checkm_completeness: checkm.checkm_completeness,
                    checkm_contamination: checkm.checkm_contamination,
                    checkm_strain_heterogeneity: checkm.checkm_strain_heterogeneity,
                    contig_count: checkm.contig_count,
                    genome_size: checkm.genome_size,
                    gtdb_genome_representative: "",
                    gtdb_representative: "",
                    gtdb_taxonomy: ""
                    ]
                ]
            }
            .set { ch_checkm_metadata }

        ch_checkm_metadata
            .map {
                accno, checkm -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 !in it }
            .map{ it[0] }
            .join(
                ch_checkm_metadata.map { accno, checkm -> [ checkm.accno, checkm ] }
            )
            .map{ it[1]}
        .set{ ch_checkm_filtered }

        ch_checkm_metadata
            .map {
                accno, checkm -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 in it }
            .map{ it[0] }
            .join(
                ch_gtdb_metadata.map { gtdb -> [ gtdb.accno, gtdb ] }
            )
            .map{ it[1] }
            .mix(ch_checkm_filtered)
        .set{ ch_metadata }

    } else if( params.gtdbtk_metadata && !params.checkm_metadata && params.gtdb_metadata) {
        ch_gtdbtk_metadata
        .map{ accno, gtdbtk -> [ accno, gtdbtk ] }
        .map { accno, gtdbtk ->
                [
                    accno,
                    [
                    accno: accno[0],
                    checkm_completeness: "",
                    checkm_contamination: "",
                    checkm_strain_heterogeneity: "",
                    contig_count: "",
                    genome_size: "",
                    gtdb_genome_representative: gtdbtk.gtdb_genome_representative,
                    gtdb_representative: gtdbtk.gtdb_representative,
                    gtdb_taxonomy: gtdbtk.gtdb_taxonomy
                    ]
                ]
            }
        .set { ch_gtdbtk_metadata }

    ch_gtdbtk_metadata
            .map {
                accno, gtdbtk -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 !in it }
            .map{ it[0] }
            .join(
                ch_gtdbtk_metadata.map { accno, gtdbtk -> [ gtdbtk.accno, gtdbtk ] }
            )
            .map{ it[1]}
        .set{ ch_gtdbtk_filtered }

        ch_gtdbtk_metadata
            .map {
                accno, gtdbtk -> accno[0]
            }
            .join(
                ch_gtdb_metadata.map { it -> [ it.accno, 1 ] }, remainder: true
            )
            .filter{ 1 in it }
            .map{ it[0] }
            .join(
                ch_gtdb_metadata.map { gtdb -> [ gtdb.accno, gtdb ] }
            )
            .map{ it[1]}
            .mix(ch_gtdbtk_filtered)
        .set{ ch_metadata }
    } else if( params.gtdbtk_metadata && params.checkm_metadata && !params.gtdb_metadata) {
        ch_gtdbtk_metadata
            .map{ accno, gtdbtk -> [ accno, gtdbtk ] }
            .join( ch_checkm_metadata
            .map{ accno, checkm -> [ accno, checkm ] }
            )
            .map { accno, gtdbtk, checkm ->
                    [
                    accno: accno[0],
                    checkm_completeness: checkm.checkm_completeness,
                    checkm_contamination: checkm.checkm_contamination,
                    checkm_strain_heterogeneity: checkm.checkm_strain_heterogeneity,
                    contig_count: checkm.contig_count,
                    genome_size: checkm.genome_size,
                    gtdb_genome_representative: gtdbtk.gtdb_genome_representative,
                    gtdb_representative: gtdbtk.gtdb_representative,
                    gtdb_taxonomy: gtdbtk.gtdb_taxonomy
                    ]
            }
        .set { ch_metadata }
    } else if( !params.gtdbtk_metadata && !params.checkm_metadata && !params.gtdb_metadata) {
        ch_metadata = Channel.empty()
    } else if( !params.gtdbtk_metadata && params.checkm_metadata && !params.gtdb_metadata) {
        ch_checkm_metadata
            .map{ accno, checkm ->
                [
                    accno: accno[0],
                    checkm_completeness: checkm.checkm_completeness,
                    checkm_contamination: checkm.checkm_contamination,
                    checkm_strain_heterogeneity: checkm.checkm_strain_heterogeneity,
                    contig_count: checkm.contig_count,
                    genome_size: checkm.genome_size,
                    gtdb_genome_representative: "",
                    gtdb_representative: "",
                    gtdb_taxonomy: ""
                ]
            }
            .set { ch_metadata }
    } else if( params.gtdbtk_metadata && !params.checkm_metadata && !params.gtdb_metadata) {
        ch_gtdbtk_metadata
        .map{ accno, gtdbtk -> [ accno, gtdbtk ] }
        .map { accno, gtdbtk ->
                [
                    accno: accno[0],
                    checkm_completeness: "",
                    checkm_contamination: "",
                    checkm_strain_heterogeneity: "",
                    contig_count: "",
                    genome_size: "",
                    gtdb_genome_representative: gtdbtk.gtdb_genome_representative,
                    gtdb_representative: gtdbtk.gtdb_representative,
                    gtdb_taxonomy: gtdbtk.gtdb_taxonomy
                ]
            }
        .set { ch_metadata }
    } else if( !params.gtdbtk_metadata && !params.checkm_metadata && params.gtdb_metadata) {
        ch_metadata = ch_gtdb_metadata
    }

    emit:
    metadata = ch_metadata
}
