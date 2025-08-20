//
// Manipulate GTDB metadata tables, GTDB-Tk and CheckM output to extract the genome information we want.
//

workflow METADATA {
    take:
    gtdb_metadata_files
    gtdbtk_metadata_files
    checkm_metadata_files

    main:

    //
    // INPUT: Parse any GTDB metadata files
    //
    ch_gtdb_metadata = gtdb_metadata_files
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

    //
    // INPUT: GTDB-Tk metadata
    //
    ch_gtdbtk_metadata = gtdbtk_metadata_files
        .splitCsv( sep: '\t', header: true)
        .map { [ [it.user_genome],
            [
                gtdb_genome_representative: it.user_genome,
                gtdb_representative: "f",
                gtdb_taxonomy: it.classification
            ] ]
        }

    //
    // INPUT: Parse CheckM metadata
    //
    ch_checkm_metadata = checkm_metadata_files
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

    //
    // gtdbtk_metadata and checkm_metadata need to be joined
    //
    if ( ch_gtdbtk_metadata && ch_checkm_metadata && ch_gtdb_metadata ) {
        ch_checkm_gtdb_metadata = ch_gtdbtk_metadata
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

        ch_gtdbtk_checkm_filtered = ch_checkm_gtdb_metadata
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

        ch_metadata = ch_checkm_gtdb_metadata
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

    } else if ( !ch_gtdbtk_metadata && ch_checkm_metadata && ch_gtdb_metadata) {
        ch_checkm_metadata = ch_checkm_metadata
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

        ch_checkm_filtered = ch_checkm_metadata
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
            .map{ it[1] }

        ch_metadata = ch_checkm_metadata
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

    } else if( ch_gtdbtk_metadata && !ch_checkm_metadata && ch_gtdb_metadata) {
        ch_gtdbtk_metadata = ch_gtdbtk_metadata
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

        ch_gtdbtk_filtered = ch_gtdbtk_metadata
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
            .map{ it[1] }

        ch_metadata = ch_gtdbtk_metadata
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
            .map{ it[1] }
            .mix(ch_gtdbtk_filtered)

    } else if( ch_gtdbtk_metadata && ch_checkm_metadata && !ch_gtdb_metadata) {
        ch_metadata = ch_gtdbtk_metadata
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

    } else if( !ch_gtdbtk_metadata && !ch_checkm_metadata && !ch_gtdb_metadata) {
        ch_metadata = Channel.empty()

    } else if( !ch_gtdbtk_metadata && ch_checkm_metadata && !ch_gtdb_metadata) {
        ch_metadata = ch_checkm_metadata
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

    } else if( ch_gtdbtk_metadata && !ch_checkm_metadata && !ch_gtdb_metadata) {
        ch_metadata = ch_gtdbtk_metadata
        .map { accno, gtdbtk -> [ accno, gtdbtk ] }
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

    } else if( !ch_gtdbtk_metadata && !ch_checkm_metadata && ch_gtdb_metadata) {
        ch_metadata = ch_gtdb_metadata
    }

    emit:
    metadata = ch_metadata
}
