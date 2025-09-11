//
// Select genomes from the downstream analysis that are identified by Sourmash
//

include { SOURMASH_GATHER                   } from '../../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH as GENOMES_SKETCH } from '../../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_INDEX  as GENOMES_INDEX  } from '../../../modules/nf-core/sourmash/index/main'
include { SOURMASH_SKETCH as SAMPLES_SKETCH } from '../../../modules/nf-core/sourmash/sketch/main'

workflow SOURMASH {
    take:
        ch_sample_reads             // Fastq files with reads for each sample [ val(meta), [ path(reads) ] ]
        ch_indexes                  // List of Sourmash indexs [ path(index) ]
        ch_user_genomeinfo          // User provided genomes [ path(genome) ]
        ch_remote_genome_sources    // Paths to genome information in NCBI format, i.e. containing at least the assembly_accession and ftp_path fields: path(csvfile)
        ksize                       // K-mere size to use: val(odd_int)
        save_unassigned             // Boolean value passed to sourmash/gather
        save_matches_sig            // Boolean value passed to sourmash/gather
        save_prefetch               // Boolean value passed to sourmash/gather
        save_prefetch_csv           // Boolean value passed to sourmash/gather

    main:
        ch_versions = Channel.empty()

        ch_ncbi_genomeinfo = ch_remote_genome_sources
                .splitCsv(skip: 1, header: true, sep: '\t')
                .map { row ->
                    [
                        accno: row["#assembly_accession"],
                        genome_fna: "${row.ftp_path}/${row.ftp_path - ~/\/$/ - ~/.*\//}_genomic.fna.gz",
                        genome_gff: ""
                    ]
                }

        GENOMES_SKETCH(ch_user_genomeinfo.map { [ [ id: it.accno ], it.genome_fna ] })
        ch_versions = ch_versions.mix(GENOMES_SKETCH.out.versions)

        SAMPLES_SKETCH(ch_sample_reads)
        ch_versions = ch_versions.mix(SAMPLES_SKETCH.out.versions)

        ch_sample_sigs = SAMPLES_SKETCH.out.signatures
            .collect { it[1] }
            .map { [ [ id: 'samples_sig' ], it ] }

        ch_genome_sigs = GENOMES_SKETCH.out.signatures
            .collect { meta, sig -> [ sig ] }
            .map { sig -> [ [ id: 'signatures' ], sig ] }

        GENOMES_INDEX(ch_genome_sigs, ksize)
        ch_versions = ch_versions.mix(GENOMES_INDEX.out.versions)

        ch_database = GENOMES_INDEX.out.signature_index
            .map{ meta, sig -> sig }
            .mix( ch_indexes )
            .collect()

        SOURMASH_GATHER(ch_sample_sigs, ch_database, save_unassigned, save_matches_sig, save_prefetch, save_prefetch_csv )
        ch_versions = ch_versions.mix(SOURMASH_GATHER.out.versions)

        ch_accnos_ncbi = SOURMASH_GATHER.out.result
            .map { meta, csv -> csv }
            .splitCsv( sep: ',', header: true, quote: '"')
            .map { it.name }
            .unique()
            .map { name ->
                def matcher = (name =~ /(GCA_[0-9]+\.[0-9]+|GCF_[0-9]+\.[0-9]+)/)
                if (matcher) {
                    return matcher[0][0] // Return the matched pattern
                }
            }

        ch_all_non_ncbi_user_accnos = SOURMASH_GATHER.out.result
            .map{ meta, csv -> csv }
            .splitCsv( sep: ',', header: true, quote: '"')
            .map { row -> row.name }
            .unique()
            .filter { name ->
                !(name =~ /(GCA_[0-9]+\.[0-9]+|GCF_[0-9]+\.[0-9]+)/)
            }


        // Subset the two genome info channels to only contain those that Sourmash identified
        // The user supplied channel takes precedence, so start with that
        ch_matching_user_non_ncbi_genomes = ch_all_non_ncbi_user_accnos
            .join(ch_user_genomeinfo.map { [ it.accno, [ it ] ]} )
            .map { it[1][0] }

        ch_matching_user_ncbi_genomes = ch_accnos_ncbi
            .join(ch_user_genomeinfo.map { [ it.accno, [ it ] ]} )
            .map { it[1][0] }

        ch_filtered_genomes = ch_accnos_ncbi
            .map { accno -> [accno, null] } // Initialize the channel with accno and null
            .join(ch_matching_user_ncbi_genomes.map { [it.accno, [ it ] ] }, remainder: true) // Perform the join
            .flatMap { tuple ->

                def accno = tuple[0] // accno from the left channel
                def matches = tuple[2] // Should be a list or null

                if (matches == null || matches.isEmpty()) {
                    return [[accno, null]]
                } else {
                    return matches.collect { match -> [accno, match] }
                }
            }
            .filter { it[1] == null } // Keep only tuples with null data
            .map { it[0] } // Extract accno
            .join(ch_ncbi_genomeinfo.map { [it.accno, it ] }) // Join with genome info
            .flatMap { tuple ->

                    def accno = tuple[0] // accno from the left channel
                    def genomeInfos = tuple[1] // Should be a list or null
                    if (genomeInfos == null || genomeInfos.isEmpty()) {
                        return [] // Discard tuples with null or empty genomeInfos
                    }

                    // Extract values from the map
                    def genomeFna = genomeInfos.genome_fna ?: ''
                    def genomeGff = genomeInfos.genome_gff ?: ''

                // Return the formatted result
                return [[accno: accno, genome_fna: genomeFna, genome_gff: genomeGff]]
            }
            .flatten()
            .mix(ch_matching_user_ncbi_genomes) // Ensure proper mixing with other data
            .mix(ch_matching_user_non_ncbi_genomes)
            .view { "ch_filtered_genomes: ${it}" }

    emit:
        gindex           = GENOMES_SKETCH.out.signatures
        sindex           = SAMPLES_SKETCH.out.signatures
        filtered_genomes = ch_filtered_genomes
        versions         = ch_versions
}
