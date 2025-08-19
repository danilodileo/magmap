//
// Create a concatenated set of gff files. Handles large number of input files by concatenating in two passes.
//

include { CAT_MANY as CAT_GFF   } from '../../../modules/local/cat_many'
include { GENOMES2ORFS          } from '../../../modules/local/genomes2orfs'
include { PROKKAGFF2TSV         } from '../../../modules/local/prokkagff2tsv'

workflow CAT_GFFS {
    take:
    ch_genome_gffs

    main:
        ch_versions = Channel.empty()

        CAT_GFF([id:'genomes'], ch_genome_gffs.collect())
        ch_versions = ch_versions.mix(CAT_GFF.out.versions)

        GENOMES2ORFS(ch_genome_gffs.collect().map { gffs -> [ [ id: 'genomes' ], gffs ] })
        ch_versions = ch_versions.mix(GENOMES2ORFS.out.versions)

        PROKKAGFF2TSV(CAT_GFF.out.concatenated_files)
        ch_versions = ch_versions.mix(PROKKAGFF2TSV.out.versions)

    emit:
    gff          = CAT_GFF.out.concatenated_files
    genomes2orfs = GENOMES2ORFS.out.genomes2orfs
    gfftsv       = PROKKAGFF2TSV.out.tsv
    versions     = ch_versions
}
