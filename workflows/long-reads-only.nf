/*
 * Modules for assembling long reads
 */

// Canu assembler
include { canu_assembly } from '../modules/LongReads/canu.nf'

// Unicycler assembler
include { unicycler_lreads_assembly } from '../modules/LongReads/unicycler_lreads.nf'

// Flye assembler
include { flye_assembly } from '../modules/LongReads/flye.nf'

// Raven assembler
include { raven_assembly } from '../modules/LongReads/raven.nf'

// wtdbg2 assembler
include { wtdbg2_assembly } from '../modules/LongReads/wtdbg2.nf'

// Shasta assembler
include { shasta_assembly } from '../modules/LongReads/shasta.nf'

/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include { nanopolish } from '../modules/LongReads/nanopolish.nf'

// Medaka (for nanopore data)
include { medaka } from '../modules/LongReads/medaka.nf'

// gcpp Pacbio
include { gcpp } from '../modules/LongReads/gcpp.nf'

/*
 * Module for assessing assembly qualities
 */
include { quast } from '../modules/QualityAssessment/quast.nf'
include { multiqc } from '../modules/QualityAssessment/multiqc.nf'

workflow lreadsonly_nf {
  take:
      reads
      fast5
      fast5_dir
      bamFile
      nBams
  main:

      /*
       * Channels for placing the assemblies
       */
      canu_ch      = Channel.empty()
      unicycler_ch = Channel.empty()
      flye_ch      = Channel.empty()
      raven_ch     = Channel.empty()
      wtdbg2_ch    = Channel.empty()
      shasta_ch    = Channel.empty()

      /*
       * Channels for placing polished assemblies
       */
      medaka_ch     = Channel.empty()
      nanopolish_ch = Channel.empty()
      gcpp_ch       = Channel.empty()


      /*
       * Canu
       */
      if (!params.skip_canu) {
        canu_assembly(reads)
        canu_ch = canu_assembly.out[1]
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        flye_assembly(reads)
        flye_ch = flye_assembly.out[1]
      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        unicycler_lreads_assembly(reads)
        unicycler_ch = unicycler_lreads_assembly.out[1]
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        raven_assembly(reads)
        raven_ch = raven_assembly.out[1]
      }

      /*
       * Shasta
       */
      if (!params.skip_shasta && params.lr_type == 'nanopore') {
        shasta_assembly(reads)
        shasta_ch = shasta_assembly.out[1]
      }

      /*
       * wtdbg2
       */
      if (!params.skip_wtdbg2) {
        wtdbg2_assembly(reads)
        wtdbg2_ch = wtdbg2_assembly.out[1]
      }

      // gather assemblies
      assemblies_ch = canu_ch.mix(unicycler_ch, flye_ch, raven_ch, wtdbg2_ch, shasta_ch)

      /*
       * Run medaka?
       */
      if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
        medaka(assemblies_ch.combine(reads))
        medaka_ch = medaka.out[1]
      }

      /*
       * Run nanopolish?
       */
      if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
        nanopolish(assemblies_ch.combine(reads).combine(fast5).combine(fast5_dir))
        nanopolish_ch = nanopolish.out[0]
      }

      /*
       * gcpp?
       */
      if (params.pacbio_bams && params.lr_type == 'pacbio') {
        gcpp(assemblies_ch.combine(bamFile.collect().toList()).combine(nBams))
        gcpp_ch = gcpp.out[1]
      }

      // Gather polishings
      polished_ch = medaka_ch.mix(nanopolish_ch, gcpp_ch)

      /*
       * Run quast
       */
      quast(assemblies_ch.mix(polished_ch).combine(reads), reads)

      /*
       * Run multiqc
       */
      multiqc(quast.out[0].collect(), quast.out[1].distinct(), quast.out[2], Channel.value("$workflow.runName"))
}
