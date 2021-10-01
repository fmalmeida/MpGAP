/*
 * DEFINITION OF MULTI-SAMPLE (BATCH) MODE
 */

/*
 * Modules for assembling long reads
 */

// Canu assembler
include { canu } from '../modules_batch/LongReads/canu.nf'

// Unicycler assembler
include { unicycler } from '../modules_batch/LongReads/unicycler.nf'

// Flye assembler
include { flye } from '../modules_batch/LongReads/flye.nf'

// Raven assembler
include { raven } from '../modules_batch/LongReads/raven.nf'

// wtdbg2 assembler
include { wtdbg2 } from '../modules_batch/LongReads/wtdbg2.nf'

// Shasta assembler
include { shasta } from '../modules_batch/LongReads/shasta.nf'

/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include { nanopolish } from '../modules_batch/LongReads/nanopolish.nf'

// Medaka (for nanopore data)
include { medaka } from '../modules_batch/LongReads/medaka.nf'

// gcpp Pacbio
include { gcpp } from '../modules_batch/LongReads/gcpp.nf'

/*
 * Module for assessing assembly qualities
 */
include { quast   } from '../modules_batch/QualityAssessment/quast.nf'
include { multiqc } from '../modules_batch/QualityAssessment/multiqc.nf'

workflow lreadsonly_nf {
  take:
      input_tuple

  main:

      /*
       * Channels for placing the assemblies
       */
      canu_ch       = Channel.empty()
      unicycler_ch  = Channel.empty()
      flye_ch       = Channel.empty()
      raven_ch      = Channel.empty()
      wtdbg2_ch     = Channel.empty()
      shasta_ch     = Channel.empty()

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
        canu(input_tuple)
        canu_ch = canu.out[1]
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        flye(input_tuple)
        flye_ch = flye.out[1]
      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        unicycler(input_tuple)
        unicycler_ch = unicycler.out[1]
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        raven(input_tuple)
        raven_ch = raven.out[1]
      }

      /*
       * Shasta
       */
      if (!params.skip_shasta) {
        shasta(input_tuple)
        shasta_ch = shasta.out[1]
      }

      /*
       * wtdbg2
       */
      if (!params.skip_wtdbg2) {
        wtdbg2(input_tuple)
        wtdbg2_ch = wtdbg2.out[1]
      }

      // gather assemblies
      assemblies_ch = canu_ch.mix(unicycler_ch, flye_ch, raven_ch, wtdbg2_ch, shasta_ch)

      /*
       * Run medaka?
       */
      medaka(assemblies_ch.combine(input_tuple, by: 0))
      medaka_ch = medaka.out[1]

      /*
       * Run nanopolish?
       */
      nanopolish(assemblies_ch.combine(input_tuple, by: 0))
      nanopolish_ch = nanopolish.out[0]

      /*
       * gcpp?
       */
      gcpp(assemblies_ch.combine(input_tuple, by: 0))
      gcpp_ch = gcpp.out[1]

      // Gather polishings
      polished_ch = medaka_ch.mix(nanopolish_ch, gcpp_ch)

      /*
       * Run quast
       */
      quast(assemblies_ch.mix(polished_ch).combine(input_tuple, by: 0))

      /*
       * Run multiqc
       */
      multiqc(quast.out[0].collect(), prefix_ch, Channel.value("$workflow.runName"))
}
