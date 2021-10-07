/*
 * Include modules
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
 * Modules for Hybrid assemblies
 */

// Unicycler hybrid
include { unicycler_hybrid } from '../modules_batch/Hybrid/unicycler_hybrid.nf'

// Unicycler hybrid
include { haslr_hybrid } from '../modules_batch/Hybrid/haslr_hybrid.nf'

// SPAdes hybrid
include { spades_hybrid } from '../modules_batch/Hybrid/spades_hybrid.nf'

// Pilon polish paired
include { pilon_polish } from '../modules_batch/Hybrid/unicycler_polish.nf'

/*
 * Module for assessing assembly qualities
 */
include { quast }   from '../modules_batch/QualityAssessment/quast.nf'
include { multiqc } from '../modules_batch/QualityAssessment/multiqc.nf'

workflow hybrid_nf {
  take:
      input_tuple
  
  main:

      /*
       * Channels for placing the assemblies
       */
       // lreads
      lreads_assemblies_ch = Channel.empty()
      canu_ch              = Channel.empty()
      unicycler_ch         = Channel.empty()
      flye_ch              = Channel.empty()
      raven_ch             = Channel.empty()
      wtdbg2_ch            = Channel.empty()
      shasta_ch            = Channel.empty()

       // hybrid
      hybrid_assemblies_ch = Channel.empty()
      spades_ch            = Channel.empty()
      unicycler_h_ch       = Channel.empty()
      haslr_ch             = Channel.empty()

      /*
       * Channels for placing polished assemblies
       */
      medaka_ch     = Channel.empty()
      nanopolish_ch = Channel.empty()
      gcpp_ch       = Channel.empty()
      pilon_ch      = Channel.empty()

      /*
       * Full (default) hybrid mode
       */

      // SPAdes
      if (!params.skip_spades) {
        spades_hybrid(input_tuple.filter { it[1] == "hybrid-strategy-1" })
        spades_ch = spades_hybrid.out[1]
      }
      // Unicycler
      if (!params.skip_unicycler) {
        unicycler_hybrid(input_tuple.filter { it[1] == "hybrid-strategy-1" })
        unicycler_h_ch = unicycler_hybrid.out[1]
      }
      // Haslr
      if (!params.skip_haslr) {
        haslr_hybrid(input_tuple.filter { it[1] == "hybrid-strategy-1" })
        haslr_ch = haslr_hybrid.out[1]
      }

      // Get hybrid assemblies
      hybrid_assemblies_ch = spades_ch.mix(unicycler_h_ch, haslr_ch)

      /*
       * Polish a long reads assembly
       */
      
      /*
       * Canu
       */
      if (!params.skip_canu) {
        canu(input_tuple.filter { it[1] == "hybrid-strategy-2" })
        canu_ch = canu.out[1]
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        flye(input_tuple.filter { it[1] == "hybrid-strategy-2" })
        flye_ch = flye.out[1]
      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        unicycler(input_tuple.filter { it[1] == "hybrid-strategy-2" })
        unicycler_ch = unicycler.out[1]
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        raven(input_tuple.filter { it[1] == "hybrid-strategy-2" })
        raven_ch = raven.out[1]
      }

      /*
       * Shasta
       */
      if (!params.skip_shasta && params.lr_type == 'nanopore') {
        shasta(input_tuple.filter { it[1] == "hybrid-strategy-2" })
        shasta_ch = shasta.out[1]
      }

      /*
       * wtdbg2
       */
      if (!params.skip_wtdbg2) {
        wtdbg2(input_tuple.filter { it[1] == "hybrid-strategy-2" })
        wtdbg2_ch = wtdbg2.out[1]
      }

      // Get long reads assemblies
      lreads_assemblies_ch = canu_ch.mix(flye_ch, unicycler_ch, raven_ch, wtdbg2_ch, shasta_ch)

      // combine again with metadata
      lreads_combined_ch = lreads_assemblies_ch.combine(input_tuple, by: 0)

      /*
       * Run medaka?
       */
      medaka(lreads_combined_ch)
      medaka_ch = medaka.out[1]

      /*
       * Run nanopolish?
       */
      nanopolish(lreads_combined_ch)
      nanopolish_ch = nanopolish.out[0]

      /*
       * gcpp?
       */
      gcpp(lreads_combined_ch)
      gcpp_ch = gcpp.out[1]

      /*
       * Finally, run pilon for all
       */
      pilon_combined_ch = lreads_assemblies_ch.mix(medaka_ch, nanopolish_ch, gcpp_ch).combine(input_tuple, by: 0)
      pilon_polish(pilon_combined_ch)
      pilon_ch = pilon_polish.out[1]
      

      // Run quast (with all)
      quast(
        hybrid_assemblies_ch.mix(
          lreads_assemblies_ch, 
          medaka_ch, nanopolish_ch, 
          gcpp_ch, 
          pilon_ch).combine(
          preads.combine(sreads).collect().toList()
        ).combine(
          prefix_ch
        )
      )

      // Run multiqc
      multiqc(quast.out[0].collect(), prefix_ch, Channel.value("$workflow.runName"))

}
