/*
 * Include modules
 */

/*
 * Module for prefix evaluation
 */
include { define_prefix } from '../modules/misc/define_prefix.nf'

/*
 * Modules for assembling long reads
 */

// Canu assembler
include { canu } from '../modules/LongReads/canu.nf'

// Unicycler assembler
include { unicycler } from '../modules/LongReads/unicycler.nf'

// Flye assembler
include { flye } from '../modules/LongReads/flye.nf'

// Raven assembler
include { raven } from '../modules/LongReads/raven.nf'

// wtdbg2 assembler
include { wtdbg2 } from '../modules/LongReads/wtdbg2.nf'

// Shasta assembler
include { shasta } from '../modules/LongReads/shasta.nf'

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
 * Modules for Hybrid assemblies
 */

// Unicycler hybrid
include { unicycler_hybrid } from '../modules/Hybrid/unicycler_hybrid.nf'

// Unicycler hybrid
include { haslr_hybrid } from '../modules/Hybrid/haslr_hybrid.nf'

// SPAdes hybrid
include { spades_hybrid } from '../modules/Hybrid/spades_hybrid.nf'

// Pilon polish paired
include { pilon_polish } from '../modules/Hybrid/unicycler_polish.nf'

/*
 * Module for assessing assembly qualities
 */
include { quast } from '../modules/QualityAssessment/quast.nf'
include { multiqc } from '../modules/QualityAssessment/multiqc.nf'

workflow hybrid_nf {
  take:
      preads
      sreads
      lreads
      fast5
      fast5_dir
      bamFile
      nBams
  
  main:

      /*
       * Check input reads in order to evaluate better prefix
       */
      define_prefix(lreads, preads, sreads)
      prefix_ch = define_prefix.out[0]

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

      if (!params.strategy_2) {

        // SPAdes
        if (!params.skip_spades) {
          spades_hybrid(lreads, preads, sreads, prefix_ch)
          spades_ch = spades_hybrid.out[1]
        }
        // Unicycler
        if (!params.skip_unicycler) {
          unicycler_hybrid(lreads, preads, sreads, prefix_ch)
          unicycler_h_ch = unicycler_hybrid.out[1]
        }
        // Haslr
        if (!params.skip_haslr) {
          haslr_hybrid(lreads, preads, sreads, prefix_ch)
          haslr_ch = haslr_hybrid.out[1]
        }

        // Get hybrid assemblies
        hybrid_assemblies_ch = spades_ch.mix(unicycler_h_ch, haslr_ch)

      }

      /*
       * Polish a long reads assembly
       */

      if (params.strategy_2) {
        /*
         * Canu
         */
        if (!params.skip_canu) {
          canu(lreads.combine(prefix_ch))
          canu_ch = canu.out[1]
        }

        /*
         * Flye
         */
        if (!params.skip_flye) {
          flye(lreads.combine(prefix_ch))
          flye_ch = flye.out[1]
        }

        /*
         * Unicycler
         */
        if (!params.skip_unicycler) {
          unicycler(lreads.combine(prefix_ch))
          unicycler_ch = unicycler.out[1]
        }

        /*
         * Raven
         */
        if (!params.skip_raven) {
          raven(lreads.combine(prefix_ch))
          raven_ch = raven.out[1]
        }

        /*
         * Shasta
         */
        if (!params.skip_shasta && params.lr_type == 'nanopore') {
          shasta(lreads.combine(prefix_ch))
          shasta_ch = shasta.out[1]
        }

        /*
         * wtdbg2
         */
        if (!params.skip_wtdbg2) {
          wtdbg2(lreads.combine(prefix_ch))
          wtdbg2_ch = wtdbg2.out[1]
        }

        // Get long reads assemblies
        lreads_assemblies_ch = canu_ch.mix(flye_ch, unicycler_ch, raven_ch, wtdbg2_ch, shasta_ch)

        /*
         * Run medaka?
         */
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka(lreads_assemblies_ch.combine(lreads).combine(prefix_ch))
          medaka_ch = medaka.out[1]
        }

        /*
         * Run nanopolish?
         */
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish(
            lreads_assemblies_ch.combine(lreads).combine(fast5).combine(fast5_dir).combine(prefix_ch),
          )
          nanopolish_ch = nanopolish.out[0]
        }

        /*
         * gcpp?
         */
        if (params.pacbio_bams && params.lr_type == 'pacbio') {
          gcpp(
            lreads_assemblies_ch.combine(bamFile.collect().toList()).combine(nBams).combine(prefix_ch),
          )
          gcpp_ch = gcpp.out[1]
        }

        /*
         * Finally, run pilon for all
         */
        pilon_polish(
          lreads_assemblies_ch.mix(
            medaka_ch, 
            nanopolish_ch, 
            gcpp_ch).combine(
            preads.combine(sreads).collect().toList()
          ).combine(
            prefix_ch
          ),
        )
        pilon_ch = pilon_polish.out[1]
      }

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
