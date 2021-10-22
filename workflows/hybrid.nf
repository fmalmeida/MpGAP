/*
 * Include modules
 */

/*
 * Modules for assembling long reads
 */

// Canu assembler
include { canu as strategy_2_canu } from '../modules/LongReads/canu.nf'

// Unicycler assembler
include { unicycler as strategy_2_unicycler } from '../modules/LongReads/unicycler.nf'

// Flye assembler
include { flye as strategy_2_flye } from '../modules/LongReads/flye.nf'

// Raven assembler
include { raven as strategy_2_raven } from '../modules/LongReads/raven.nf'

// wtdbg2 assembler
include { wtdbg2 as strategy_2_wtdbg2 } from '../modules/LongReads/wtdbg2.nf'

// Shasta assembler
include { shasta as strategy_2_shasta } from '../modules/LongReads/shasta.nf'

/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include { nanopolish as strategy_2_nanopolish } from '../modules/LongReads/nanopolish.nf'

// Medaka (for nanopore data)
include { medaka as strategy_2_medaka } from '../modules/LongReads/medaka.nf'

// gcpp Pacbio
include { gcpp as strategy_2_gcpp } from '../modules/LongReads/gcpp.nf'

/*
 * Modules for Hybrid assemblies
 */

// Unicycler hybrid
include { unicycler_hybrid as strategy_1_unicycler } from '../modules/Hybrid/unicycler_hybrid.nf'

// Unicycler hybrid
include { haslr_hybrid as strategy_1_haslr } from '../modules/Hybrid/haslr_hybrid.nf'

// SPAdes hybrid
include { spades_hybrid as strategy_1_spades } from '../modules/Hybrid/spades_hybrid.nf'

// Pilon polish paired
include { pilon_polish as strategy_2_pilon } from '../modules/Hybrid/unicycler_polish.nf'

/*
 * Module for assessing assembly qualities
 */
include { quast }   from '../modules/QualityAssessment/quast.nf'
include { multiqc } from '../modules/QualityAssessment/multiqc.nf'

workflow HYBRID {
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
       * create branches
       */
      input_tuple.branch{
        main: it[1] == "hybrid_strategy_1"
        secondary: it[1] == "hybrid_strategy_2"
      }.set { input_branches }

      /*
       * Full (default) hybrid mode
       */
      
      // SPAdes
      if (!params.skip_spades) {
        strategy_1_spades(input_branches.main)
        spades_ch =strategy_1_spades.out[1]
      }
      // Unicycler
      if (!params.skip_unicycler) {
        strategy_1_unicycler(input_branches.main)
        unicycler_h_ch =strategy_1_unicycler.out[1]
      }
      // Haslr
      if (!params.skip_haslr) {
        strategy_1_haslr(input_branches.main)
        haslr_ch =strategy_1_haslr.out[1]
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
        strategy_2_canu(input_branches.secondary)
        canu_ch = strategy_2_canu.out[1]
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        strategy_2_flye(input_branches.secondary)
        flye_ch = strategy_2_flye.out[1]
      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        strategy_2_unicycler(input_branches.secondary)
        unicycler_ch = strategy_2_unicycler.out[1]
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        strategy_2_raven(input_branches.secondary)
        raven_ch = strategy_2_raven.out[1]
      }

      /*
       * Shasta
       */
      if (!params.skip_shasta && params.lr_type == 'nanopore') {
        strategy_2_shasta(input_branches.secondary)
        shasta_ch = strategy_2_shasta.out[1]
      }

      /*
       * wtdbg2
       */
      if (!params.skip_wtdbg2) {
        strategy_2_wtdbg2(input_branches.secondary)
        wtdbg2_ch = strategy_2_wtdbg2.out[1]
      }

      // Get long reads assemblies
      lreads_assemblies_ch = canu_ch.mix(flye_ch, unicycler_ch, raven_ch, wtdbg2_ch, shasta_ch)

      // combine again with metadata
      lreads_combined_ch = lreads_assemblies_ch.combine(input_tuple, by: 0)

      /*
       * Run medaka?
       */
      strategy_2_medaka(lreads_combined_ch)
      medaka_ch = strategy_2_medaka.out[1]

      /*
       * Run nanopolish?
       */
      strategy_2_nanopolish(lreads_combined_ch)
      nanopolish_ch = strategy_2_nanopolish.out[0]

      /*
       * gcpp?
       */
      strategy_2_gcpp(lreads_combined_ch)
      gcpp_ch = strategy_2_gcpp.out[1]

      /*
       * Finally, run pilon for all
       */
      pilon_combined_ch = lreads_assemblies_ch.mix(medaka_ch, nanopolish_ch, gcpp_ch)
      strategy_2_pilon(pilon_combined_ch.combine(input_tuple, by: 0))
      pilon_ch = strategy_2_pilon.out[1]

      // Run quast (with all)
      quast(hybrid_assemblies_ch.mix(pilon_combined_ch).mix(pilon_ch).combine(input_tuple, by: 0))

      // Run multiqc
      multiqc(quast.out[0].groupTuple(), quast.out[1], Channel.value("$workflow.runName"))

}
