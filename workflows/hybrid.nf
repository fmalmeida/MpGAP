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
include { pilon_polish as strategy_2_pilon } from '../modules/Hybrid/pilon_polish.nf'

workflow HYBRID {
  take:
      input_tuple
  
  main:

      // Define default output channes
      // default must be a empty channel that
      // will be overwritten if assembler is used
      def LONGREADS_OUTPUTS = [:]
      LONGREADS_OUTPUTS['CANU']        = Channel.empty()
      LONGREADS_OUTPUTS['UNICYCLER']   = Channel.empty()
      LONGREADS_OUTPUTS['FLYE']        = Channel.empty()
      LONGREADS_OUTPUTS['RAVEN']       = Channel.empty()
      LONGREADS_OUTPUTS['WTDBG2']      = Channel.empty()
      LONGREADS_OUTPUTS['SHASTA']      = Channel.empty()
      LONGREADS_OUTPUTS['MEDAKA']      = Channel.empty()
      LONGREADS_OUTPUTS['NANOPOLISH']  = Channel.empty()
      LONGREADS_OUTPUTS['GCPP']        = Channel.empty()

      def HYBRID_OUTPUTS = [:]
      HYBRID_OUTPUTS['UNICYCLER'] = Channel.empty()
      HYBRID_OUTPUTS['SPADES']    = Channel.empty()
      HYBRID_OUTPUTS['HASLR']     = Channel.empty()
      HYBRID_OUTPUTS['PILON']     = Channel.empty()

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
        HYBRID_OUTPUTS['SPADES'] = strategy_1_spades.out[1]
      }
      // Unicycler
      if (!params.skip_unicycler) {
        strategy_1_unicycler(input_branches.main)
        HYBRID_OUTPUTS['UNICYCLER'] = strategy_1_unicycler.out[1]
      }
      // Haslr
      if (!params.skip_haslr) {
        strategy_1_haslr(input_branches.main)
        HYBRID_OUTPUTS['HASLR'] = strategy_1_haslr.out[1]
      }

      // Get hybrid assemblies
      HYBRID_OUTPUTS['ASSEMBLIES'] = HYBRID_OUTPUTS['SPADES']
                                     .mix(HYBRID_OUTPUTS['UNICYCLER'], 
                                          HYBRID_OUTPUTS['HASLR'])
                                     .combine(input_tuple, by: 0)

      /*
       * Polish a long reads assembly
       */
      
      /*
       * Canu
       */
      if (!params.skip_canu) {
        strategy_2_canu(input_branches.secondary)
        LONGREADS_OUTPUTS['CANU'] = strategy_2_canu.out[1]
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        strategy_2_flye(input_branches.secondary)
        LONGREADS_OUTPUTS['FLYE'] = strategy_2_flye.out[1]
      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        strategy_2_unicycler(input_branches.secondary)
        LONGREADS_OUTPUTS['UNICYCLER'] = strategy_2_unicycler.out[1]
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        strategy_2_raven(input_branches.secondary)
        LONGREADS_OUTPUTS['RAVEN'] = strategy_2_raven.out[1]
      }

      /*
       * Shasta
       */
      if (!params.skip_shasta) {
        strategy_2_shasta(input_branches.secondary)
        LONGREADS_OUTPUTS['SHASTA'] = strategy_2_shasta.out[1]
      }

      /*
       * wtdbg2
       */
      if (!params.skip_wtdbg2) {
        strategy_2_wtdbg2(input_branches.secondary)
        LONGREADS_OUTPUTS['WTDBG2'] = strategy_2_wtdbg2.out[1]
      }

      // Get long reads assemblies
      LONGREADS_OUTPUTS['RAW_ASSEMBLIES'] = LONGREADS_OUTPUTS['CANU']
                                            .mix(LONGREADS_OUTPUTS['FLYE'], 
                                                 LONGREADS_OUTPUTS['UNICYCLER'], 
                                                 LONGREADS_OUTPUTS['RAVEN'], 
                                                 LONGREADS_OUTPUTS['WTDBG2'], 
                                                 LONGREADS_OUTPUTS['SHASTA'])
                                            .combine(input_tuple, by: 0)

      /*
       * Run medaka?
       */
      strategy_2_medaka(LONGREADS_OUTPUTS['RAW_ASSEMBLIES'])
      LONGREADS_OUTPUTS['MEDAKA'] = strategy_2_medaka.out[1]

      /*
       * Run nanopolish?
       */
      strategy_2_nanopolish(LONGREADS_OUTPUTS['RAW_ASSEMBLIES'])
      LONGREADS_OUTPUTS['NANOPOLISH'] = strategy_2_nanopolish.out[0]

      /*
       * gcpp?
       */
      strategy_2_gcpp(LONGREADS_OUTPUTS['RAW_ASSEMBLIES'])
      LONGREADS_OUTPUTS['GCPP'] = strategy_2_gcpp.out[1]

      // Gather long reads assemblies polished
      LONGREADS_OUTPUTS['POLISHED_ASSEMBLIES'] = LONGREADS_OUTPUTS['MEDAKA']
                                                 .mix(LONGREADS_OUTPUTS['NANOPOLISH'],
                                                      LONGREADS_OUTPUTS['GCPP'])
                                                 .combine(input_tuple, by: 0)

      /*
       * Finally, run pilon for all
       */
      
      if (params.skip_raw_assemblies_polishing) {
        ch_pilon = LONGREADS_OUTPUTS['POLISHED_ASSEMBLIES']
      } else {
        ch_pilon = LONGREADS_OUTPUTS['RAW_ASSEMBLIES'].mix(LONGREADS_OUTPUTS['POLISHED_ASSEMBLIES'])
      }
      strategy_2_pilon(ch_pilon)
      HYBRID_OUTPUTS['PILON'] = strategy_2_pilon.out[1].combine(input_tuple, by: 0)

      // Gather assemblies for qc
      HYBRID_OUTPUTS['ALL_RESULTS'] = HYBRID_OUTPUTS['ASSEMBLIES']
                                      .mix(LONGREADS_OUTPUTS['RAW_ASSEMBLIES'],
                                           LONGREADS_OUTPUTS['POLISHED_ASSEMBLIES'],
                                           HYBRID_OUTPUTS['PILON'])
  
  emit:
    HYBRID_OUTPUTS['ALL_RESULTS']

}
