/*
 * DEFINITION OF MULTI-SAMPLE (BATCH) MODE
 */

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

// Hifiasm assembler
include { hifiasm } from '../modules/LongReads/hifiasm.nf'

/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include { nanopolish } from '../modules/LongReads/nanopolish.nf'

// Medaka (for nanopore data)
include { medaka } from '../modules/LongReads/medaka.nf'

// gcpp Pacbio
include { gcpp } from '../modules/LongReads/gcpp.nf'

workflow LONGREADS_ONLY {
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
      LONGREADS_OUTPUTS['HIFIASM']     = Channel.empty()

      /*
       * Canu
       */
      if (!params.skip_canu) {
        canu(input_tuple)
        LONGREADS_OUTPUTS['CANU']  = canu.out[1]
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        flye(input_tuple)
        LONGREADS_OUTPUTS['FLYE']  = flye.out[1]
      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        unicycler(input_tuple)
        LONGREADS_OUTPUTS['UNICYCLER'] = unicycler.out[1]
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        raven(input_tuple)
        LONGREADS_OUTPUTS['RAVEN'] = raven.out[1]
      }

      /*
       * Shasta
       */
      if (!params.skip_shasta) {
        shasta(input_tuple)
        LONGREADS_OUTPUTS['SHASTA'] = shasta.out[1]
      }

      /*
       * wtdbg2
       */
      if (!params.skip_wtdbg2) {
        wtdbg2(input_tuple)
        LONGREADS_OUTPUTS['WTDBG2'] = wtdbg2.out[1]
      }

      /*
       *Hifiasm
       */
      if (!params.skip_hifiasm) {
        hifiasm(input_tuple)
        LONGREADS_OUTPUTS['HIFIASM'] = hifiasm.out[1]
      }

      // gather assemblies for polishing steps
      LONGREADS_OUTPUTS['RAW_ASSEMBLIES'] = LONGREADS_OUTPUTS['CANU']
                                            .mix(LONGREADS_OUTPUTS['UNICYCLER'], 
                                                 LONGREADS_OUTPUTS['FLYE'] , 
                                                 LONGREADS_OUTPUTS['RAVEN'], 
                                                 LONGREADS_OUTPUTS['WTDBG2'], 
                                                 LONGREADS_OUTPUTS['SHASTA'],
                                                 LONGREADS_OUTPUTS['HIFIASM'])
                                            .combine(input_tuple, by: 0)

      /*
       * Run medaka?
       */
      medaka(LONGREADS_OUTPUTS['RAW_ASSEMBLIES'])
      LONGREADS_OUTPUTS['MEDAKA'] = medaka.out[1]

      /*
       * Run nanopolish?
       */
      nanopolish(LONGREADS_OUTPUTS['RAW_ASSEMBLIES'])
      LONGREADS_OUTPUTS['NANOPOLISH'] = nanopolish.out[0]

      /*
       * gcpp?
       */
      gcpp(LONGREADS_OUTPUTS['RAW_ASSEMBLIES'])
      LONGREADS_OUTPUTS['GCPP'] = gcpp.out[1]

      // Gather polished assemblies
      LONGREADS_OUTPUTS['POLISHED_ASSEMBLIES'] = LONGREADS_OUTPUTS['MEDAKA']
                                                 .mix(LONGREADS_OUTPUTS['NANOPOLISH'], 
                                                      LONGREADS_OUTPUTS['GCPP'])
                                                  .combine(input_tuple, by: 0)

      // Gather all assemblies for qc
      LONGREADS_OUTPUTS['ALL_RESULTS'] = LONGREADS_OUTPUTS['RAW_ASSEMBLIES']
                                         .mix(LONGREADS_OUTPUTS['POLISHED_ASSEMBLIES'])
  
  emit:
    LONGREADS_OUTPUTS['ALL_RESULTS']
}
