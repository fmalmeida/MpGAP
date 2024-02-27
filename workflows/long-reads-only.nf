/*
 * DEFINITION OF MULTI-SAMPLE (BATCH) MODE
 */

/*
 * Modules for assembling long reads
 */

// Canu assembler
include { canu } from '../modules/local/LongReads/canu.nf'

// Unicycler assembler
include { unicycler } from '../modules/local/LongReads/unicycler.nf'

// Flye assembler
include { flye } from '../modules/local/LongReads/flye.nf'

// Raven assembler
include { raven } from '../modules/local/LongReads/raven.nf'

// wtdbg2 assembler
include { wtdbg2 } from '../modules/local/LongReads/wtdbg2.nf'

// Shasta assembler
include { shasta } from '../modules/local/LongReads/shasta.nf'

// Hifiasm assembler
include { hifiasm } from '../modules/LongReads/hifiasm.nf'

/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include { nanopolish } from '../modules/local/LongReads/nanopolish.nf'

// Medaka (for nanopore data)
include { medaka } from '../modules/local/LongReads/medaka.nf'

// gcpp Pacbio
include { gcpp } from '../modules/local/LongReads/gcpp.nf'

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

      ch_versions_lr = Channel.empty()

      /*
       * Canu
       */
      if (!params.skip_canu) {
        canu(input_tuple)
        LONGREADS_OUTPUTS['CANU']  = canu.out[1]
        ch_versions_lr = ch_versions_lr.mix(canu.out.versions.first())
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        flye(input_tuple)
        LONGREADS_OUTPUTS['FLYE']  = flye.out[1]
        ch_versions_lr = ch_versions_lr.mix(flye.out.versions.first())
      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        unicycler(input_tuple)
        LONGREADS_OUTPUTS['UNICYCLER'] = unicycler.out[1]
        ch_versions_lr = ch_versions_lr.mix(unicycler.out.versions.first())
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        raven(input_tuple)
        LONGREADS_OUTPUTS['RAVEN'] = raven.out[1]
        ch_versions_lr = ch_versions_lr.mix(raven.out.versions.first())
      }

      /*
       * Shasta
       */
      if (!params.skip_shasta) {
        shasta(input_tuple)
        LONGREADS_OUTPUTS['SHASTA'] = shasta.out[1]
        ch_versions_lr = ch_versions_lr.mix(shasta.out.versions.first())
      }

      /*
       * wtdbg2
       */
      if (!params.skip_wtdbg2) {
        wtdbg2(input_tuple)
        LONGREADS_OUTPUTS['WTDBG2'] = wtdbg2.out[1]
        ch_versions_lr = ch_versions_lr.mix(wtdbg2.out.versions.first())
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
      ch_versions_lr = ch_versions_lr.mix(medaka.out.versions.first())

      /*
       * Run nanopolish?
       */
      nanopolish(LONGREADS_OUTPUTS['RAW_ASSEMBLIES'])
      LONGREADS_OUTPUTS['NANOPOLISH'] = nanopolish.out[0]
      ch_versions_lr = ch_versions_lr.mix(nanopolish.out.versions.first())

      /*
       * gcpp?
       */
      gcpp(LONGREADS_OUTPUTS['RAW_ASSEMBLIES'])
      LONGREADS_OUTPUTS['GCPP'] = gcpp.out[1]
      ch_versions_lr = ch_versions_lr.mix(gcpp.out.versions.first())

      // Gather polished assemblies
      LONGREADS_OUTPUTS['POLISHED_ASSEMBLIES'] = LONGREADS_OUTPUTS['MEDAKA']
                                                 .mix(LONGREADS_OUTPUTS['NANOPOLISH'], 
                                                      LONGREADS_OUTPUTS['GCPP'])
                                                  .combine(input_tuple, by: 0)

      // Gather all assemblies for qc
      LONGREADS_OUTPUTS['ALL_RESULTS'] = LONGREADS_OUTPUTS['RAW_ASSEMBLIES']
                                         .mix(LONGREADS_OUTPUTS['POLISHED_ASSEMBLIES'])
  
  emit:
  results  = LONGREADS_OUTPUTS['ALL_RESULTS']
  versions = ch_versions_lr
}
