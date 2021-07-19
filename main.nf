#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Generic multiplatform genome assembly pipeline (MpGAP)
 */

/*
 * Include functions
 */
include { helpMessage } from './nf_functions/help.nf'
include { exampleMessage } from './nf_functions/examples.nf'
include { paramsCheck } from './nf_functions/paramsCheck.nf'
include { logMessage } from './nf_functions/logMessages.nf'

/*
 * Check parameters
 */
paramsCheck()
params.help = false
if (params.help){
  helpMessage()
  exit 0
}
params.examples = false
if (params.examples){
  exampleMessage()
  exit 0
}

 /*
           Download configuration file, if necessary.
 */
 params.get_hybrid_config = false
 if (params.get_hybrid_config) {
   new File("hybrid.config").write(new URL ("https://github.com/fmalmeida/mpgap/raw/master/configuration_example/hybrid.config").getText())
   println ""
   println "hybrid.config file saved in working directory"
   println "After configuration, run:"
   println "nextflow run fmalmeida/mpgap -c ./hybrid.config"
   println "Nice code!\n"

   exit 0
 }
 params.get_lreads_config = false
 if (params.get_lreads_config) {
   new File("lreads-only.config").write(new URL ("https://github.com/fmalmeida/mpgap/raw/master/configuration_example/lreads.config").getText())
   println ""
   println "lreads.config file saved in working directory"
   println "After configuration, run:"
   println "nextflow run fmalmeida/mpgap -c ./lreads.config"
   println "Nice code!\n"

   exit 0
 }
 params.get_sreads_config = false
 if (params.get_sreads_config) {
   new File("sreads-only.config").write(new URL ("https://github.com/fmalmeida/mpgap/raw/master/configuration_example/sreads.config").getText())
   println ""
   println "sreads.config file saved in working directory"
   println "After configuration, run:"
   println "nextflow run fmalmeida/mpgap -c ./sreads.config"
   println "Nice code!\n"

   exit 0
 }

 /*
  * Load general parameters and establish defaults
  */

// General
params.outdir  = 'output'
params.threads = 4
params.cpus    = 2

// Assemblers?
params.skip_flye      = false
params.skip_spades    = false
params.skip_shovill   = false
params.skip_canu      = false
params.skip_unicycler = false
params.skip_haslr     = false
params.skip_raven     = false

// Additional parameters for assemblers
params.genomeSize = ''
params.canu_additional_parameters      = ''
params.unicycler_additional_parameters = ''
params.flye_additional_parameters      = ''
params.spades_additional_parameters    = ''
params.shovill_additional_parameters   = ''
params.haslr_additional_parameters     = ''
params.raven_additional_parameters     = ''

// Short reads
params.shortreads_paired = ''
params.shortreads_single = ''

// Long reads
params.corrected_lreads = false
params.longreads = ''
params.lr_type = ''
params.medaka_sequencing_model = ''
params.nanopolish_fast5Path = ''
params.nanopolish_max_haplotypes = 1000
params.pacbio_bams = ''

// Hybrid strategy 2
params.strategy_2 = false
params.pilon_memory_limit = 50

/*
 * Define log message
 */
logMessage()

/*
 * Define custom workflows
 */

// Short reads only
include { sreads_only_nf } from './workflows/short-reads-only.nf'

// Long reads only
include { lreadsonly_nf } from './workflows/long-reads-only.nf'

// Hybrid
include { hybrid_nf } from './workflows/hybrid.nf'

                                  /*
                                   * DEFINE (RUN) MAIN WORKFLOW
                                   */

workflow {

  /*
   * Long reads only assembly
   */

  if (!params.shortreads_paired && !params.shortreads_single &&
      params.longreads && params.lr_type) {

    // Giving inputs
    lreadsonly_nf(
      // Longreads - required
      Channel.fromPath(params.longreads),

      // Will run Nanopolish?
      (params.nanopolish_fast5Path && params.lr_type == 'nanopore') ? Channel.fromPath(params.nanopolish_fast5Path) : Channel.value(''),
      (params.nanopolish_fast5Path && params.lr_type == 'nanopore') ? Channel.fromPath(params.nanopolish_fast5Path, type: 'dir') : Channel.value(''),

      // Will run gcpp?
      (params.pacbio_bams && params.lr_type == 'pacbio') ? Channel.fromPath(params.pacbio_bams).collect() : Channel.value(''),
      (params.pacbio_bams && params.lr_type == 'pacbio') ? Channel.fromPath(params.pacbio_bams).count().subscribe { println it } : Channel.value('')
    )
  }


  /*
   * Short reads only assembly
   */

   if (!params.longreads &&
       (params.shortreads_paired || params.shortreads_single)) {

     // Giving inputs
     sreads_only_nf(
       // Have paired end reads?
       (params.shortreads_paired) ? Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ) : Channel.value(['', '', '']),

       // Have unpaired reads?
       (params.shortreads_single) ? Channel.fromPath(params.shortreads_single, hidden: true) : Channel.value('')
     )
   }

  /*
   * Hybrid assembly
   */

   if (params.longreads && params.lr_type &&
       (params.shortreads_paired || params.shortreads_single)) {

     // Giving inputs
     hybrid_nf(
      // Have paired end reads?
      (params.shortreads_paired) ? Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ) : Channel.value(['', '', '']),

      // Have unpaired reads?
      (params.shortreads_single) ? Channel.fromPath(params.shortreads_single, hidden: true) : Channel.value(''),

      // Long reads - required
      Channel.fromPath(params.longreads),

      // Will run Nanopolish?
      (params.nanopolish_fast5Path && params.lr_type == 'nanopore') ? Channel.fromPath(params.nanopolish_fast5Path) : Channel.value(''),
      (params.nanopolish_fast5Path && params.lr_type == 'nanopore') ? Channel.fromPath(params.nanopolish_fast5Path, type: 'dir') : Channel.value(''),

      // Will run Arrow?
      (params.pacbio_bams && params.lr_type == 'pacbio') ? Channel.fromPath(params.pacbio_bams).collect() : Channel.value(''),
      (params.pacbio_bams && params.lr_type == 'pacbio') ? Channel.fromPath(params.pacbio_bams).count().subscribe { println it } : Channel.value('')
     )
   }

}

/*
 * Completition message
 */
 workflow.onComplete {
     println "Pipeline completed at: $workflow.complete"
     println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
     println "Execution duration: $workflow.duration"
     println ""
     println "${ workflow.success ? 'I wish you nice results!' : 'Do not give up, we can fix it!' }"
     println "${ workflow.success ? 'Thank you for using fmalmeida/mpgap pipeline!' : '' }"
     println ""
 }
