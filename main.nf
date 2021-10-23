#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import org.yaml.snakeyaml.Yaml

/*
 * Generic multiplatform genome assembly pipeline (MpGAP)
 */

/*
 * Include functions
 */
include { helpMessage         } from './nf_functions/help.nf'
include { helpMessageAdvanced } from './nf_functions/help.nf'
include { exampleMessage      } from './nf_functions/examples.nf'
include { paramsCheck         } from './nf_functions/paramsCheck.nf'
include { logMessage          } from './nf_functions/logMessages.nf'

/*
 * Check parameters
 */
paramsCheck()
params.help = false
if (params.help){
  helpMessage()
  exit 0
}
params.show_advanced_parameters = false
if (params.show_advanced_parameters){
  helpMessageAdvanced()
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
params.prefix  = ''
params.threads = 4
params.cpus    = 2
params.in_yaml = ''

// Assemblers?
params.skip_flye      = false
params.skip_spades    = false
params.skip_shovill   = false
params.skip_canu      = false
params.skip_unicycler = false
params.skip_haslr     = false
params.skip_raven     = false
params.skip_wtdbg2    = false
params.skip_shasta    = false

// Additional parameters for assemblers and quast
params.genomeSize = ''
params.quast_additional_parameters     = ''
params.canu_additional_parameters      = ''
params.unicycler_additional_parameters = ''
params.flye_additional_parameters      = ''
params.spades_additional_parameters    = ''
params.shovill_additional_parameters   = ''
params.haslr_additional_parameters     = ''
params.raven_additional_parameters     = ''
params.wtdbg2_additional_parameters    = ''
params.wtdbg2_technology               = 'ont'
params.shasta_additional_parameters    = ''

// Short reads
params.shortreads_paired = ''
params.shortreads_single = ''

// Long reads
params.corrected_lreads = false
params.longreads = ''
params.lr_type = 'nanopore'
params.medaka_sequencing_model = 'r941_min_high_g360'
params.nanopolish_fast5Path = ''
params.nanopolish_max_haplotypes = 1000
params.pacbio_bam = ''

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

// misc
include { parse_samplesheet } from './workflows/parse_samples.nf'
include { define_prefix } from './modules/misc/define_prefix.nf'

// Short reads only
include { SHORTREADS_ONLY } from './workflows/short-reads-only.nf'

// Long reads only
include { LONGREADS_ONLY } from './workflows/long-reads-only.nf'

// Hybrid
include { HYBRID } from './workflows/hybrid.nf'

                                  /*
                                   * DEFINE (RUN) MAIN WORKFLOW
                                   */

workflow {

  // Message to user
  println("""
    Launching defined workflows!
    By default, all workflows will appear in the console "log" message.
    However, the processes of each workflow will be launched based on the inputs received.
    You can see that processes that were not launched have an empty [-       ].
  """)

  // with samplesheet?
  if (params.in_yaml) {

     // MULTI SAMPLE ANALYSIS --- WITH SAMPLESHEET

    parameter_yaml = new FileInputStream(new File(params.in_yaml))
    new Yaml().load(parameter_yaml).each { k, v -> params[k] = v }

    // Read YAML file
    parse_samplesheet(params.samplesheet)

    // short reads only samples
    SHORTREADS_ONLY(parse_samplesheet.out[0])
    
    // long reads only samples
    LONGREADS_ONLY(parse_samplesheet.out[1])

    // hybrid samples
    HYBRID(parse_samplesheet.out[2])

  } else {

    // SINGLE SAMPLE ANALYSIS --- WITHOUT SAMPLESHEET

  /*
   * Parse inputs
   */
    // load long reads
    lreads = (params.longreads) ? Channel.fromPath(params.longreads) : Channel.value('missing')

    // load short reads paired
    sreads_paired = (params.shortreads_paired) ? Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ) : Channel.value(['missing', 'missing', 'missing'])

    // separate paired end reads to create input_tuple
    sreads_paired.multiMap {
      id:  it[0]
      fwd: it[1]
      rev: it[2]
    }.set { parsed_paired }

    // load short reads unpaired
    sreads_single = (params.shortreads_single) ? Channel.fromPath(params.shortreads_single, hidden: true) : Channel.value('missing')

    // define prefixes
    define_prefix(lreads, sreads_paired, sreads_single)

    // check desired workflow
    if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
      desired_workflow = "longreads_only"
    } else if (!params.longreads && (params.shortreads_paired || params.shortreads_single)) {
      desired_workflow = "shortreads_only"
    } else if (params.longreads && params.lr_type && (params.shortreads_paired || params.shortreads_single)) {
      desired_workflow = (params.strategy_2) ? "hybrid_strategy_2" : "hybrid_strategy_1"
    }

    // create input tuple
    // for concat, every "entry" must be a channel
    input_tuple = define_prefix.out.sample_name.concat(
      Channel.from(desired_workflow),
      (params.shortreads_paired) ? parsed_paired.fwd : Channel.from("missing_pairFWD"),
      (params.shortreads_paired) ? parsed_paired.rev : Channel.from("missing_pairREV"),
      (params.shortreads_single) ? sreads_single : Channel.from("missing_single"),
      (params.longreads) ? lreads : Channel.from("missing_lreads"),
      Channel.from(
        params.lr_type,
        params.wtdbg2_technology,
        params.genomeSize,
        (params.corrected_lreads) ? 'true' : 'false',
        params.medaka_sequencing_model
      ),
      (params.nanopolish_fast5Path && params.lr_type == 'nanopore') ? Channel.fromPath(params.nanopolish_fast5Path) : Channel.from("missing_fast5"),
      (params.pacbio_bam && params.lr_type == 'pacbio') ? Channel.fromPath(params.pacbio_bam).toList() : Channel.from("missing_pacbio_bam"),
      define_prefix.out.prefix_dir
    ).toList()

    // long reads only workflow
    if (desired_workflow == "longreads_only") {
      LONGREADS_ONLY(input_tuple)
    }

    // short reads only workflow
    if (desired_workflow == "shortreads_only") {
      SHORTREADS_ONLY(input_tuple)
    }

    // hybrid workflow
    if (desired_workflow =~ /hybrid/) {
      HYBRID(input_tuple)
    }

  } // end of else statement -- single genome workflow

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
