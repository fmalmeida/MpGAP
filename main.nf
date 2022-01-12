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
include { logMessage          } from './nf_functions/logMessages.nf'

 /*
  * Help message
  */
params.help = false
if (params.help) {
  helpMessage
  exit 1
}

 /*
           Download configuration file, if necessary.
 */
 params.get_config = false
 if (params.get_config) {
   new File("MPGAP.config").write(new URL ("https://github.com/fmalmeida/mpgap/raw/master/nextflow.config").getText())
   println ""
   println "MPGAP.config file saved in working directory"
   println "After configuration, run:"
   println "nextflow run fmalmeida/mpgap -c ./MPGAP.config"
   println "Nice code!\n"

   exit 0
 }

 /*
           Download samplesheet, if necessary.
 */
 params.get_samplesheet = false
 if (params.get_samplesheet) {
   new File("MPGAP_samplesheet.yml").write(new URL ("https://github.com/fmalmeida/mpgap/raw/master/example_samplesheet.yml").getText())
   println ""
   println "Samplesheet (MPGAP_samplesheet.yml) file saved in working directory"
   println "Nice code!\n"

   exit 0
 }

 /*
  * Load general parameters and establish defaults
  */

// General
params.output  = 'output'
params.threads = 3
params.input   = ''

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

// shasta default configuration
params.shasta_config = 'Nanopore-Oct2021'

// medaka model default configuration
params.medaka_model = 'r941_min_high_g360'

// genome size for canu and wtdbg2 and haslr
params.genome_size = ''

// Additional parameters for assemblers and quast
params.quast_additional_parameters     = ''
params.canu_additional_parameters      = ''
params.unicycler_additional_parameters = ''
params.flye_additional_parameters      = ''
params.spades_additional_parameters    = ''
params.shovill_additional_parameters   = ''
params.haslr_additional_parameters     = ''
params.raven_additional_parameters     = ''
params.wtdbg2_additional_parameters    = ''
params.shasta_additional_parameters    = ''

// Long reads
params.corrected_long_reads = false
params.nanopolish_max_haplotypes = 1000

// Hybrid strategies default configuration
params.hybrid_strategy = 1
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
include { ASSEMBLY_QC } from './workflows/assembly_qc.nf'

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
  if (params.input) {

    // Load YAML
    samplesheet_yaml = file(params.input)
    parameter_yaml   = samplesheet_yaml.readLines().join("\n")
    new Yaml().load(parameter_yaml).each { k, v -> params[k] = v }

    // Copy YAML samplesheet to output directory so user has a copy of it
    file(params.output).mkdir()
    samplesheet_yaml.copyTo(params.output + "/" + "${samplesheet_yaml.getName()}")

    // Parse YAML file
    parse_samplesheet(params.samplesheet)

    // short reads only samples
    SHORTREADS_ONLY(parse_samplesheet.out[0])
    
    // long reads only samples
    LONGREADS_ONLY(parse_samplesheet.out[1])

    // hybrid samples
    HYBRID(parse_samplesheet.out[2])

    // QC
    ASSEMBLY_QC(SHORTREADS_ONLY.out.mix(LONGREADS_ONLY.out, HYBRID.out))

  } else {

    // Message to user
    println("""
    ERROR!
    A major error has occurred!
      ==> A samplesheet has not been provided. Please, provide a samplesheet to run the analysis. Online documentation is available at: https://mpgap.readthedocs.io/en/latest/
    Please, read the docs.
    Cheers.
    """)
  
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
