#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import org.yaml.snakeyaml.Yaml

/*
 * Generic multiplatform genome assembly pipeline (MpGAP)
 */

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/
WorkflowMpGAP.initialise(params, log)
WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    LOAD WORKFLOWS
========================================================================================
*/

include { PARSE_SAMPLESHEET } from './workflows/parse_samples.nf'
include { ASSEMBLY_QC       } from './workflows/assembly_qc.nf'
include { SHORTREADS_ONLY   } from './workflows/short-reads-only.nf'
include { LONGREADS_ONLY    } from './workflows/long-reads-only.nf'
include { HYBRID            } from './workflows/hybrid.nf'

/*
========================================================================================
    DEFINE MAIN WORKFLOW
========================================================================================
*/

workflow {

  // Message to user
  println("""
    Launching defined workflows!
    By default, all workflows will appear in the console "log" message.
    However, the processes of each workflow will be launched based on the inputs received.
    You can see that processes that were not launched have an empty [-       ].
  """)

  // Load YAML
  samplesheet_yaml = file(params.input)
  parameter_yaml   = samplesheet_yaml.readLines().join("\n")
  new Yaml().load(parameter_yaml).each { k, v -> params[k] = v }

  // Copy YAML samplesheet to output directory so user has a copy of it
  file(params.output).mkdir()
  samplesheet_yaml.copyTo(params.output + "/" + "${samplesheet_yaml.getName()}")

  // ch for versions
  ch_versions = Channel.empty()

  // Parse YAML file
  PARSE_SAMPLESHEET( params.samplesheet )

  // short reads only samples
  SHORTREADS_ONLY( PARSE_SAMPLESHEET.out.shortreads )
  ch_versions = ch_versions.mix(SHORTREADS_ONLY.out.versions)
    
  // long reads only samples
  LONGREADS_ONLY( PARSE_SAMPLESHEET.out.longreads )
  ch_versions = ch_versions.mix(LONGREADS_ONLY.out.versions)

  // hybrid samples
  HYBRID( PARSE_SAMPLESHEET.out.hybrid )
  ch_versions = ch_versions.mix(HYBRID.out.versions)

  // QC
  ASSEMBLY_QC( 
    SHORTREADS_ONLY.out.results.mix( LONGREADS_ONLY.out.results, HYBRID.out.results ),
    ch_versions
  )
    
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
