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

  // Parse YAML file
  PARSE_SAMPLESHEET( params.samplesheet )

  // short reads only samples
  SHORTREADS_ONLY( PARSE_SAMPLESHEET.out.shortreads )
    
  // long reads only samples
  LONGREADS_ONLY( PARSE_SAMPLESHEET.out.longreads )

  // hybrid samples
  HYBRID( PARSE_SAMPLESHEET.out.hybrid )

  // QC
  ch_all_assemblies = SHORTREADS_ONLY.out.mix( LONGREADS_ONLY.out, HYBRID.out )
  ASSEMBLY_QC( ch_all_assemblies )

  // generate bacannot samplesheet
  def final_outdir = file(params.output).toUriString()
  Channel.value( 'samplesheet:' )
  .mix( 
    ch_all_assemblies
    .map{ 
      def sample   = it[0].toString()
      def asm_type = it[2].toString()
      def assembly = it[1].toString().split('/')[-1]
      def asm_path = "${final_outdir}/final_assemblies/${sample}_${assembly}"

      def final_string = "\s\s- id: ${sample}_${asm_type}\n\s\s\s\sassembly: ${asm_path}\n"
    },

    Channel.value("\n")
  )
  .collectFile( name: 'bacannot_samplesheet.yml', storeDir: params.output, sort: false, newLine: false )
    
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
