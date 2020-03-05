#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
 * Generic multiplatform genome assembly pipeline (MpGAP)
 */

/*
 * Define help message
 */
 def helpMessage() {
    log.info """
    Usage:
    nextflow run fmalmeida/MpGAP [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]

    Comments:
    This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would
    cause the command to be huge. Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make
    parameterization easier and more readable.

    Creating a configuration file:
    nextflow run fmalmeida/MpGAP [--get_hybrid_config] [--get_lreads_config] [--get_sreads_config]

    Show command line examples:
    nextflow run fmalmeida/MpGAP --show

    Execution Reports:
    nextflow run fmalmeida/MpGAP [ -c nextflow.config ] -with-report
    nextflow run fmalmeida/MpGAP [ -c nextflow.config ] -with-trace
    nextflow run fmalmeida/MpGAP [ -c nextflow.config ] -with-timeline

    OBS: These reports can also be enabled through the configuration file.

    OPTIONS:
             General Parameters - Mandatory

     --outdir <string>                      Output directory name
     --prefix <string>                      Set prefix for output files
     --threads <int>                        Number of threads to use
     --assembly_type <string>               Selects assembly mode: hybrid, illumina-only or longreads-only
     --try_canu                             Execute assembly with Canu. Multiple assemblers can be chosen.
     --canu_additional_parameters           Give additional parameters to Canu assembler. Must be in quotes
                                            and separated by one space. Must be given as shown in Canu manual.
                                            E.g. 'correctedErrorRate=0.075 corOutCoverage=200'
     --try_unicycler                        Execute assembly with Unicycler. Multiple assemblers can be chosen.
     --unicycler_additional_parameters      Give additional parameters to Unicycler assembler. Must be in quotes
                                            and separated by one space. Must be given as shown in Unicycler manual.
                                            E.g. '--mode conservative --no_correct'
     --try_flye                             Execute assembly with Flye. Multiple assemblers can be chosen.
     --flye_additional_parameters           Give additional parameters to Flye assembler. Must be in quotes
                                            and separated by one space. Must be given as shown in Unicycler manual.
                                            E.g. '--meta --iterations 4'
     --try_spades                           Execute assembly with Spades. Multiple assemblers can be chosen.
     --spades_additional_parameters         Give additional parameters to Spades assembler. Must be in quotes
                                            and separated by one space. Must be given as shown in Spades manual.
                                            E.g. '--meta --plasmids'


             Parameters for illumina-only mode. Can be executed by SPAdes and Unicycler assemblers.
             Users can use paired or single end reads. If both types are given at once, assemblers
             will be executed with a mix of both.

     --shortreads_paired <string>           Path to Illumina paired end reads.
     --shortreads_single <string>           Path to Illumina single end reads.
     --ref_genome <string>                  Path to reference genome for guided assembly. Used only by SPAdes.

             Parameters for hybrid mode. Can be executed by spades and unicycler assemblers.

     --shortreads_paired <string>           Path to Illumina paired end reads.
     --shortreads_single <string>           Path to Illumina single end reads.
     --ref_genome <string>                  Path to reference genome for guided assembly. Used only by SPAdes.
     --longreads <string>                   Path to longreads in FASTA or FASTQ formats.
     --lr_type <string>                     Sets wich type of long reads are being used: pacbio or nanopore

             Parameters for longreads-only mode. Can be executed by canu, flye and unicycler assemblers.
             In the end, long reads only assemblies can be polished with illumina reads through pilon.

     --longreads <string>                   Path to longreads in FASTA or FASTQ formats.
     --fast5Path <string>                   Path to directory containing FAST5 files for given reads.
                                            Whenever set, the pipeline will execute a polishing step
                                            with Nanopolish. This makes the pipeline extremely SLOW!!
     --pacbio_all_baxh5_path <string>       Path to all bax.h5 files for given reads. Whenever set, the pipeline
                                            will execute a polishing step with VarianCaller.
     --pacbio_all_bam_path <string>         Path to all subreads bam files for given reads. Whenever set, the pipeline
                                            will execute a polishing step with VarianCaller.
     --genomeSize                           Canu and Flye require an estimative of genome size in order
                                            to be executed. Examples: 5.6m; 1.2g
     --lr_type <string>                     Sets wich type of long reads are being used: pacbio or nanopore
     --illumina_polish_longreads_contigs    This tells the pipeline to polish long reads only assemblies
                                            with Illumina reads through Pilon. This is another hybrid methodology.
                                            For that, users have to set path to Illumina reads through
                                            --shortreads_paired or --shortreads_single.



    """.stripIndent()
 }

 def exampleMessage() {
    log.info """
    """.stripIndent()
 }

 /*
           Display Help Message
 */
 params.help = false
  // Show help emssage
  if (params.help){
    helpMessage()
    //file('work').deleteDir()
    exit 0
 }

 /*
           Display CLI examples
 */
 params.show = false
  // Show help emssage
  if (params.show){
    exampleMessage()
    exit 0
 }

 /*
           Download configuration file, if necessary.
 */
 params.get_hybrid_config = false
 if (params.get_hybrid_config) {
   new File("hybrid.config") << new URL ("https://github.com/fmalmeida/MpGAP/raw/master/configuration_example/hybrid.config").getText()
   println ""
   println "hybrid.config file saved in working directory"
   println "After configuration, run:"
   println "nextflow run fmalmeida/MpGAP -c ./hybrid.config"
   println "Nice code!\n"

   exit 0
 }
 params.get_lreads_config = false
 if (params.get_lreads_config) {
   new File("lreads-only.config") << new URL ("https://github.com/fmalmeida/MpGAP/raw/master/configuration_example/lreads.config").getText()
   println ""
   println "lreads.config file saved in working directory"
   println "After configuration, run:"
   println "nextflow run fmalmeida/MpGAP -c ./lreads.config"
   println "Nice code!\n"

   exit 0
 }
 params.get_sreads_config = false
 if (params.get_sreads_config) {
   new File("sreads-only.config") << new URL ("https://github.com/fmalmeida/MpGAP/raw/master/configuration_example/sreads.config").getText()
   println ""
   println "sreads.config file saved in working directory"
   println "After configuration, run:"
   println "nextflow run fmalmeida/MpGAP -c ./sreads.config"
   println "Nice code!\n"

   exit 0
 }

 /*
  * Load general parameters and establish defaults
  */
  params.longreads = ''
  params.fast5Path = ''
  params.pacbio_all_baxh5_path = ''
  params.pacbio_all_bam_path = ''
  params.lr_type = ''
  params.shortreads_paired = ''
  params.shortreads_single = ''
  params.ref_genome = ''
  params.assembly_type = ''
  params.illumina_polish_longreads_contigs = false
  params.pilon_memmory_limit = 50
  params.try_canu = false
  params.canu_additional_parameters = ''
  params.try_unicycler = false
  params.unicycler_additional_parameters = ''
  params.try_flye = false
  params.flye_additional_parameters = ''
  params.try_spades = false
  params.spades_additional_parameters = ''
  params.genomeSize = ''
  params.outdir = 'output'
  params.prefix = 'out'
  params.threads = 3
  params.cpus = 2

/*
 * Define log message
 */
log.info "================================================================="
log.info " Docker-based, fmalmeida/mpgap, generic genome assembly Pipeline "
log.info "================================================================="
def summary = [:]
summary['Long Reads']   = params.longreads
summary['Fast5 files dir']   = params.fast5Path
summary['Long Reads']   = params.longreads
summary['Short single end reads']   = params.shortreads_single
summary['Short paired end reads']   = params.shortreads_paired
summary['Fasta Ref']    = params.ref_genome
summary['Output dir']   = params.outdir
summary['Assembly assembly_type chosen'] = params.assembly_type
summary['Long read sequencing technology'] = params.lr_type
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Command used']   = "$workflow.commandLine"
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 * Include modules
 */

/*
 * Modules for assembling long reads
 */

// Canu assembler
include canu_assembly from './modules/canu.nf' params(outdir: params.outdir, lr_type: params.lr_type,
  canu_additional_parameters: params.canu_additional_parameters, threads: params.threads,
  genomeSize: params.genomeSize, prefix: params.prefix)

// Unicycler assembler
include unicycler_lreads from './modules/unicycler_lreads.nf' params(outdir: params.outdir,
  unicycler_additional_parameters: params.unicycler_additional_parameters, threads: params.threads)

// Flye assembler
include flye_assembly from './modules/flye.nf' params(outdir: params.outdir, lr_type: params.lr_type,
  flye_additional_parameters: params.flye_additional_parameters, threads: params.threads,
  genomeSize: params.genomeSize, prefix: params.prefix)

/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include nanopolish from './modules/nanopolish.nf' params(outdir: params.outdir,
  cpus: params.cpus, threads: params.threads, prefix: params.prefix)

// Medaka (for nanopore and pacbio? data)

/*
 * Define custom workflows
 */

workflow lreads_only_nf {
  take:
      reads
      fast5
      fast5_dir
  main:
      // User wants to use Canu
      if (params.try_canu) {
        canu_assembly(reads)
        if (params.fast5Path && params.lr_type == 'nanopore') {
          nanopolish(canu_assembly.out[1], reads, fast5, fast5_dir)
        }
      }

      // User wants to use Flye
      if (params.try_flye) {
        flye_assembly(reads)
        if (params.fast5Path && params.lr_type == 'nanopore') {
          nanopolish(flye_assembly.out[1], reads, fast5, fast5_dir)
        }
      }

      // User wants to use Unicycler
      if (params.try_unicycler) {
        unicycler_lreads(reads)
        if (params.fast5Path && params.lr_type == 'nanopore') {
          nanopolish(unicycler_lreads.out[1], reads, fast5, fast5_dir)
        }
      }
}

workflow shortreads_only_nf {

}

/*
 * Define main workflow
 */

workflow {

  /*
   * Long reads only assembly
   */
  if (params.assembly_type == 'longreads-only') {
    lreads_only_nf(Channel.fromPath(params.longreads), Channel.fromPath(params.fast5Path), Channel.fromPath(params.fast5Path, type: 'dir'))
  }

  /*
   * Short reads only assembly
   */
   if (params.assembly_type == 'shortreads-only') {

   }

}

/*
 * Completition message
 */
workflow.onComplete {
    println ""
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
    println "Thank you for using fmalmeida/mpgap pipeline!"
}
