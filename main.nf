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

     --outdir <string>                              Output directory name

     --threads <int>                                Number of threads to use

     --assembly_type <string>                       Selects assembly mode: hybrid, illumina-only or longreads-only

     --try_canu                                     Execute assembly with Canu. Multiple assemblers can be chosen.

     --canu_additional_parameters <string>          Give additional parameters to Canu assembler. Must be in quotes
                                                    and separated by one space. Must be given as shown in Canu manual.
                                                    E.g. 'correctedErrorRate=0.075 corOutCoverage=200'

     --try_unicycler                                Execute assembly with Unicycler. Multiple assemblers can be chosen.

     --unicycler_additional_parameters <string>     Give additional parameters to Unicycler assembler. Must be in quotes
                                                    and separated by one space. Must be given as shown in Unicycler manual.
                                                    E.g. '--mode conservative --no_correct'

     --try_flye                                     Execute assembly with Flye. Multiple assemblers can be chosen.

     --flye_additional_parameters <string>          Give additional parameters to Flye assembler. Must be in quotes
                                                    and separated by one space. Must be given as shown in Flye manual.
                                                    E.g. '--meta --iterations 4'

     --try_spades                                   Execute assembly with Spades. Multiple assemblers can be chosen.

     --spades_additional_parameters <string>        Give additional parameters to Spades assembler. Must be in quotes
                                                    and separated by one space. Must be given as shown in Spades manual.
                                                    E.g. '--meta --plasmids'


                       Parameters for illumina-only mode. Can be executed by SPAdes and Unicycler assemblers.
                       Users can use paired or single end reads. If both types are given at once, assemblers
                       will be executed with a mix of both.

     --shortreads_paired <string>                   Path to Illumina paired end reads.

     --shortreads_single <string>                   Path to Illumina single end reads.

                       Parameters for hybrid mode. Can be executed by spades and unicycler assemblers.

     --shortreads_paired <string>                   Path to Illumina paired end reads.

     --shortreads_single <string>                   Path to Illumina single end reads.

     --longreads <string>                           Path to longreads in FASTA or FASTQ formats.

     --lr_type <string>                             Sets wich type of long reads are being used: pacbio or nanopore

                       Parameters for longreads-only mode. Can be executed by canu, flye and unicycler assemblers.
                       In the end, long reads only assemblies can be polished with illumina reads through pilon.

     --longreads <string>                           Path to longreads in FASTA or FASTQ formats.

     --medaka_sequencing_model <string>             Tells Medaka polisher which model to use according to the basecaller
                                                    used. For example the model named r941_min_fast_g303 should be used
                                                    with data from MinION (or GridION) R9.4.1 flowcells using the fast
                                                    Guppy basecaller version 3.0.3. Where a version of Guppy has been
                                                    used without an exactly corresponding medaka model, the medaka model
                                                    with the highest version equal to or less than the guppy version
                                                    should be selected. Models available: r941_min_fast_g303,
                                                    r941_min_high_g303, r941_min_high_g330, r941_min_high_g344,
                                                    r941_prom_fast_g303, r941_prom_high_g303, r941_prom_high_g344,
                                                    r941_prom_high_g330, r10_min_high_g303, r10_min_high_g340,
                                                    r103_min_high_g345, r941_prom_snp_g303, r941_prom_variant_g303,
                                                    r941_min_high_g340_rle.


     --nanopolish_fast5Path <string>                Path to directory containing FAST5 files for given reads.
                                                    Whenever set, the pipeline will execute a polishing step
                                                    with Nanopolish. This makes the pipeline extremely SLOW!!

     --pacbio_all_bam_path <string>                 Path to all subreads bam files for given reads. Whenever set, the pipeline
                                                    will execute a polishing step with VarianCaller with arrow.
                                                    Arrow is supported for PacBio Sequel data and RS data with the P6-C4 chemistry

     --genomeSize <string>                          Canu and Flye require an estimative of genome size in order
                                                    to be executed. Examples: 5.6m; 1.2g

     --lr_type <string>                             Sets wich type of long reads are being used: pacbio or nanopore

     --illumina_polish_longreads_contigs            This tells the pipeline to polish long reads only assemblies
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
  params.medaka_sequencing_model = ''
  params.nanopolish_fast5Path = ''
  params.pacbio_all_bam_path = ''
  params.lr_type = ''
  params.shortreads_paired = ''
  params.shortreads_single = ''
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
summary['Fast5 files dir']   = params.nanopolish_fast5Path
summary['Long Reads']   = params.longreads
summary['Short single end reads']   = params.shortreads_single
summary['Short paired end reads']   = params.shortreads_paired
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
 * Modules for assembling short reads
 */

// SPAdes paired
include spades_sreads_paired_assembly from './modules/spades_sreads_paired.nf' params(outdir: params.outdir,
  threads: params.threads, spades_additional_parameters: params.spades_additional_parameters,
  shortreads_paired: params.shortreads_paired)

// SPAdes single
include spades_sreads_single_assembly from './modules/spades_sreads_single.nf' params(outdir: params.outdir,
  threads: params.threads, spades_additional_parameters: params.spades_additional_parameters,
  shortreads_paired: params.shortreads_single)

/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include nanopolish from './modules/nanopolish.nf' params(outdir: params.outdir,
  cpus: params.cpus, threads: params.threads, prefix: params.prefix)

// Medaka (for nanopore data)
include medaka from './modules/medaka.nf' params(medaka_sequencing_model: params.medaka_sequencing_model,
  threads: params.threads, outdir: params.outdir, prefix: params.prefix)

// VariantCaller Pacbio
include variantCaller from './modules/variantCaller.nf' params(threads: params.threads,
  outdir: params.outdir, prefix: params.prefix)

/*
 * Define custom workflows
 */

                                /*
                                 * LONG READS ONLY WORKFLOWS
                                 */

workflow lreadsonly_nf {
  take:
      reads
  main:
      // User wants to use Canu
      if (params.try_canu) {
        canu_assembly(reads)
        }

      // User wants to use Flye
      if (params.try_flye) {
        flye_assembly(reads)
        }

      // User wants to use Unicycler
      if (params.try_unicycler) {
        unicycler_lreads(reads)
        }
}

// With Nanopolish
workflow lreadsonly_nanopolish_nf {
  take:
      reads
      fast5
      fast5_dir
  main:
      // User wants to use Canu
      if (params.try_canu) {
        canu_assembly(reads)
        nanopolish(canu_assembly.out[1], reads, fast5, fast5_dir)
        }

      // User wants to use Flye
      if (params.try_flye) {
        flye_assembly(reads)
        nanopolish(flye_assembly.out[1], reads, fast5, fast5_dir)
        }

      // User wants to use Unicycler
      if (params.try_unicycler) {
        unicycler_lreads(reads)
        nanopolish(unicycler_lreads.out[1], reads, fast5, fast5_dir)
        }
}

// With Medaka
workflow lreadsonly_medaka_nf {
  take:
      reads
  main:
      // User wants to use Canu
      if (params.try_canu) {
        canu_assembly(reads)
        medaka(canu_assembly.out[1], reads)
      }

      // User wants to use Flye
      if (params.try_flye) {
        flye_assembly(reads)
        medaka(flye_assembly.out[1], reads)
      }

      // User wants to use Unicycler
      if (params.try_unicycler) {
        unicycler_lreads(reads)
        medaka(unicycler_lreads.out[1], reads)
      }
}

// With Arrow (PacBio)
workflow lreadsonly_variantCaller_nf {
  take:
      reads
      bamFile
      nBams
  main:
      // User wants to use Canu
      if (params.try_canu) {
        canu_assembly(reads)
        variantCaller(canu_assembly.out[1], bamFile, nBams)
        }

      // User wants to use Flye
      if (params.try_flye) {
        flye_assembly(reads)
        variantCaller(flye_assembly.out[1], bamFile, nBams)
        }

      // User wants to use Unicycler
      if (params.try_unicycler) {
        unicycler_lreads(reads)
        variantCaller(unicycler_lreads.out[1], bamFile, nBams)
        }
}

                                  /*
                                   * SHORT READS ONLY WORKFLOWS
                                   */

// Paired end reads
workflow sreads_only_paired_nf {
  take:
      paired_reads
  main:
      // User wants to use SPAdes
      if (params.try_spades) {
        spades_sreads_paired_assembly(paired_reads)
      }
}

// Single end reads
workflow sreads_only_single_nf {
  take:
      single_reads
  main:
      // User wants to use SPAdes
      if (params.try_spades) {
        spades_sreads_single_assembly(single_reads)
      }
}

                                  /*
                                   * DEFINE MAIN WORKFLOW
                                   */

workflow {

                            /*
                             * Long reads only assembly without polish
                             */

  if (params.assembly_type == 'longreads-only' && params.medaka_sequencing_model == '' &&
      params.nanopolish_fast5Path == '' && params.pacbio_all_bam_path == '') {
    lreadsonly_nf(Channel.fromPath(params.longreads))
  }

                            /*
                             * Long reads only assembly with polish
                             */

  // With Medaka
  if (params.assembly_type == 'longreads-only' && params.lr_type == 'nanopore' &&
      params.medaka_sequencing_model) {
    lreadsonly_medaka_nf(Channel.fromPath(params.longreads))
  }

  // With Nanopolish
  if (params.assembly_type == 'longreads-only' && params.lr_type == 'nanopore' &&
      params.nanopolish_fast5Path) {
    lreadsonly_nanopolish_nf(Channel.fromPath(params.longreads),
                    Channel.fromPath(params.nanopolish_fast5Path),
                    Channel.fromPath(params.nanopolish_fast5Path, type: 'dir'))
  }

  // With Pacbio Genomic Consensus with bax files
  if (params.assembly_type == 'longreads-only' && params.lr_type == 'pacbio' &&
      params.pacbio_all_bam_path) {
    lreadsonly_variantCaller_nf(Channel.fromPath(params.longreads),
                    Channel.fromPath(params.pacbio_all_bam_path),
                    Channel.fromPath(params.pacbio_all_bam_path).count().subscribe { println it })
  }

                          /*
                           * Short reads only assembly
                           */

   if (params.assembly_type == 'illumina-only') {
     // Using paired end reads
     if (params.shortreads_paired) {
       sreads_only_paired_nf(Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ))
     }

     // Using single end reads
     if (params.shortreads_single) {
       sreads_only_single_nf(Channel.fromPath(params.shortreads_single))
     }
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
