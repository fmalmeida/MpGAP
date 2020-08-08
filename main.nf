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

     --outdir <string>                                                          Output directory name. Outputs are prefixed with reads IDs (basenames).

     --threads <int>                                                            Number of threads to use.

     --assembly_type <string>                                                   Selects assembly mode: hybrid, illumina-only or longreads-only

     --try_canu                                                                 Execute assembly with Canu. Multiple assemblers can be chosen.

     --canu_additional_parameters <string>                                      Give additional parameters to Canu assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Canu manual.
                                                                                E.g. 'correctedErrorRate=0.075 corOutCoverage=200'

     --try_unicycler                                                            Execute assembly with Unicycler. Multiple assemblers can be chosen.

     --unicycler_additional_parameters <string>                                 Give additional parameters to Unicycler assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Unicycler manual.
                                                                                E.g. '--mode conservative --no_correct'

     --try_flye                                                                 Execute assembly with Flye. Multiple assemblers can be chosen.

     --flye_additional_parameters <string>                                      Give additional parameters to Flye assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Flye manual.
                                                                                E.g. '--meta --iterations 4'

     --try_spades                                                               Execute assembly with Spades. Multiple assemblers can be chosen.

     --spades_additional_parameters <string>                                    Give additional parameters to Spades assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Spades manual.
                                                                                E.g. '--meta --plasmids'



                                                          Short reads only assembly



                       It can be executed by SPAdes and Unicycler assemblers. Users can use paired or single end reads.
                       If both types are given at once, assemblers will be executed using both. Remember to always
                       write the paths with regex (*, {1,2}, etc.) inside double quotes.



     --shortreads_paired <string>                                               Path to Illumina paired end reads. E.g. "read_pair_{1,2}.fastq"

     --shortreads_single <string>                                               Path to Illumina single end reads. E.g. "reads*.fastq"



                                                          Hybrid Assembly


                       Parameters for hybrid mode. Can be executed by spades and unicycler assemblers.

     --shortreads_paired <string>                                               Path to Illumina paired end reads.

     --shortreads_single <string>                                               Path to Illumina single end reads.

     --longreads <string>                                                       Path to longreads in FASTA or FASTQ formats.

     --lr_type <string>                                                         Sets wich type of long reads are being used: pacbio or nanopore




                                                          Long reads only assembly



                       Parameters for longreads-only mode. Can be executed by canu, flye and unicycler assemblers.
                       In the end, long reads only assemblies can be polished with illumina reads through pilon.

     --longreads <string>                                                       Path to longreads in FASTA or FASTQ formats.

     --lr_type <string>                                                         Sets wich type of long reads are being used: pacbio or nanopore

     --medaka_sequencing_model <string>                                         Tells Medaka polisher which model to use according to the basecaller
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


     --nanopolish_fast5Path <string>                                            Path to directory containing FAST5 files for given reads.
                                                                                Whenever set, the pipeline will execute a polishing step
                                                                                with Nanopolish. This makes the pipeline extremely SLOW!!

     --nanopolish_max_haplotypes                                                This parameter sets to nanopolish the max number of haplotypes to be considered.
                                                                                Sometimes the pipeline may crash because to much variation was found exceeding the
                                                                                limit. Try augmenting this value (Default: 1000)

     --pacbio_all_bam_path <string>                                             Path to all subreads bam files for given reads. Whenever set, the pipeline
                                                                                will execute a polishing step with VarianCaller with arrow.
                                                                                Arrow is supported for PacBio Sequel data and RS data with the P6-C4 chemistry

     --genomeSize <string>                                                      Canu and Flye require an estimative of genome size in order
                                                                                to be executed. Examples: 5.6m; 1.2g

     --illumina_polish_longreads_contigs                                        This tells the pipeline to polish long reads only assemblies
                                                                                with Illumina reads through Pilon. This is another hybrid methodology.
                                                                                For that, users have to set path to Illumina reads through
                                                                                --shortreads_paired or --shortreads_single.



    """.stripIndent()
 }

 def exampleMessage() {
    log.info """

    Examplification on how to run fmalmeida/MpGAP pipeline using the CLI configuration

    Short reads only - PAIRED:
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type illumina-only --try_spades --try_unicycler --shortreads_paired "dataset_1/sampled/illumina_R{1,2}.fastq"

    Short reads only - SINGLE:
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type illumina-only --try_spades --try_unicycler --shortreads_single "dataset_1/sampled/illumina_single.fastq"

    Short reads only - Both PAIRED and SINGLE:
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type illumina-only --try_spades --try_unicycler --shortreads_paired "dataset_1/sampled/illumina_R{1,2}.fastq" --shortreads_single "dataset_1/sampled/illumina_single.fastq"

    Long reads only - ONT (Using both Nanopolish and Medaka):
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type longreads-only --try_canu --try_flye --try_unicycler --medaka_sequencing_model r941_min_fast_g303 \
--nanopolish_fast5Path "dataset_1/ont/fast5_pass" --nanopolish_max_haplotypes 2000 --genomeSize 2m --lr_type nanopore --longreads "dataset_1/ont/ont_reads.fastq"

    Long reads only - Pacbio (Using VariantCaller):
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type "longreads-only" --try_unicycler --try_flye --genomeSize "2m" --lr_type "pacbio" \
--longreads "E01_1/Analysis_Results/preprocessed/longreads/pacbio/m141013_011508_sherri_c100709962550000001823135904221533_s1_p0.subreads.subset.fastq" --pacbio_all_bam_path "E01_1/Analysis_Results/preprocessed/longreads/pacbio/*.subreads.bam"

    Hybrid assembly - Using both paired and single end short reads
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type hybrid --try_unicycler --shortreads_paired "dataset_1/sampled/illumina_R{1,2}.fastq" \
--shortreads_single "dataset_1/sampled/illumina_single.fastq" --lr_type nanopore --longreads "dataset_1/ont/ont_reads.fastq"
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
   new File("hybrid.config").write(new URL ("https://github.com/fmalmeida/MpGAP/raw/master/configuration_example/hybrid.config").getText())
   println ""
   println "hybrid.config file saved in working directory"
   println "After configuration, run:"
   println "nextflow run fmalmeida/MpGAP -c ./hybrid.config"
   println "Nice code!\n"

   exit 0
 }
 params.get_lreads_config = false
 if (params.get_lreads_config) {
   new File("lreads-only.config").write(new URL ("https://github.com/fmalmeida/MpGAP/raw/master/configuration_example/lreads.config").getText())
   println ""
   println "lreads.config file saved in working directory"
   println "After configuration, run:"
   println "nextflow run fmalmeida/MpGAP -c ./lreads.config"
   println "Nice code!\n"

   exit 0
 }
 params.get_sreads_config = false
 if (params.get_sreads_config) {
   new File("sreads-only.config").write(new URL ("https://github.com/fmalmeida/MpGAP/raw/master/configuration_example/sreads.config").getText())
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

// General
params.outdir = 'output'
params.threads = 4
params.cpus = 2
params.assembly_type = ''

// Assemblers?
params.try_flye = false
params.try_spades = false
params.try_canu = false
params.try_unicycler = false

// Additional parameters for assemblers
params.genomeSize = ''
params.canu_additional_parameters = ''
params.unicycler_additional_parameters = ''
params.flye_additional_parameters = ''
params.spades_additional_parameters = ''

// Short reads
params.shortreads_paired = ''
params.shortreads_single = ''

// Long reads
params.longreads = ''
params.lr_type = ''
params.medaka_sequencing_model = ''
params.nanopolish_fast5Path = ''
params.nanopolish_max_haplotypes = 1000
params.pacbio_all_bam_path = ''

// Hybrid plus
params.illumina_polish_longreads_contigs = false
params.pilon_memmory_limit = 50


/*
 * Define log message
 */
log.info "================================================================="
log.info " Docker-based, fmalmeida/mpgap, generic genome assembly Pipeline "
log.info "================================================================="
def summary = [:]
if (params.longreads) { summary['Long Reads']   = params.longreads }
if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') { summary['Fast5 files dir']   = params.nanopolish_fast5Path }
if (params.shortreads_single) { summary['Short single end reads']   = params.shortreads_single }
if (params.shortreads_paired) { summary['Short paired end reads']   = params.shortreads_paired }
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
include canu_assembly from './modules/LongReads/canu.nf' params(outdir: params.outdir, lr_type: params.lr_type,
  canu_additional_parameters: params.canu_additional_parameters, threads: params.threads,
  genomeSize: params.genomeSize)

// Unicycler assembler
include unicycler_lreads from './modules/LongReads/unicycler_lreads.nf' params(outdir: params.outdir,
  unicycler_additional_parameters: params.unicycler_additional_parameters, threads: params.threads)

// Flye assembler
include flye_assembly from './modules/LongReads/flye.nf' params(outdir: params.outdir, lr_type: params.lr_type,
  flye_additional_parameters: params.flye_additional_parameters, threads: params.threads,
  genomeSize: params.genomeSize)



/*
 * Modules for assembling short reads
 */

// SPAdes paired
include spades_sreads_paired_assembly from './modules/ShortReads/spades_sreads_paired.nf' params(outdir: params.outdir,
  threads: params.threads, spades_additional_parameters: params.spades_additional_parameters)

// SPAdes single
include spades_sreads_single_assembly from './modules/ShortReads/spades_sreads_single.nf' params(outdir: params.outdir,
  threads: params.threads, spades_additional_parameters: params.spades_additional_parameters)

// SPAdes both
include spades_sreads_both_assembly from './modules/ShortReads/spades_sreads_both.nf' params(outdir: params.outdir,
  threads: params.threads, spades_additional_parameters: params.spades_additional_parameters)

// Unicycler paired
include unicycler_sreads_paired_assembly from './modules/ShortReads/unicycler_sreads_paired.nf' params(outdir: params.outdir,
  threads: params.threads, unicycler_additional_parameters: params.unicycler_additional_parameters)

// Unicycler single
include unicycler_sreads_single_assembly from './modules/ShortReads/unicycler_sreads_single.nf' params(outdir: params.outdir,
  threads: params.threads, unicycler_additional_parameters: params.unicycler_additional_parameters)

// Unicycler both
include unicycler_sreads_both_assembly from './modules/ShortReads/unicycler_sreads_both.nf' params(outdir: params.outdir,
  threads: params.threads, unicycler_additional_parameters: params.unicycler_additional_parameters)


/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include { nanopolish as nanopolish_canu;
          nanopolish as nanopolish_unicycler;
          nanopolish as nanopolish_flye } from './modules/LongReads/nanopolish.nf' params(outdir: params.outdir, cpus: params.cpus, threads: params.threads,
                                                                                          nanopolish_max_haplotypes: params.nanopolish_max_haplotypes)

// Medaka (for nanopore data)
include { medaka as medaka_canu;
          medaka as medaka_unicycler;
          medaka as medaka_flye } from './modules/LongReads/medaka.nf' params(medaka_sequencing_model: params.medaka_sequencing_model, threads: params.threads, outdir: params.outdir)

// VariantCaller Pacbio
include { variantCaller as variantCaller_canu;
          variantCaller as variantCaller_unicycler;
          variantCaller as variantCaller_flye } from './modules/LongReads/variantCaller.nf' params(threads: params.threads, outdir: params.outdir)

/*
 * Modules for Hybrid assemblies
 */

// Unicycler hybrid
include unicycler_hybrid from './modules/Hybrid/unicycler_hybrid.nf' params(outdir: params.outdir,
  threads: params.threads, unicycler_additional_parameters: params.unicycler_additional_parameters,
  shortreads_single: params.shortreads_single, shortreads_paired: params.shortreads_paired)


// SPAdes hybrid
include spades_hybrid from './modules/Hybrid/spades_hybrid.nf' params(outdir: params.outdir,
  threads: params.threads, spades_additional_parameters: params.spades_additional_parameters,
  shortreads_single: params.shortreads_single, shortreads_paired: params.shortreads_paired,
  lr_type: params.lr_type)

// Pilon polish paired
include { pilon_polish as pilon_polish_flye; pilon_polish as pilon_polish_flye_nanopolish;
          pilon_polish as pilon_polish_flye_medaka; pilon_polish as pilon_polish_flye_variantCaller;
          pilon_polish as pilon_polish_canu; pilon_polish as pilon_polish_canu_nanopolish;
          pilon_polish as pilon_polish_canu_medaka; pilon_polish as pilon_polish_canu_variantCaller;
          pilon_polish as pilon_polish_unicycler; pilon_polish as pilon_polish_unicycler_nanopolish;
          pilon_polish as pilon_polish_unicycler_medaka; pilon_polish as pilon_polish_unicycler_variantCaller } \
\
          from './modules/Hybrid/unicycler_polish.nf' params(outdir: params.outdir, threads: params.threads)



/*
 * Module for assessing assembly qualities
 */
include { quast as quast_sreads_spades;     quast as quast_sreads_unicycler;
          quast as quast_lreads_canu;       quast as quast_lreads_flye; quast as quast_lreads_unicycler;
          quast as quast_nanopolish_canu;   quast as quast_nanopolish_flye; quast as quast_nanopolish_unicycler;
          quast as quast_medaka_canu;       quast as quast_medaka_flye; quast as quast_medaka_unicycler;
          quast as quast_variantcaller_canu; quast as quast_variantcaller_flye; quast as quast_variantcaller_unicycler;
          quast as quast_hybrid_unicycler; quast as quast_hybrid_spades } \
\
          from './modules/QualityAssessment/quast.nf' params(threads: params.threads, outdir: params.outdir, assembly_type: params.assembly_type,
            shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)




/*
 * Define custom workflows
 */




                                 /*
                                  * WORKFLOW: SHORT READS ONLY
                                  */

workflow sreads_only_nf {
  take:
      preads
      sreads
  main:
      /*
       * Paired end reads
       */

       if (params.shortreads_paired && !params.shortreads_single) {
         // SPAdes
         if (params.try_spades) {
           spades_sreads_paired_assembly(preads)
           quast_sreads_spades(spades_sreads_paired_assembly.out[1], preads)
         }
         // Unicycler
         if (params.try_unicycler) {
           unicycler_sreads_paired_assembly(preads)
           quast_sreads_unicycler(unicycler_sreads_paired_assembly.out[1], preads)
         }
       }


      /*
       * Single end reads
       */

      if (params.shortreads_single && !params.shortreads_paired) {
        // SPAdes
        if (params.try_spades) {
          spades_sreads_single_assembly(sreads)
          quast_sreads_spades(spades_sreads_single_assembly.out[1], sreads)
        }
        // Unicycler
        if (params.try_unicycler) {
          unicycler_sreads_single_assembly(sreads)
          quast_sreads_unicycler(unicycler_sreads_single_assembly.out[1], sreads)
        }
      }

      /*
       * Both Library Types
       */

      if (params.shortreads_paired && params.shortreads_single) {
        // SPAdes
        if (params.try_spades) {
          spades_sreads_both_assembly(preads, sreads)
          quast_sreads_spades(spades_sreads_both_assembly.out[1], preads.concat(sreads).collect())
        }
        // Unicycler
        if (params.try_unicycler) {
          unicycler_sreads_both_assembly(preads, sreads)
          quast_sreads_unicycler(unicycler_sreads_both_assembly.out[1], preads.concat(sreads).collect())
        }
      }

}



                                /*
                                 * LONG READS ONLY WORKFLOWS
                                 */

workflow lreadsonly_nf {
  take:
      reads
      fast5
      fast5_dir
      bamFile
      nBams
  main:
      /*
       * Canu
       */
      if (params.try_canu) {
        canu_assembly(reads)
        quast_lreads_canu(canu_assembly.out[1], reads)

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_canu(canu_assembly.out[1], reads, fast5, fast5_dir)
          quast_nanopolish_canu(nanopolish_canu.out[0], reads)
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_canu(canu_assembly.out[1], reads)
          quast_medaka_canu(medaka_canu.out[1], reads)
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_canu(canu_assembly.out[1], bamFile, nBams)
          quast_variantcaller_canu(variantCaller_canu.out[1], reads)
        }
      }

      /*
       * Flye
       */
      if (params.try_flye) {
        flye_assembly(reads)
        quast_lreads_flye(flye_assembly.out[1], reads)

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_flye(flye_assembly.out[1], reads, fast5, fast5_dir)
          quast_nanopolish_flye(nanopolish_flye.out[0], reads)
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_flye(flye_assembly.out[1], reads)
          quast_medaka_flye(medaka_flye.out[1], reads)
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_flye(flye_assembly.out[1], bamFile, nBams)
          quast_variantcaller_flye(variantCaller_flye.out[1], reads)
        }

      }

      /*
       * Unicycler
       */
      if (params.try_unicycler) {
        unicycler_lreads(reads)
        quast_lreads_unicycler(unicycler_lreads.out[1], reads)

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_unicycler(unicycler_lreads.out[1], reads, fast5, fast5_dir)
          quast_nanopolish_unicycler(nanopolish_unicycler.out[0], reads)
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_unicycler(unicycler_lreads.out[1], reads)
          quast_medaka_unicycler(medaka_unicycler.out[1], reads)
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_unicycler(unicycler_lreads.out[1], bamFile, nBams)
          quast_variantcaller_unicycler(variantCaller_unicycler.out[1], reads)
        }
      }
}

                                  /*
                                   * WORKFLOW: HYBRID
                                   */

workflow hybrid_nf {
  take:
      preads
      sreads
      lreads
  main:

      /*
       * Full hybrid mode
       */
      // SPAdes
      if (params.try_spades) {
        spades_hybrid(lreads, preads, sreads)
        quast_hybrid_spades(spades_hybrid.out[1], preads.concat(sreads).collect())
      }
      // Unicycler
      if (params.try_unicycler) {
        unicycler_hybrid(lreads, preads, sreads)
        quast_hybrid_unicycler(unicycler_hybrid.out[1], preads.concat(sreads).collect())
      }

      /*
       * Polish a long reads assembly
       */
      if (params.illumina_polish_longreads_contigs){
        /*
         * Canu
         */
        if (params.try_canu) {
          canu_assembly(lreads)
          pilon_polish_canu(canu_assembly.out[1], preads.concat(sreads).collect())
          quast_lreads_canu(pilon_polish_canu.out[1], preads.concat(sreads).collect())

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_canu(canu_assembly.out[1], lreads, fast5, fast5_dir)
            pilon_polish_canu_nanopolish(nanopolish_canu.out[0], preads.concat(sreads).collect())
            quast_nanopolish_canu(pilon_polish_canu_nanopolish.out[1], preads.concat(sreads).collect())
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_canu(canu_assembly.out[1], lreads)
            pilon_polish_canu_medaka(medaka_canu.out[1], preads.concat(sreads).collect())
            quast_medaka_canu(pilon_polish_canu_medaka.out[1], preads.concat(sreads).collect())
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_canu(canu_assembly.out[1], bamFile, nBams)
            pilon_polish_canu_variantCaller(variantCaller_canu.out[1], preads.concat(sreads).collect())
            quast_variantcaller_canu(pilon_polish_canu_variantCaller.out[1], preads.concat(sreads).collect())
          }
        }

        /*
         * Flye
         */
        if (params.try_flye) {
          flye_assembly(lreads)
          pilon_polish_flye(flye_assembly.out[1], preads.concat(sreads).collect())
          quast_lreads_flye(pilon_polish_flye.out[1], preads.concat(sreads).collect())

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_flye(flye_assembly.out[1], lreads, fast5, fast5_dir)
            pilon_polish_flye_nanopolish(nanopolish_flye.out[0], preads.concat(sreads).collect())
            quast_nanopolish_flye(pilon_polish_flye_nanopolish.out[1], preads.concat(sreads).collect())
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_flye(flye_assembly.out[1], lreads)
            pilon_polish_flye_medaka(medaka_flye.out[1], preads.concat(sreads).collect())
            quast_medaka_flye(pilon_polish_flye_medaka.out[1], preads.concat(sreads).collect())
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_flye(flye_assembly.out[1], bamFile, nBams)
            pilon_polish_flye_variantCaller(variantCaller_flye.out[1], preads.concat(sreads).collect())
            quast_variantcaller_flye(pilon_polish_flye_variantCaller.out[1], preads.concat(sreads).collect())
          }

        }

        /*
         * Unicycler
         */
        if (params.try_unicycler) {
          unicycler_lreads(lreads)
          pilon_polish_unicycler(unicycler_lreads.out[1], preads.concat(sreads).collect())
          quast_lreads_unicycler(pilon_polish_unicycler.out[1], preads.concat(sreads).collect())

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_unicycler(unicycler_lreads.out[1], lreads, fast5, fast5_dir)
            pilon_polish_unicycler_nanopolish(nanopolish_unicycler.out[0], preads.concat(sreads).collect())
            quast_nanopolish_unicycler(pilon_polish_unicycler_nanopolish.out[1], preads.concat(sreads).collect())
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_unicycler(unicycler_lreads.out[1], lreads)
            pilon_polish_unicycler_medaka(medaka_unicycler.out[1], preads.concat(sreads).collect())
            quast_medaka_unicycler(pilon_polish_unicycler_medaka.out[1], preads.concat(sreads).collect())
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_unicycler(unicycler_lreads.out[1], bamFile, nBams)
            pilon_polish_unicycler_variantCaller(variantCaller_unicycler.out[1], preads.concat(sreads).collect())
            quast_variantcaller_unicycler(pilon_polish_unicycler_variantCaller.out[1], preads.concat(sreads).collect())
          }
        }
      }

}

                                  /*
                                   * DEFINE (RUN) MAIN WORKFLOW
                                   */

workflow {

  /*
   * Long reads only assembly
   */

  if (params.assembly_type == 'longreads-only') {

    // Without any polishing (ONT and Pacbio)
    if (!params.medaka_sequencing_model && !params.nanopolish_fast5Path && !params.pacbio_all_bam_path) {
      lreadsonly_nf(Channel.fromPath(params.longreads), Channel.value(''), Channel.value(''), Channel.value(''), Channel.value(''))
    }

    // Use Nanopolish? (ONT)
    if (params.nanopolish_fast5Path && !params.pacbio_all_bam_path && params.lr_type == 'nanopore') {
      lreadsonly_nf(Channel.fromPath(params.longreads), Channel.fromPath(params.nanopolish_fast5Path),
                    Channel.fromPath(params.nanopolish_fast5Path, type: 'dir'), Channel.value(''), Channel.value(''))
    }

    // Use VariantCaller (Pacbio)
    if (!params.nanopolish_fast5Path && params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
      lreadsonly_nf(Channel.fromPath(params.longreads), Channel.value(''), Channel.value(''),
                    Channel.fromPath(params.pacbio_all_bam_path).collect(), Channel.fromPath(params.pacbio_all_bam_path).count().subscribe { println it })
    }

  }


  /*
   * Short reads only assembly
   */

   if (params.assembly_type == 'illumina-only') {

     // Using paired end reads
     if (params.shortreads_paired && !params.shortreads_single) {
       sreads_only_nf(Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ),
                      Channel.empty())
     }

     // Using single end reads
     if (params.shortreads_single && !params.shortreads_paired) {
       sreads_only_nf(Channel.empty(),
                      Channel.fromPath(params.shortreads_single, hidden: true))
     }

     // Using both paired and single end reads
     if (params.shortreads_single && params.shortreads_paired) {
       sreads_only_nf(Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ),
                      Channel.fromPath(params.shortreads_single, hidden: true))
     }
   }

  /*
   * Hybrid assembly
   */

   if (params.assembly_type == 'hybrid') {

     // Using paired end reads
     if (!params.shortreads_single && params.shortreads_paired) {
       hybrid_nf(Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ),
                 Channel.value(''),
                 Channel.fromPath(params.longreads))
     }

     // Using single end reads
     if (params.shortreads_single && !params.shortreads_paired) {
       hybrid_nf(Channel.value(''),
                 Channel.fromPath(params.shortreads_single, hidden: true),
                 Channel.fromPath(params.longreads))
     }

     // Using both paired and single end reads
     if (params.shortreads_single && params.shortreads_paired) {
       hybrid_nf(Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ),
                 Channel.fromPath(params.shortreads_single, hidden: true),
                 Channel.fromPath(params.longreads))
     }
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
