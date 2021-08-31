/*
 * Define help message
 */
 def helpMessage() {
    log.info """
    Usage:
    nextflow run fmalmeida/mpgap [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]

    Comments:
    This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would
    cause the command to be huge. Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make
    parameterization easier and more readable.

    Creating a configuration file:
    nextflow run fmalmeida/mpgap [--get_hybrid_config] [--get_lreads_config] [--get_sreads_config]

    Show command line examples:
    nextflow run fmalmeida/mpgap --examples

    Execution Reports:
    nextflow run fmalmeida/mpgap [OPTIONS] [-with-report] [-with-trace] [-with-timeline]

    OPTIONS:

            # General Parameters


     --outdir <string>                                                          Main output directory. Results will be written in a sub-folder
                                                                                in this directory, using the input reads basename: For shortreads
                                                                                assemblies we will use the shortreads basename, for the others,
                                                                                either hybrid and longreads-only, we will use the longreads basename.

     --threads <int>                                                            Number of threads to use.

     --parallel_jobs <int>                                                      Number of jobs to run in parallel. Each job can consume up
                                                                                to N threads (--threads). Default: 1.

     --genomeSize <string>                                                      Canu, wtdbg2 and Haslr require an estimative of genome size in order
                                                                                to be executed. It is optional for Flye. Examples: 5.6m; 1.2g

            # Input parameters.
            # The pipeline will choose between: hybrid, shortreads or longreads only assemblies
            # based on the combination of input files given
            # Remember to always quote file paths.

     --shortreads_paired <string>                                               Path to Illumina paired end reads. E.g. "read_pair_{1,2}.fastq"

     --shortreads_single <string>                                               Path to Illumina single end reads.

     --longreads <string>                                                       Path to longreads in FASTA or FASTQ formats.

     --lr_type <string>                                                         Sets wich type of long reads are being used: pacbio or nanopore

     --wtdbg2_technology <string>                                               When assembling pacbio long reads with wtdbg2, it is necessary to 
                                                                                tell the pipeline whether reads are rs, sq or ccs, so it is properly 
                                                                                passed to the assembler. Which could take value "rs" for 
                                                                                PacBio RSII, "sq" for PacBio Sequel, "ccs" for PacBio CCS reads                                             

     --corrected_lreads                                                         Tells the pipeline to interpret the long reads as "corrected" long reads.
                                                                                This will activate (if available) the options for corrected reads in the
                                                                                assemblers: -corrected (in canu), --pacbio-corr|--nano-corr (in flye), etc.
                                                                                Be cautious when using this parameter. If your reads are not corrected, and
                                                                                you use this parameter, you will probably do not generate any contig.

            # Activation of alternative hybrid assembly mode
            # Only useful when giving both short and long reads

     --strategy_2                                                               This tells the pipeline to execute an alternative hybrid method
                                                                                instead of running Unicycler/SPAdes default hybrid workflows.
                                                                                This creates a longreads-only assembly with Canu, Unicycler, Raven
                                                                                or Flye and polish it with shortreads using Pilon. This represents
                                                                                another hybrid methodology.

                                                                                It is also possible to polish the longreads-only assemblies with Nanopolish,
                                                                                Medaka or gcpp (depending on the sequencing technology) before polishing
                                                                                it with shortreads. For that, users must check the longreads parameters:
                                                                                --medaka_sequencing_model, --nanopolish_fast5Path and --pacbio_bams

            # Assembly polishing using long reads raw data
            # Parameters useful for polishing longreads-only assemblies
            # Polishers ==> ONT: Nanopolish or Medaka; Pacbio: gcpp.

    --medaka_sequencing_model <string>                                          Tells Medaka polisher which model to use according to the basecaller
                                                                                used. For example the model named r941_min_fast_g303 should be used
                                                                                with data from MinION (or GridION) R9.4.1 flowcells using the fast
                                                                                Guppy basecaller version 3.0.3. Where a version of Guppy has been
                                                                                used without an exactly corresponding medaka model, the medaka model
                                                                                with the highest version equal to or less than the guppy version
                                                                                should be selected. [ Default: r941_min_high_g360 ]

                                                                                Models available: r103_min_high_g345, r103_min_high_g360, r103_prom_high_g360,
                                                                                r103_prom_snp_g3210, r103_prom_variant_g3210, r10_min_high_g303, r10_min_high_g340,
                                                                                r941_min_fast_g303, r941_min_high_g303, r941_min_high_g330, r941_min_high_g340_rle,
                                                                                r941_min_high_g344, r941_min_high_g351, r941_min_high_g360, r941_prom_fast_g303,
                                                                                r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344, r941_prom_high_g360,
                                                                                r941_prom_high_g4011, r941_prom_snp_g303, r941_prom_snp_g322, r941_prom_snp_g360,
                                                                                r941_prom_variant_g303, r941_prom_variant_g322, r941_prom_variant_g360


     --nanopolish_fast5Path <string>                                            Path to directory containing FAST5 files for given reads.
                                                                                Whenever set, the pipeline will execute a polishing step
                                                                                with Nanopolish. This makes the pipeline extremely SLOW!!

     --cpus <int>                                                               Number of cores to run nanopolish in parallel.
                                                                                Beware of your system limits. Default: 2.

     --nanopolish_max_haplotypes <int>                                          This parameter sets to nanopolish the max number of haplotypes to be considered.
                                                                                Sometimes the pipeline may crash because to much variation was found exceeding the
                                                                                limit. Try augmenting this value (Default: 1000)

     --pacbio_bams <string>                                                     Path to all subreads bam files for given reads. Whenever set, the pipeline
                                                                                will execute a polishing step with gcpp. GCpp is the machine-code successor 
                                                                                of the venerable GenomicConsensus suite which has reached EOL, with the 
                                                                                exception of not supporting Quiver/RSII anymore. In order to nextflow properly
                                                                                use it, one needs to store all the data, from all the cells in one single 
                                                                                directory and set the filepath as "some/data/*bam".

            # Advanced parameters
            # Controlling the execution of assemblers and QUAST
            # Also adding the possibility to pass additional parameters to them

     --skip_spades                                                              Skip assembly with Spades (hybrid and shortreads only assembler)

     --skip_shovill                                                             Skip assembly with Shovill (paired shortreads only assembler)

     --skip_unicycler                                                           Skip assembly with Unicycler (hybrid, long and short reads only assembler)

     --skip_haslr                                                               Skip assembly with Haslr (hybrid assembler)

     --skip_canu                                                                Skip assembly with Canu (longreads only assembler)

     --skip_flye                                                                Skip assembly with Flye (longreads only assembler)

     --skip_raven                                                               Skip assembly with Raven (longreads only assembler)

     --skip_wtdbg2                                                              Skip assembly with wtdbg2 (longreads only assembler)

     --skip_shasta                                                              Skip assembly with Shasta (nanopore longreads only assembler)

     --quast_additional_parameters <string>                                     Give additional parameters to Quast while assessing assembly metrics.
                                                                                Must be in quotes and separated by spaces. 
                                                                                Must be given as shown in Quast manual. 
                                                                                E.g. " --large --eukaryote ".

     --spades_additional_parameters <string>                                    Give additional parameters to Spades assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Spades manual.
                                                                                E.g. " --meta --plasmids "

     --shovill_additional_parameters <string>                                   Give additional parameters to Shovill assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Shovill manual.
                                                                                E.g. " --depth 15 --assembler skesa "

     --unicycler_additional_parameters <string>                                 Give additional parameters to Unicycler assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Unicycler manual.
                                                                                E.g. " --mode conservative --no_correct "

     --haslr_additional_parameters <string>                                     Give additional parameters to Haslr assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Haslr manual.
                                                                                E.g. " --cov-lr 30 "

     --canu_additional_parameters <string>                                      Give additional parameters to Canu assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Canu manual.
                                                                                E.g. " correctedErrorRate=0.075 corOutCoverage=200 "

     --flye_additional_parameters <string>                                      Give additional parameters to Flye assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Flye manual.
                                                                                E.g. " --meta --iterations 4 "

     --raven_additional_parameters <string>                                     Give additional parameters to Raven assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in Raven manual.
                                                                                E.g. " --polishing-rounds 4 "

     --wtdbg2_additional_parameters <string>                                    Give additional parameters to wtdbg2 assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in wtdbg2 manual.
                                                                                E.g. " --tidy-reads 5000 "

     --shasta_additional_parameters <string>                                    Give additional parameters to shasta assembler. Must be in quotes
                                                                                and separated by one space. Must be given as shown in shasta manual.
                                                                                E.g. " --Reads.minReadLength 5000 "  
    """.stripIndent()
 }
