.. _config:

Configuration File
******************

To download a configuration file template users just need to run ``nextflow run fmalmeida/MpGAP [--get_hybrid_config] [--get_lreads_config] [--get_sreads_config] [--get_yaml]``

Using a config file your code is lot more clean and concise: ``nextflow run fmalmeida/MpGAP -c [path-to-config]``

Check out some `templates <https://github.com/fmalmeida/MpGAP/tree/master/configuration_example>`_.

Example of Hybrid assembly config file:
"""""""""""""""""""""""""""""""""""""""

.. code-block:: groovy

  /*
  * This is the assembly.config file. Here are the required parameters in order to execute the
  * assembly.nextflow pipeline. It is important to keep in mind that some parameters are incompatible
  * with each other such as: params.shortreads.paired and params.shortreads.single. These type of parameters
  * can have some empty values.
  *
  *                         It has the necessary parameters for a hybrid assembly.
  */

  /*
  * Long reads input files. Just needed to be specified in case of hybrid or longreads-only assembly.
  * If none of these are wanted it may be left in blank. Fast5 path (set the path to the directory containing
  * the fast5 files) are needed to perform the polish step with nanopolish (Used with Canu software). Also, it
  * is important to keep in mind that, nanopolish only works with long reads in FASTA format, not in FASTQ.
  * In order to use nanopolish, the user might give the directory path to the FAST5 files and the full path to
  * the FASTA file (the pipeline checks the format, and, if in FASTQ, it converts to FASTA).
  */
  params.longreads = '' // Already extracted in fasta or fastq
  /*
  * This parameter is used to specify the long read sequencing technology used.
  * It might be set as one of both: nanopore ; pacbio
  * If no long read is used, it is not needed.
  */
  params.lr_type = ''
  /*
  * Short reads input files. They need to be specified in case of hybrid or shortreads-only assembly.
  * If none of these are wnated it may be left in blank. The files might be single or paired ended. They just
  * need to be properly identified as the examples below.
  * Examples for illumina reads:
  * Paired: params.shortreads.paired = 'SRR6307304_{1,2}.fastq'
  * Single: params.shortreads.single = 'SRR7128258*'
  */
  params.shortreads_paired = ''
  params.shortreads_single = ''
  /*
  * Parameter for reference genome. It is not required and just used in Spades assembly pipeline.
  * It may be left in blank.
  */
  params.ref_genome = ''
  /*
  * Here we chose the assembly type wanted. This is required to perform the assembly type wanted.
  * It must be set as one of these posibilities: longreads-only ; hybrid ; illumina-only
  */
  params.assembly_type = 'hybrid'
  /*
  * Here it is set the software wanted to perform the assembly with.
  * It must be set as true to use the software and false to skip it.
  * Users must note that only unicycler has a intrinsic hybrid
  * methodology that polishes the final assembly. Also, it is and
  * SPAdes amelioration. Therefore, users can run just unicycler for that.
  */
  params.try_unicycler = true
  params.try_spades = true
  /*
  * Here are some other parameters that do not drastically influence the workflow.
  */
  //Output folder name
  params.outDir = 'output'
  //Prefix for output files
  params.prefix = 'output'
  //Number of threads to be used by each software.
  params.threads = 3
  //Number of cores to run nanopolish in parallel
  params.cpus = 2
  //Path to the ?.yaml file containing additional parameters for the software. It may be left blank.
  params.yaml = ''
