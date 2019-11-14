.. _config:

Configuration File
""""""""""""""""""

To download a configuration file template users just need to run ``nextflow run fmalmeida/ngs-preprocess [--get_illumina_config] [--get_ont_config] [--get_pacbio_config]``

Using a config file your code is lot more clean and concise: ``nextflow run fmalmeida/ngs-preprocess -c [path-to-config]``

Default configuration:

.. code-block:: groovy

  /*
   * Configuration File to run NGS-PreProcess pipeline.
   */

  /*
   * General Parameters
   */
      // Output folder name
  params.outDir = 'output'
      // Number of threads to be used
  params.threads = 2
      // Set true or false to run or not longreads and shortreads pipeline.
  params.run_shortreads_pipeline = false
  params.run_longreads_pipeline  = false

  /*
   * Short Reads Parameters
   */
      // Loading input reads
      // Remember to use wildcards in order to allow the pipeline to get reads ID.
      // Loading examples:
      // For Paired end reads: params.shortreads = 'SRR6307304_{1,2}.fastq' & params.reads.size = '2'
      // For Single end: params.shortreads = 'SRR7128258*' & params.reads.size = '1'
  params.shortreads = ''
      // 1 for single end, 2 for paired ends. Let it blank if not needed.
  params.reads_size = 2
      //
      // TrimGalore Parameters
      //
      // These are the parameters used to set the number of bases to clip from
      // 5' end and 3' end of paired end reads in TrimGalore. 0 < value < read length.
      // Optional. Quality default is 20 (phred)
      // Clip from 5' end
  params.clip_r1 = 0
  params.clip_r2 = 0
      // Clip from 3' end
  params.three_prime_clip_r1 = 0
  params.three_prime_clip_r2 = 0
      // This one might me left blank to use 20 as default or set a integer
  params.quality_trim = 20
      //
      // Lighter error correction parameters
      // Set wheter to run or not lighter correction step.
  params.lighter_execute = false
      //
      // Which k-mer to use. Check Ligther's manual (https://github.com/mourisl/Lighter)
  params.lighter_kmer = 21
      // Bacterial Genome Size
  params.lighter_genomeSize =
      // Lighter alpha paramter. Rule of thumb: (7/C) where C is coverage.
      // If left blank, Lighter will automatically calculate the best value.
  params.lighter_alpha =
      //
      // PEAR - Software used to merge reads if wanted.
      //
      // Set to true or false if you want to execute its process.
  params.pear_execute = false

  /*
   * Long Reads Parameters
   */
      // General
      // Set which technology is your data from, pacbio or nanopore.
      // use lowercase.
  params.lreads_type = ''
      // Loading input files (Full path)
      // params.longReads must be set when extracted reads are already available. (Path to file)
      // params.fast5Path must be set when extracted reads are not available. (Path to dir)
      //
      // Set if fasta is already extracted (basecalled).
  params.longReads = ''
      // Set whether your .fastq data is barcoded or not.
  params.lreads_is_barcoded = false
      //
      // PACBIO only
      //
      // See manual (https://github.com/PacificBiosciences/pbh5tools/blob/master/doc/index.rst)
      //
      // Set the subreads bam path to be extracted
      // IF PACBIO DATA IS NOT ALREADY EXTRACTED. Path to file.
      // Always set the path to .subreads.bam and always keep
      // the relative .subreads.bam.pbi in the same directory as input
  params.pacbio_bamPath = ''
      // Set path to pacbio legacy .bas.h5 data
  params.pacbio_h5Path = ''

  /*
   * Configuring Nextflow Scopes.
   * Enable or not the production of Nextflow Reports
   */

  //Trace Report
  trace {
      enabled = false
      file = "${params.outDir}" + "/annotation_pipeline_trace.txt"
      fields = 'task_id,name,status,exit,realtime,cpus,%cpu,memory,%mem,rss'
  }

  //Timeline Report
  timeline {
      enabled = false
      file = "${params.outDir}" + "/annotation_pipeline_timeline.html"
  }

  //Complete Report
  report {
      enabled = false
      file = "${params.outDir}" + "/annotation_pipeline_nextflow_report.html"
  }

  // DO NOT CHANGE
  //Queue limit
  executor.$local.queueSize = 1
  //Docker usage
  docker.enabled = true
  docker.runOptions = '-u $(id -u):root'
