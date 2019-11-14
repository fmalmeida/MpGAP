.. _manual:

Manual
======

Overview
""""""""

An overview of all annotation steps automatically taken by the pipeline.


Input
"""""

    * path to fastq files containing sequencing reads (Illumina, Nanopore or Pacbio)
    * path to Pacbio .bam or .h5 files containing raw data

.. note::

   Users must **never** use hard or symbolic links. This will make nextflow fail.
   When setting the parameters, please **always** give full path to a hard file,
   not to a link. This will prevent file access fail.

.. tip::

  Users must choose between a pipeline for short reads or for long reads using one
  of the following parameters: ``--run_shortreads_pipeline`` or ``--run_longreads_pipeline``

  The pipeline only runs either long or short reads workflows at once.

Usage example
"""""""""""""

::

   nextflow run fmalmeida/ngs-preprocess [OPTIONS]

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--outDir``
     - Y
     - output
     - Name of directory to store output values

   * - ``--threads``
     - N
     - 2
     - Number of threads to use

   * - ``--run_shortreads_pipeline``
     - Y
     - False
     - Tells the pipeline to run the short reads pre-processing workflow

   * - ``--run_longreads_pipeline``
     - Y
     - False
     - Tells the pipeline to run the long reads pre-processing workflow

   * - ``--shortreads``
     - Y (if ``--run_shortreads_pipeline``)
     - NA
     - String Pattern to find short reads. Example: "SRR6307304_{1,2}.fastq"

   * - ``--reads_size``
     - Y (if ``--run_shortreads_pipeline``)
     - NA
     - Tells wheter input is unpaired or paired end. 1 is unpaired. 2 is paired

   * - ``--clip_r1``
     - N
     - 0
     - Number of bases to always remove from 5' of read pair 1 or from unpaired read.

   * - ``--clip_r2``
     - N
     - 0
     - Number of bases to always remove from 5' of read pair 2.

   * - ``--three_prime_clip_r1``
     - Y (if ``--run_shortreads_pipeline``)
     - 0
     - Number of bases to always remove from 3' of read pair 1 or from unpaired read

   * - ``--three_prime_clip_r2``
     - Y (if ``--run_shortreads_pipeline``)
     - 0
     - Number of bases to always remove from 3' of read pair 2.

   * - ``--quality_trim``
     - N
     - 20
     - Phred quality threshold for trimming.

   * - ``--lighter_execute``
     - N
     - False
     - Tells wheter to run or not Lighter correction tool

   * - ``--lighter_kmer``
     - Y (If ``--lighter_execute``)
     - 21
     - Lighter k-mer to use in correction step.

   * - ``--lighter_genomeSize``
     - Y (If ``--lighter_execute``)
     - NA
     - Approximate genome size

   * - ``--lighter_alpha``
     - N
     - NA
     - Lighter sample rate alpha parameter. If empty, Lighter will automatically calculate its value.

   * - ``--flash_execute``
     - N
     - False
     - If set, PEAR will be executed to merge paired end reads

   * - ``--longReads``
     - Y (If ``--run_longreads_pipeline``)
     - NA
     - Sets path to long reads fastq files (Nanopore or Pacbio). Pre-processes basecalled long reads.

   * - ``--lreads_type``
     - Y (If ``--run_longreads_pipeline``)
     - NA
     - Tells wheter your input is nanopore or pacbio data. Possibilities: pacbio | nanopore

   * - ``--lreads_is_barcoded``
     - N
     - False
     - Tells wheter your data (Nanopore or Pacbio) is barcoded or not. It will split barcodes into single files. Users with legacy pacbio data need to first produce a new barcoded_subreads.bam file.

   * - ``--pacbio_bamPath``
     - Y (If input is pacbio .bam)
     - NA
     - Path to Pacbio subreads.bam. Only used if user wants to basecall subreads.bam to FASTQ. Always keep subreads.bam and its relative subreads.bam.pbi files in the same directory

   * - ``--pacbio_h5Path``
     - Y (If input is legacy pacbio)
     - NA
     - Path to legacy .bas.h5 data. It will be used to extract reads in FASTQ file. All related .bas.h5 and .bax.h5 files MUST be in the SAME dir.


All this parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Examples
""""""""

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
