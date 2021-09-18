.. _manual:

******
Manual
******

Input files
===========

* path to fastq files containing sequencing reads (Illumina, Nanopore or Pacbio)
* path to Pacbio subreads.bam file containing raw data (Optional)
* path to Nanopore FAST5 files containing raw data (Optional)

.. note::

  Users must **never** use hard or symbolic links. This will probably make nextflow fail. Remember to **always** write input paths inside double quotes.

.. note::

  When using paired end reads it is **required** that input reads are set with the "{1,2}" pattern. For example: "SRR6307304_{1,2}.fastq". This will properly load reads "SRR6307304_1.fastq" and "SRR6307304_2.fastq".

.. warning::

  When running hybrid assemblies or mixing short read types it is advised to **avoid not required REGEX** and write the full file path, using only the required REGEX for paired end reads when applicable. Since nextflow randomly loads inputs, this is said to avoid unwanted combination of inputs while loading all reads that match the REGEX.

  We are currently working in provinding a way to run multiple samples at once avoinding unwanted combination.

Assembly possibilities
======================

The pipeline is capable of assembling Illumina, ONT and Pacbio reads in three main ways:

1. Short reads only assemblies

   + Unicycler
   + SPAdes
   + Shovill (for paired reads only)

2. Long reads only assemblies

   + Unicycler
   + Canu
   + Flye
   + Raven
   + Shasta
   + wtdbg2

3. Hybrid (both short and long reads)

   + Unicycler
   + SPAdes
   + Haslr
   + Use short reads to correct errors (polish) in long reads assemblies. See `hybrid assembly strategy 2 <https://mpgap.readthedocs.io/en/latest/manual.html#strategy-2>`_.

Parameters documentation
========================

General parameters
------------------

.. list-table::
   :widths: 15 15 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--outdir``
     - Y
     - output
     - Name of directory to store output values. Input reads basenames will be used to create sub-folder under this directory.
   
   * - ``--prefix``
     - N
     - Input reads names
     - Gives a custom prefix for sample results. By default the pipeline creates one using the input reads names

   * - ``--genomeSize``
     - | Y
       | (for Canu, wtdbg2 and Haslr assemblers)
     - NA
     - Sets expected genome size. E.g. 5.6m; 1.2g.

   * - ``--threads``
     - N
     - 3
     - Number of threads to use

   * - ``--parallel_jobs``
     - N
     - 1
     - Number of jobs to run in parallel. Each job can consume up to N threads (``--threads``)

Input files
-----------

.. list-table::
   :widths: 20 25 10 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--shortreads_paired``
     - | Y
       | (for hybrid and illumina-only modes)
     - NA
     - Path to Illumina paired end reads. E.g. "read_pair\_{1,2}.fastq".

   * - ``--shortreads_single``
     - | Y
       | (for hybrid and illumina-only modes)
     - NA
     - Path to Illumina single end reads. E.g. "reads\*.fastq".

   * - ``--longreads``
     - | Y
       | (for hybrid and longreads-only modes)
     - NA
     - Path to longreads in FASTA or FASTQ formats.

   * - ``--lr_type``
     - | Y
       | (for hybrid and longreads-only modes)
     - nanopore
     - Tells whether input longreads are: pacbio or nanopore.
   
   * - ``--wtdbg2_technology``
     - | Y
       | (when running wtdbg2 longreads-only assembly with pacbio)
     - ont
     - | When assembling pacbio long reads with wtdbg2, it is necessary to tell the pipeline
       | whether reads are "rs" for PacBio RSII, "sq" for PacBio Sequel, "ccs" for PacBio CCS reads.
       | With do not want it, consider using ``--skip_wtdbg2``.

   * - ``--corrected_lreads``
     - N
     - False
     - | Tells the pipeline to interpret the long reads as "corrected" long reads.
       |
       | This will activate (if available) the options for corrected reads in the
       | assemblers: ``-corrected`` (in canu), ``--pacbio-corr|--nano-corr`` (in flye), etc. Be cautious when using this parameter. If your reads are not corrected, and you use this parameter, you will probably do not generate any contig.

Hybrid assembly strategy
------------------------

Hybrid assemblies can be produced using one of two available strategies:

Strategy 1
^^^^^^^^^^

By using `Unicycler <https://github.com/rrwick/Unicycler#method-hybrid-assembly>`_, `Haslr <https://github.com/vpc-ccg/haslr>`_ and/or `SPAdes <https://pubmed.ncbi.nlm.nih.gov/26589280/>`_ specialized hybrid assembly modules.

.. note::

  It is achieved when not using the parameter ``--strategy_2``

Strategy 2
^^^^^^^^^^

By polishing (correcting errors) a long reads only assembly with Illumina reads. For that, users will have to use the parameter ``--strategy_2``. This will tell the pipeline to produce a long reads only assembly (with canu, raven, flye or unicycler) and polish it with Pilon (for unpaired reads) or with `Unicycler-polish program <https://github.com/rrwick/Unicycler/blob/master/docs/unicycler-polish.md>`_ (for paired end reads).

.. note::

  Note that, ``--strategy_2`` parameter is an alternative workflow, when used, it will execute ONLY strategy 2 and not both strategies. When false, only strategy 1 will be executed.

.. list-table::
   :widths: 20 30 10 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--strategy_2``
     - N
     - False
     - | Tells the pipeline to create a long reads only assembly and polish it with short reads.
       |
       | By default, the hybrid modes of Unicycler, Haslr and SPAdes are executed. This parameter tells to excute the hybrid strategy 2 (longreads -> polish) instead of Unicycler/Haslr/SPAdes hybrid modes.

Long reads assembly polishing parameters (also used for hybrid strategy 2)
--------------------------------------------------------------------------

Long reads only assemblies can also be polished with Nanopolish or Racon+Medaka tools for nanopore reads and gcpp for Pacbio reads. For that, users must properly set the parameters. given below.

.. note::

	 For assembly polishing with medaka models, the assembly is first polished one time with racon using the ``-m 8 -x -6 -g -8 -w 500`` as this is the dataset in which Medaka has been trained on. Therefore, the medaka polishing in this pipeline mean Racon 1X + Medaka.

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--medaka_sequencing_model``
     - N
     - r941_min_high_g360
     - | Used to polish a longreads-only assembly with Medaka. It selects a Medaka ONT sequencing model for polishing.
       | Please read `medaka manual <https://github.com/nanoporetech/medaka#models>`_ for more instructions.

   * - ``--nanopolish_fast5Path``
     - N
     - NA
     - | Used to polish a longreads-only assembly with Nanopolish.
       | It sets path to the directory containing all the FAST5 files containing the raw data.

   * - ``--nanopolish_max_haplotypes``
     - N
     - 1000
     - It sets the max number of haplotypes to be considered by Nanopolish. Sometimes the pipeline may crash because to much variation was found exceeding the limit.

   * - ``--pacbio_bams``
     - N
     - NA
     - | Path to all subreads.bam files for the given reads (can be '\*.bam')
       | In order to nextflow properly use it, one needs to store all the data, from all the cells in one single directory and set the filepath as "some/data/\*bam".
       |
       | Whenever set, the pipeline will execute a polishing step with gcpp. GCpp is the machine-code successor of the venerable GenomicConsensus suite which has reached EOL, with the exception of not supporting Quiver/RSII anymore.

Advanced assembler customization options
----------------------------------------

.. note::

  Additional parameters must be given inside double quotes separated by blank spaces.

.. list-table::
   :widths: 30 10 10 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--quast_additional_parameters``
     - N
     - NA
     - | Give additional parameters to Quast while assessing assembly metrics. Must be given as shown in Quast manual. E.g. ``' --large --eukaryote '``.

   * - ``--skip_canu``
     - N
     - False
     - Skip the execution of Canu

   * - ``--canu_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Canu assembler. E.g. ``' correctedErrorRate=0.075 corOutCoverage=200 '``. Must be given as shown in Canu's manual.

   * - ``--skip_flye``
     - N
     - False
     - Skip the execution of Flye

   * - ``--flye_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Flye assembler. E.g. ``' --meta --iterations 4 '``. Must be given as shown in Flye's manual.

   * - ``--skip_raven``
     - N
     - False
     - Skip the execution of Raven

   * - ``--raven_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Raven assembler. E.g. ``' --polishing-rounds 4 '``. Must be given as shown in Raven's manual.
   
   * - ``--skip_shasta``
     - N
     - False
     - Skip the execution of Shasta

   * - ``--shasta_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Raven assembler. E.g. ``' --Assembly.detangleMethod 1 '``. Must be given as shown in Shasta's manual.
   
   * - ``--skip_wtdbg2``
     - N
     - False
     - Skip the execution of Raven

   * - ``--wtdbg2_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for wtdbg2 assembler. E.g. ``' -k 250 '``. Must be given as shown in wtdbg2's manual. Remember, the script called for wtdbg2 is ``wtdbg2.pl`` thus you must give the parameters used by it.

   * - ``--skip_unicycler``
     - N
     - False
     - Skip the execution of Unicycler

   * - ``--unicycler_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Unicycler assembler. E.g. ``' --mode conservative --no_correct '``. Must be given as shown in Unicycler's manual.

   * - ``--skip_spades``
     - N
     - False
     - Skip the execution of SPAdes

   * - ``--spades_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for SPAdes assembler. E.g. ``' --meta --plasmids '``. Must be given as shown in Spades' manual.

   * - ``--skip_haslr``
     - N
     - False
     - Skip the execution of Haslr

   * - ``--haslr_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Haslr assembler. E.g. ``' --cov-lr 30 '``. Must be given as shown in Haslr' manual.

   * - ``--skip_shovill``
     - N
     - False
     - Skip the execution of Shovill

   * - ``--shovill_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Shovill assembler. E.g. ``' --depth 15 --assembler skesa '``. Must be given as shown in Shovill' manual.

.. tip::

  All these parameters are configurable through a configuration file. We encourage users to use the configuration file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Usage examples
==============

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
