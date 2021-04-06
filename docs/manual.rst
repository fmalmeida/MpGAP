.. _manual:

******
Manual
******

Input files
===========

    + path to fastq files containing sequencing reads (Illumina, Nanopore or Pacbio)
    + path to Pacbio subreads.bam file containing raw data (Optional)
    + path to Nanopore FAST5 files containing raw data (Optional)

.. note::

  Users must **never** use hard or symbolic links. This will make nextflow fail. When setting the parameters, please **always** give full path to a hard file, not to a link. This will prevent file access fail. Remember to **always** write input paths inside double quotes.

.. warning::

  Remember: the pipeline does not concatenate the reads. Whenever you use a pattern such as \* the pipeline will assemble each pair separately. When doing hybrid assemblies or mixing read types it is advised to **not use REGEX** and instead write the full file path.

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
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--outdir``
     - Y
     - output
     - Name of directory to store output values. Input reads basenames will be used to create sub-folder under this directory.

   * - ``--genomeSize``
     - Y (for Canu and Haslr assemblers)
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
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--shortreads_paired``
     - Y (for hybrid and illumina-only modes)
     - NA
     - Path to Illumina paired end reads. E.g. "read_pair\_{1,2}.fastq".

   * - ``--shortreads_single``
     - Y (for hybrid and illumina-only modes)
     - NA
     - Path to Illumina single end reads. E.g. "reads\*.fastq".

   * - ``--longreads``
     - Y (for hybrid and longreads-only modes)
     - NA
     - Path to longreads in FASTA or FASTQ formats.

   * - ``--lr_type``
     - Y (for hybrid and longreads-only modes)
     - nanopore
     - Tells whether input longreads are: pacbio or nanopore.

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
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--strategy_2``
     - N
     - False
     - | Tells the pipeline to create a long reads only assembly and polish it with short reads.
       | By default, the hybrid modes of Unicycler, Haslr and SPAdes are executed.
       | This parameter tells to excute the hybrid strategy 2 (longreads -> polish) instead of Unicycler/Haslr/SPAdes hybrid modes.

Long reads assembly polishing parameters (also used for hybrid strategy 2)
--------------------------------------------------------------------------

Long reads only assemblies can also be polished with Nanopolish or Racon+Medaka tools for nanopore reads and Arrow for Pacbio reads. For that, users must properly set the parameters. given below.

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

   * - ``--pacbio_all_bam_path``
     - N
     - NA
     - | Path to all subreads.bam files for the given reads (can be '\*.bam')
       | In order to nextflow properly use it, one needs to store all the data, from all the cells in one single directory and set the filepath as "some/data/\*bam".
       | Whenever set, the pipeline will execute a polishing step with VarianCaller through arrow.
       | Arrow is supported for PacBio Sequel data and RS data with the P6-C4 chemistry.

Advanced assembler customization options
----------------------------------------

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--skip_canu``
     - N
     - False
     - Skip the execution of Canu

   * - ``--canu_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Canu assembler. E.g. 'correctedErrorRate=0.075 corOutCoverage=200'.
       | Must be given as shown in Canu's manual.

   * - ``--skip_flye``
     - N
     - False
     - Skip the execution of Flye

   * - ``--flye_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Flye assembler. E.g. '--meta --iterations 4'.
       | Must be given as shown in Flye's manual.

   * - ``--skip_raven``
     - N
     - False
     - Skip the execution of Raven

   * - ``--raven_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Raven assembler. E.g. '--polishing-rounds 4'.
       | Must be given as shown in Raven's manual.

   * - ``--skip_unicycler``
     - N
     - False
     - Skip the execution of Unicycler

   * - ``--unicycler_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Unicycler assembler. E.g. '--mode conservative --no_correct'.
       | Must be given as shown in Unicycler's manual.

   * - ``--skip_spades``
     - N
     - False
     - Skip the execution of SPAdes

   * - ``--spades_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for SPAdes assembler. E.g. '--meta --plasmids'.
       | Must be given as shown in Spades' manual.

   * - ``--skip_haslr``
     - N
     - False
     - Skip the execution of Haslr

   * - ``--haslr_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Haslr assembler. E.g. '--cov-lr 30'.
       | Must be given as shown in Haslr' manual.

   * - ``--skip_shovill``
     - N
     - False
     - Skip the execution of Shovill

   * - ``--shovill_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Shovill assembler. E.g. '--depth 15 --assembler skesa'.
       | Must be given as shown in Shovill' manual.

Container manager
"""""""""""""""""

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--singularity``
     - N
     - False
     - Use Singularity instead of Docker to manage containers?

.. tip::

  All these parameters are configurable through a configuration file. We encourage users to use the configuration file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Usage examples
==============

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
