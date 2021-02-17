.. _manual:

******
Manual
******

Input
=====

    * path to fastq files containing sequencing reads (Illumina, Nanopore or Pacbio)
    * path to Pacbio subreads.bam file containing raw data
    * path to Nanopore FAST5 files containing raw data

.. note::

  Users must **never** use hard or symbolic links. This will make nextflow fail.
  When setting the parameters, please **always** give full path to a hard file,
  not to a link. This will prevent file access fail. Remember to **always** write input paths inside double quotes.

.. warning::

  Remember: the pipeline does not concatenate the reads. Whenever you use a pattern such as \* the pipeline will assemble each pair
  separately. When doing hybrid assemblies or mixing read types it is advised to **not use REGEX** and instead write the full file
  path.

.. tip::

  The parameters ``--use_tower`` and ``--tower_token`` allows the user to launch the pipeline via `nextflow tower <https://tower.nf/>`_ in order to visualize its execution.

Hybrid assembly strategies
==========================

Hybrid assemblies can be produced using one of two available strategies:

Strategy 1
----------

By using Unicycler and/or SPAdes hybrid assembly modes. For instance, it will use Unicycler hybrid mode which will first assemble a high quality assembly graph with Illumina
data and then it will use long reads to bridge the gaps. More information about Unicycler Hybrid mode can be found `here <https://github.com/rrwick/Unicycler#method-hybrid-assembly>`_.

.. note::

  It is achieved when not using the parameter ``--illumina_polish_longreads_contigs``

Strategy 2
----------

By polishing a long reads only assembly with Illumina reads. For that, users will have to use the parameter ``--illumina_polish_longreads_contigs``. This will tell the pipeline to
produce a long reads only assembly (with canu, flye or unicycler) and polish it with Pilon (for unpaired reads) or with `Unicycler-polish program <https://github.com/rrwick/Unicycler/blob/master/docs/unicycler-polish.md>`_ (for paired end reads).

.. note::

  Note that, ``--illumina_polish_longreads_contigs`` parameter is an alternative workflow, when used, it will execute ONLY strategy 2 and not both strategies.
  When false, only strategy 1 will be executed. Remember to select the desired assemblers to run with ``--try_canu``, ``--try_flye`` and/or ``--try_unicycler``

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

   * - ``--use_tower``
     - N
     - False
     - Triggers the pipeline to be launched via nextflow tower

   * - ``--tower_token``
     - Y (if ``--use_tower``)
     - NA
     - Your nextflow tower token. Used to launch the pipeline in your nextflow tower account

   * - ``--outdir``
     - Y
     - output
     - Name of directory to store output values

   * - ``--threads``
     - N
     - 3
     - Number of threads to use

   * - ``--parallel_jobs``
     - N
     - 1
     - Number of jobs to run in parallel. Each job can consume up to N threads (``--threads``)

   * - ``--assembly_type``
     - Y
     - NA
     - Selects assembly mode: hybrid; illumina-only; longreads-only

Assembler options
-----------------

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description


   * - ``--try_canu``
     - N
     - False
     - Try to assemble data with Canu

   * - ``--canu_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Canu assembler. E.g. 'correctedErrorRate=0.075 corOutCoverage=200'.
       | Must be given as shown in Canu's manual.

   * - ``--try_unicycler``
     - N
     - False
     - Try to assemble data with Unicycler

   * - ``--unicycler_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Unicycler assembler. E.g. '--mode conservative --no_correct'.
       | Must be given as shown in Unicycler's manual.

   * - ``--try_flye``
     - N
     - False
     - Try to assemble data with Flye

   * - ``--flye_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for Flye assembler. E.g. '--meta --iterations 4'.
       | Must be given as shown in Flye's manual.

   * - ``--try_spades``
     - N
     - False
     - Try to assemble data with SPAdes

   * - ``--spades_additional_parameters``
     - N
     - NA
     - | Passes additional parameters for SPAdes assembler. E.g. '--meta --plasmids'.
       | Must be given as shown in Spades' manual.

Short reads parameters (also used for hybrid)
---------------------------------------------

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

Long reads parameters (also used for hybrid)
---------------------------------------------

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--longreads``
     - Y (for hybrid and longreads-only modes)
     - NA
     - Path to longreads in FASTA or FASTQ formats.

   * - ``--lr_type``
     - Y (for hybrid and longreads-only modes)
     - nanopore
     - Tells whether input longreads are: pacbio or nanopore.

   * - ``--medaka_sequencing_model``
     - N
     - r941_min_fast_g303
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

   * - ``--genomeSize``
     - Y (for Canu and Flye assemblers)
     - NA
     - Sets expected genome size. E.g. 5.6m; 1.2g.

   * - ``--illumina_polish_longreads_contigs``
     - N
     - False
     - | Tells the pipeline to create a long reads only assembly and polish it with short reads.
       | By default, the hybrid modes of Unicycler and SPAdes are executed.
       | This parameter tells to excute the alternative hybrid method (longreads -> polish) instead of Unicycler/SPAdes hybrid modes.
       | If used, users must remember to select which assemblers to use for a long reads only assembly first: ``--try_unicycler``, ``--try_canu`` or ``--try_flye``.

All these parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Usage examples
==============

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
