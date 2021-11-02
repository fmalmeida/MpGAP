.. _manual:

Manual
======

Input files
-----------

* path to fastq files containing sequencing reads (Illumina, Nanopore or Pacbio)
* path to Pacbio subreads.bam file containing raw data (Optional)
* path to Nanopore FAST5 files containing raw data (Optional)

The input data must be provided via a samplesheet in YAML format given via the ``--input`` parameter. Please read the :ref:`samplesheet` reference page to understand how to properly create one.

.. tip::

  A samplesheet template can be downloaded with: ``nextflow run fmalmeida/mpgap --get_samplesheet``

Assembly possibilities
----------------------

The pipeline is capable of assembling Illumina, ONT and Pacbio reads in three main ways:

1. Short reads only assemblies

   + Unicycler
   + SPAdes
   + Shovill (for paired reads only).

.. note::

  `Shovill <https://github.com/tseemann/shovill>`_ is a software that can work with different assemblers as its core. The pipeline executes shovill with both ``spades``, ``skesa`` and ``megahit``, so user can compare the results.

2. Long reads only assemblies

   + Unicycler
   + Canu
   + Flye
   + Raven
   + Shasta
   + wtdbg2

3. Hybrid assemblies (using both short and long reads)

   + Unicycler
   + SPAdes
   + Haslr
   + Use short reads to correct errors (polish) in long reads assemblies.

.. note::
  
  Please read the section below to understand hybrid assembly strategies.

Hybrid assembly strategies
--------------------------

Hybrid assemblies can be produced with two available strategies that are described below. To choose the strategies adopted, users must set the ``hybrid_strategy`` parameter either from inside the YAML file as described in the :ref:`samplesheet` reference page or from the command line with ``--hybrid-strategy`` (described below in this page) which will overwrite any value set in the samplesheet.

Valid options are: ``1``, ``2`` or ``both``.

Strategy 1
"""""""""""

By using `Unicycler <https://github.com/rrwick/Unicycler#method-hybrid-assembly>`_, `Haslr <https://github.com/vpc-ccg/haslr>`_ and/or `SPAdes <https://pubmed.ncbi.nlm.nih.gov/26589280/>`_ specialized hybrid assembly modules.

.. note::

  It is achieved when not using the parameter ``--hybrid_strategy``

Strategy 2
""""""""""

By polishing (correcting errors) a long reads only assembly with Illumina reads. This will tell the pipeline to produce a long reads only assembly (with canu, wtdbg2, shasta, raven, flye or unicycler) and polish it with Pilon (for unpaired reads) or with `Unicycler-polish program <https://github.com/rrwick/Unicycler/blob/master/docs/unicycler-polish.md>`_ (for paired end reads).

Additionally, these long reads only assemblies can also be polished with Nanopolish or Racon+Medaka tools for nanopore reads and gcpp for Pacbio reads, before polishing with short reads. For that, users must properly set the parameters (``medaka_model``, ``nanopolish_fast5`` and/or ``pacbio_bam``).

Parameters documentation
------------------------

Please note that the parameters that are boolean (true or false) do not expect any value to be given for them. They must be used by itself, for example: ``--skip_spades --skip_flye``.

General parameters
""""""""""""""""""

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--output``
     - output
     - Name of directory to store assemblers results. The sample ids will be used to create sub-folder under this directory.

   * - ``--threads``
     - 3
     - Number of threads to use per process.

   * - ``--parallel_jobs``
     - NA
     - Number of processes to run in parallel. Each job can consume up to N threads (``--threads``). If not given, let's nextflow automatically handle it.

Input files
"""""""""""

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--input``
     - NA
     - Path to input samplesheet in YAML format. It is required. Please read the :ref:`samplesheet` reference page to understand how to properly create one.


Assemblies configuration
""""""""""""""""""""""""

Hybrid assembly strategies
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: 

  This overwrites any related value set inside the YAML samplesheet

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--hybrid_strategy``
     - NA
     - It tells the pipeline which hybrid assembly strategy to adopt for **all** samples. Options are: ``1``, ``2`` or ``both``. Please read the description of the hybrid assembly strategies above to better choose the right strategy.

Long reads characteristics
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: 

  These overwrite any related value set inside the YAML samplesheet

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--wtdbg2_technology``
     - NA
     - It tells the pipeline which technology the long reads are, which is required for wtdbg2. It will set a value for **all** samples. Options are: ``ont`` for Nanopore reads, ``rs`` for PacBio RSII, ``sq`` for PacBio Sequel, ``ccs`` for PacBio CCS reads. With not wanted, consider using ``--skip_wtdbg2``.
   
   * - ``--shasta_config``
     - NA
     - It tells the pipeline which shasta pre-set configuration to use. It will set a value for **all** samples. Please read the `shasta configuration manual page <https://chanzuckerberg.github.io/shasta/Configurations.html>`_ to know the available models. 

   * - ``--corrected_long_reads``
     - false
     - It tells the pipeline to interpret the long reads of **all** samples as "corrected" long reads. This will activate (if available) the options for corrected reads in the assemblers. For example: ``-corrected`` (in canu), ``--pacbio-corr|--nano-corr`` (in flye), etc. Be cautious when using this parameter. If your reads are not corrected, and you use this parameter, you will probably do not generate any contig.

Long reads polishers
""""""""""""""""""""

They are also useful for strategy 2 hybrid assemblies.

.. note:: 

  These overwrite any related value set inside the YAML samplesheet

.. list-table::
   :widths: 30 10 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--medaka_model``
     - NA
     - It tells the pipeline which available medaka model to use for **all** samples. Please read `medaka manual <https://github.com/nanoporetech/medaka#models>`_ to see available models.

   * - ``--nanopolish_max_haplotypes``
     - 1000
     - It sets the maximum number of haplotypes to be considered by Nanopolish for **all** samples. Sometimes the pipeline may crash because to much variation was found exceeding the limit.

.. note::

	 For assembly polishing with medaka models, the assembly is first polished one time with racon using the ``-m 8 -x -6 -g -8 -w 500`` as this is the dataset in which Medaka has been trained on. Therefore, the medaka polishing in this pipeline mean Racon 1X + Medaka.

Advanced assembler customization options
""""""""""""""""""""""""""""""""""""""""

.. note::

  Additional parameters must be given inside double quotes separated by blank spaces.

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--quast_additional_parameters``
     - NA
     - | Give additional parameters to Quast while assessing assembly metrics. Must be given as shown in Quast manual. E.g. ``" --large --eukaryote "``.

   * - ``--skip_canu``
     - false
     - Skip the execution of Canu

   * - ``--canu_additional_parameters``
     - NA
     - | Passes additional parameters for Canu assembler. E.g. ``" correctedErrorRate=0.075 corOutCoverage=200 "``. Must be given as shown in Canu's manual.

   * - ``--skip_flye``
     - false
     - Skip the execution of Flye

   * - ``--flye_additional_parameters``
     - NA
     - | Passes additional parameters for Flye assembler. E.g. ``" --meta --iterations 4 "``. Must be given as shown in Flye's manual.

   * - ``--skip_raven``
     - false
     - Skip the execution of Raven

   * - ``--raven_additional_parameters``
     - NA
     - | Passes additional parameters for Raven assembler. E.g. ``" --polishing-rounds 4 "``. Must be given as shown in Raven's manual.
   
   * - ``--skip_shasta``
     - false
     - Skip the execution of Shasta

   * - ``--shasta_additional_parameters``
     - NA
     - | Passes additional parameters for Raven assembler. E.g. ``" --Assembly.detangleMethod 1 "``. Must be given as shown in Shasta's manual.
   
   * - ``--skip_wtdbg2``
     - false
     - Skip the execution of Raven

   * - ``--wtdbg2_additional_parameters``
     - NA
     - | Passes additional parameters for wtdbg2 assembler. E.g. ``" -k 250 "``. Must be given as shown in wtdbg2's manual. Remember, the script called for wtdbg2 is ``wtdbg2.pl`` thus you must give the parameters used by it.

   * - ``--skip_unicycler``
     - false
     - Skip the execution of Unicycler

   * - ``--unicycler_additional_parameters``
     - NA
     - | Passes additional parameters for Unicycler assembler. E.g. ``" --mode conservative --no_correct "``. Must be given as shown in Unicycler's manual.

   * - ``--skip_spades``
     - false
     - Skip the execution of SPAdes

   * - ``--spades_additional_parameters``
     - NA
     - | Passes additional parameters for SPAdes assembler. E.g. ``" --meta --plasmids "``. Must be given as shown in Spades' manual.

   * - ``--skip_haslr``
     - false
     - Skip the execution of Haslr

   * - ``--haslr_additional_parameters``
     - NA
     - | Passes additional parameters for Haslr assembler. E.g. ``" --cov-lr 30 "``. Must be given as shown in Haslr' manual.

   * - ``--skip_shovill``
     - false
     - Skip the execution of Shovill

   * - ``--shovill_additional_parameters``
     - NA
     - | Passes additional parameters for Shovill assembler. E.g. ``" --depth 15 "``. Must be given as shown in Shovill' manual.
       | The pipeline already executes shovill with spades, skesa and megahit, so please, do not use it with shovill's ``--assembler`` parameter.

.. tip::

  All these parameters are configurable through a configuration file. We encourage users to use the configuration file since it will keep your execution cleaner and more readable. See a :ref:`config` example.
