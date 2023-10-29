.. _manual:

Manual
======

Input files
-----------

* path to fastq files containing sequencing reads (Illumina, Nanopore or Pacbio)
* path to Pacbio subreads.bam file containing raw data (Optional)
* path to Nanopore FAST5 files containing raw data (Optional)

The input data must be provided via a samplesheet in YAML format given via the ``--input`` parameter. Please read the :ref:`samplesheet reference page<samplesheet>` to understand how to properly create one.

.. tip::

  A samplesheet template can be downloaded with: ``nextflow run fmalmeida/mpgap --get_samplesheet``

Assembly possibilities
----------------------

The pipeline is capable of assembling Illumina, ONT and Pacbio reads in three main ways:

1. Short reads only assemblies

   + Unicycler
   + SPAdes
   + Megahit
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

Hybrid assembly strategies
--------------------------

Hybrid assemblies can be produced with two available strategies that are described below. To choose the strategies adopted, users must set the ``hybrid_strategy`` parameter either from inside the YAML file (which will overwrite, for that sample, any value set) as described in the :ref:`samplesheet reference page <samplesheet>` or with the ``--hybrid-strategy`` parameter to set a new default value for all samples.

Valid options are: ``1``, ``2`` or ``both``.

Strategy 1
""""""""""

By using `Unicycler <https://github.com/rrwick/Unicycler#method-hybrid-assembly>`_, `Haslr <https://github.com/vpc-ccg/haslr>`_ and/or `SPAdes <https://pubmed.ncbi.nlm.nih.gov/26589280/>`_ specialized hybrid assembly modules.

.. note::

  It is achieved when using ``--hybrid_strategy 1`` or ``--hybrid_strategy both``

Strategy 2
""""""""""

By polishing (correcting errors) a long reads only assembly with Illumina reads. This will tell the pipeline to produce a long reads only assembly (with canu, wtdbg2, shasta, raven, flye or unicycler) and polish it with `Pilon <https://github.com/broadinstitute/pilon>`_. By default, it runs 4 rounds of polishing (``params.pilon_polish_rounds``).

.. note::

  It is achieved when using ``--hybrid_strategy 2`` or ``--hybrid_strategy both``

Additionally, these long reads only assemblies can also be polished with Nanopolish or Racon+Medaka tools for nanopore reads and gcpp for Pacbio reads, before polishing with short reads. For that, users must properly set the parameters (``medaka_model``, ``nanopolish_fast5`` and/or ``pacbio_bam``).

Parameters documentation
------------------------

Please note that, through the command line, the parameters that are boolean (true or false) do not expect any value to be given for them. They must be used by itself, for example: ``--skip_spades --skip_flye``.

.. tip::

  All parameters described can be configured through a configuration file. We encourage users to use it since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Input/Output parameters
"""""""""""""""""""""""

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--input``
     - NA
     - Path to input samplesheet in YAML format. It is required. Please read the :ref:`samplesheet` reference page to understand how to properly create one.
   
   * - ``--output``
     - output
     - Name of directory to store assemblers results. The sample ids will be used to create sub-folder under this directory.

Max job request
""""""""""""""""

Please be aware that, these parameters set the maximum amount of resources a single job can request to avoid surpassing your system limits. These resources are requested dynamically by nextflow https://nf-co.re/usage/configuration#tuning-workflow-resources (we implemented the same principle).

Also, by its nature, nextflow tries to execute as much as he can in parallel. So even if you turn these parameters to 4 cpus, for example, if your machine has 12 threads, nextflow may try to run 3 jobs in parallel, each with 4 threads. To limit the number of parallel jobs that nextflow launches at a time, you must configure it in a custom config file as shown here: https://www.nextflow.io/docs/latest/config.html#scope-executor.

.. note::
  
  The "local" executor is nextflow's default, which will be the reality of most users. Limiting parallel jobs in local executor would look like this ``executor.$local.queueSize = 1``. Please, read NF manual about it.

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--max_cpus``
     - N
     - 4
     - Max number of threads to use in parallel
   
   * - ``--max_memory``
     - N
     - 6.GB
     - Max amount of memory to be used by pipeline
   
   * - ``--max_time``
     - N
     - 40.h
     - Max time for a job to run

Assemblies configuration
""""""""""""""""""""""""

All these parameters listed below (for genome size, assembly strategy, long reads characteristics and for long reads polishers) if used via the command line or from the NF config file, they will set values in a global manner for all the samples.

However, they can also be set in a sample-specific manner. If a sample has a value for one of these parameters in the samplesheet, it will overwrite the "global/default" value **for that specific sample** and use the one provided inside the YAML.

Please, refer to the :ref:`samplesheet reference page<samplesheet>` to better understand how properly set up the samplesheet.

Genome size
^^^^^^^^^^^

A few assemblers expect you to provide an expected genome size for your assembly.

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--genome_size``
     - NA
     - This sets the expected genome sizes for canu, wtdbg2 and haslr assemblers, which require this value. Options are estimatives with common suffices, for example: ``3.7m``, ``2.8g``, etc.

Hybrid assembly strategies
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--hybrid_strategy``
     - 1
     - It tells the pipeline which hybrid assembly strategy to adopt. Options are: ``1``, ``2`` or ``both``. Please read the description of the hybrid assembly strategies above to better choose the right strategy.

Long reads characteristics
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--wtdbg2_technology``
     - The pipeline will use ``ont`` for nanopore reads and ``sq`` for pacbio reads
     - It tells the pipeline which technology the long reads are, which is required for wtdbg2. Options are: ``ont`` for Nanopore reads, ``rs`` for PacBio RSII, ``sq`` for PacBio Sequel, ``ccs`` for PacBio CCS reads. With not wanted, consider using ``--skip_wtdbg2``.
   
   * - ``--shasta_config``
     - Nanopore-Oct2021
     - It tells the pipeline which shasta pre-set configuration to use when assembling nanopore reads. Please read the `shasta configuration manual page <https://chanzuckerberg.github.io/shasta/Configurations.html>`_ to know the available models. 

   * - ``--corrected_longreads``
     - false
     - It tells the pipeline to interpret the input long reads as "corrected". This will activate (if available) the options for corrected reads in the assemblers. For example: ``-corrected`` (in canu), ``--pacbio-corr|--nano-corr`` (in flye), etc. Be cautious when using this parameter. If your reads are not corrected, and you use this parameter, you will probably do not generate any contig.

Long reads polishers
^^^^^^^^^^^^^^^^^^^^

Useful for long reads only and strategy 2 hybrid assemblies.

.. list-table::
   :widths: 30 10 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--medaka_model``
     - r941_min_high_g360
     - It tells the pipeline which available medaka model to use to polish nanopore long reads assemblies. Please read `medaka manual <https://github.com/nanoporetech/medaka#models>`_ to see available models.

   * - ``--nanopolish_max_haplotypes``
     - 1000
     - It sets the maximum number of haplotypes to be considered by Nanopolish. Sometimes the pipeline may crash because to much variation was found exceeding the limit.

.. note::

	 For assembly polishing with medaka models, the assembly is first polished one time with racon using the ``-m 8 -x -6 -g -8 -w 500`` as this is the dataset in which Medaka has been trained on. Therefore, the medaka polishing in this pipeline mean Racon 1X + Medaka.

Advanced assembler customization options
""""""""""""""""""""""""""""""""""""""""

.. note::

  Additional parameters must be set inside double quotes separated by blank spaces.

.. list-table::
   :widths: 30 10 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--quast_additional_parameters``
     - NA
     - | Give additional parameters to Quast while assessing assembly metrics. Must be given as shown in Quast manual. E.g. ``" --large --eukaryote "``.
   
   * - ``--skip_raw_assemblies_polishing``
     - false
     - | This will make the pipeline not polish raw assemblies on hybrid strategy 2.
       | For example, if a sample is assembled with flye and polished with medaka, by default, both assemblies will be passed to pilon so you can compare them.
       | If you don't need this comparison and don't want to polish the raw assembly, use this parameter.

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
     - | Passes additional parameters for Shovill assembler. E.g. ``" --depth 15 "``. Must be given as shown in Shovill's manual.
       | The pipeline already executes shovill with spades, skesa and megahit, so please, do not use it with shovill's ``--assembler`` parameter.
   
   * - ``--skip_megahit``
     - false
     - Skip the execution of Megahit

   * - ``--megahit_additional_parameters``
     - NA
     - | Passes additional parameters for Megahit assembler. E.g. ``" --presets meta-large "``. Must be given as shown in Megahit's manual.
