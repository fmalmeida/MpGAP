.. _samplesheet:

Samplesheet
===========

The samplesheet is a required YAML document that is used to describe the input samples and, if desired, its "sample-specific" configuration. The input samplesheet is given using the ``--input`` parameter.

.. tip::

  A samplesheet template can be downloaded with: ``nextflow run fmalmeida/mpgap --get_samplesheet``
    
A guide on how to proper configure it is shown below:

Samplesheet header
""""""""""""""""""

The first line of the file must be the header followed by an indentation (two white spaces):

.. code-block:: yaml

  samplesheet:
    - ...:

Sample identification
"""""""""""""""""""""

Each sample must be identified by the tag ``id`` in the YAML file, followed by the sample's input tags (YAML keys) that shall be used by the pipeline:

.. warning::

  This value will be used to create sub-directories in the output directory. Thus, to not use white spaces.

.. code-block:: yaml

  samplesheet:
    - id: sample_1
      ...:
      ...:
    - id: sample_2
      ...:
      ...:

YAML keys related to input files
""""""""""""""""""""""""""""""""

These are the tags that are used to represent/set the input files that shall be used for each sample. The available tags are:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Input tags
     - Description

   * - ``illumina``
     - Used to set path to illumina raw reads (paired, unpaired or both).

   * - ``pacbio``
     - Used to set path to pacbio raw reads (mutually excludable with ``nanopore``)
   
   * - ``pacbio_bam``
     - Used to set path to pacbio bam file (used in conjunction with ``pacbio`` for long reads assembly polishing with gcpp)

   * - ``nanopore``
     - Used to set path to nanopore raw reads (mutually excludable with ``pacbio``)

   * - ``nanopolish_fast5``
     - Used to set path to nanopore raw FAST5 data (used in conjunction with ``nanopore`` for long reads assembly polishing with Nanopolish)

.. note::

  Note for the illumina tag/key.

  * When using both paired and unpaired reads, the paired reads must be given first, in the order\: pair 1, pair 2, unpaired.
  * Otherwise, if using only paired reads, they must be given in the order\: pair 1, pair 2.
  * If using only unpaired reads, only one entry is expected. Check samples in the template to 1, 4 and 5 to understand it.
  * The illumina tag is the only one that **must** be set in indented newlines
      * two white spaces relative to the
      * one line per read as shown in the complete samplesheet example.

.. warning::

  All the other input tags **must** be set in the same line, right after the separator (":"), without quotations, white spaces or signs.

YAML keys related to configuration
""""""""""""""""""""""""""""""""""

These are the tags that are used to represent/set the "sample-specific" assembly configuration that shall be used for each sample.

By default, if not set inside the samplesheet, the pipeline will use the configurations set via the "nextflow config file" or via the command line. Otherwise, if set inside the samplesheet, they will overwrite the pipeline's configuration for that specific sample. 

Please, the :ref:`manual reference page<manual>` the global/defaults configurations.

The available tags are:

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Input tags
     - Description
   
   * - ``hybrid_strategy``
     - This sets which strategy to run when performing hybrid assemblies. Please read the :ref:`manual` reference page to understand the adopted strategies. Options are: ``1``, ``2``, or ``both``.

   * - ``corrected_longreads``
     - Tells whether the long reads used are corrected or not. Options: ``true``, ``false``.

   * - ``nanopolish_max_haplotypes``
     - It sets the max number of haplotypes to be considered by Nanopolish. Sometimes the pipeline may crash because to much variation was found exceeding the limit. Options: any integer value.
   
   * - ``medaka_model``
     - Used to polish a longreads-only assembly with Medaka. It selects a Medaka ONT sequencing model for polishing. Please read `medaka manual <https://github.com/nanoporetech/medaka#models>`_ to know the available models.
   
   * - ``shasta_config``
     - This selects the shasta configuration file to be used when assembling reads. It is now mandatory for shasta since its v0.8 release. Please read the `shasta configuration manual page <https://chanzuckerberg.github.io/shasta/Configurations.html>`_ to know the available models.
   
   * - ``genome_size``
     - This sets the expected genome sizes for canu, wtdbg2 and haslr assemblers, which require this value. Options are estimatives with common suffices, for example: ``3.7m``, ``2.8g``, etc.
   
   * - ``wtdbg2_technology``
     - This sets the technology of input reads. It is required by wtdbg2. Options are: ``ont`` for Nanopore reads, ``rs`` for PacBio RSII, ``sq`` for PacBio Sequel, ``ccs`` for PacBio CCS reads. 

Complete samplesheet example
""""""""""""""""""""""""""""

.. literalinclude:: ../example_samplesheet.yml
   :language: yaml