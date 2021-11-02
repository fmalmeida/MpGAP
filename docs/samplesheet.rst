.. _samplesheet:

Samplesheet
===========

The samplesheet is a YAML document that is used to describe the input samples and its basic configuration. It is required. The input samplesheet is given using the ``--input`` parameter.

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
     - | Used to set path to illumina raw reads (paired, unpaired or both).
       | -> When using both paired and unpaired reads, the paired reads must be given first, in the order\: pair 1, pair 2, unpaired.
       | -> Otherwise, if using only paired reads, they must be given in the order\: pair 1, pair 2.
       | -> If using only unpaired reads, only one entry is expected. Check samples in the template to 1, 4 and 5 to understand it.

   * - ``pacbio``
     - Used to set path to pacbio raw reads (mutually excludable with ``nanopore``)
   
   * - ``pacbio_bam``
     - Used to set path to pacbio bam file (used in conjunction with ``pacbio`` for long reads assembly polishing with gcpp)

   * - ``nanopore``
     - Used to set path to nanopore raw reads (mutually excludable with ``pacbio``)

   * - ``nanopolish_fast5``
     - Used to set path to nanopore raw FAST5 data (used in conjunction with ``nanopore`` for long reads assembly polishing with Nanopolish)


.. note::

  The illumina tag is the only one that **must** be set in indented newlines (one line per read) as shown in the complete samplesheet example. All the other input tags **must** be set in the same line, right after the separator (":"), without quotations, white spaces or signs

YAML keys related to configuration
""""""""""""""""""""""""""""""""""

These are the tags that are used to represent/set the assembly configuration that shall be used for each sample. They must be set for each one. 

However, all these configuration can also be set via the command line. When used through the command line, it overwrites any configuration found in the YAML. Please, the :ref:`manual` reference page to understand them.

The available tags are:

.. list-table::
   :widths: 30 30 40
   :header-rows: 1

   * - Input tags
     - Default (if not set)
     - Description
   
   * - ``hybrid_strategy``
     - 1
     - This sets which strategy to run when performing hybrid assemblies. Please read the :ref:`manual` reference page to understand the adopted strategies. Options are: ``1``, ``2``, or ``both``.

   * - ``corrected_long_reads``
     - false
     - Tells whether the long reads used are corrected or not. Options: ``true``, ``false``.

   * - ``nanopolish_max_haplotypes``
     - 1000
     - It sets the max number of haplotypes to be considered by Nanopolish. Sometimes the pipeline may crash because to much variation was found exceeding the limit. Options: any integer value.
   
   * - ``medaka_model``
     - r941_min_high_g360
     - Used to polish a longreads-only assembly with Medaka. It selects a Medaka ONT sequencing model for polishing. Please read `medaka manual <https://github.com/nanoporetech/medaka#models>`_ to know the available models.
   
   * - ``shasta_config``
     - Nanopore-Oct2021
     - This selects the shasta configuration file to be used when assembling reads. It is now mandatory for shasta since its v0.8 release. Please read the `shasta configuration manual page <https://chanzuckerberg.github.io/shasta/Configurations.html>`_ to know the available models.
   
   * - ``genome_size``
     - NA
     - This sets the expected genome sizes for canu and haslr assemblers, which require this value. Options are estimatives with common suffixes, for example: ``3.7m``, ``2.8g``, etc.
   
   * - ``wtdbg2_technology``
     - | ``ont`` if input is nanopore
       | ``sq`` if input is pacbio 
     - This sets the technology of input reads. It is required by wtdbg2. Options are: ``ont`` for Nanopore reads, ``rs`` for PacBio RSII, ``sq`` for PacBio Sequel, ``ccs`` for PacBio CCS reads. 

Complete samplesheet example
""""""""""""""""""""""""""""

.. literalinclude:: ../example_samplesheet.yml
   :language: yaml