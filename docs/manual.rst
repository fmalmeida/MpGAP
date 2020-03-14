.. _manual:

******
Manual
******

Input
=====

    * path to fastq files containing sequencing reads (Illumina, Nanopore or Pacbio)
    * path to Pacbio .bam or .h5 files containing raw data

.. note::

   Users must **never** use hard or symbolic links. This will make nextflow fail.
   When setting the parameters, please **always** give full path to a hard file,
   not to a link. This will prevent file access fail.

.. tip::

  Hybrid assemblies can be achieved with Unicycler hybrid mode, by giving long and
  short reads or, by polishing an long reads only assembly. For this, users will have
  to set ``assembly_type = 'longreads-only'``, set path to Illumina reads and use the
  ``illumina_polish_longreads_contigs`` parameter.

Usage example
=============

::

   nextflow run fmalmeida/MpGAP [OPTIONS]

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

   * - ``--prefix``
     - Y
     - NA
     - Prefix for output files

   * - ``--yaml``
     - N
     - NA
     - Path to yaml file containing assemblers additional parameters

   * - ``--assembly_type``
     - Y
     - NA
     - Selects assembly mode: hybrid; illumina-only; longreads-only

   * - ``--try_canu``
     - N
     - False
     - Try to assemble data with Canu

   * - ``--canu_additional_parameters``
     - N
     - NA
     - If user wants Canu additional parameters

   * - ``--try_unicycler``
     - N
     - False
     - Try to assemble data with Unicycler

   * - ``--unicycler_additional_parameters``
     - N
     - NA
     - If user wants Unicycler additional parameters

   * - ``--try_flye``
     - N
     - False
     - Try to assemble data with Flye

   * - ``--flye_additional_parameters``
     - N
     - NA
     - If user wants Flye additional parameters

   * - ``--try_spades``
     - N
     - False
     - Try to assemble data with SPAdes

   * - ``--spades_additional_parameters``
     - N
     - NA
     - If user wants SPAdes additional parameters

   * - ``--shortreads_paired``
     - (if assembly mode is hybrid or illumina-only)
     - NA
     - Path to Illumina paired end reads

   * - ``--shortreads_single``
     - (if assembly mode is hybrid or illumina-only)
     - NA
     - Path to Illumina unpaired reads

   * - ``--ref_genome``
     - N (Only used by SPAdes to guide assembly)
     - NA
     - Path to reference genome

   * - ``--longreads``
     - (if assembly mode is hybrid or longreads-only)
     - NA
     - Path to long reads file

   * - ``--medaka_sequencing_model``
     - (if user wants to polish longreads-only assembly with Medaka)
     - NA
     - Selects ONT sequencing model for Medaka Polish step

   * - ``--lr_type``
     - (if assembly mode is hybrid or longreads-only)
     - nanopore
     - Tells whether input longreads are: pacbio or nanopore

   * - ``--nanopolish_fast5Path``
     - (if user wants to polish longreads-only assembly with Nanopolish)
     - NA
     - Sets path to dir containing FAST5 data for nanopolish step

   * - ``--pacbio_all_bam_path``
     - N
     - NA
     - Sets path to Pacbio .bam subreads file (files .bai mus be in the same directory)

   * - ``--genomeSize``
     - (If ``--try_canu`` or ``--try_flye``)
     - NA
     - Sets expected genome size

   * - ``--illumina_polish_longreads_contigs``
     - N
     - False
     - Tells the pipeline to create a long reads only assembly and polish it with short reads. By default, only
     the hybrid mode of Unicycler and SPAdes are executed. If used, users must remember to select which assemblers
     to use for a long reads only assembly first: ``--try_unicycler``, ``--try_canu`` or ``--try_flye``.

All this parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Examples
========

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
