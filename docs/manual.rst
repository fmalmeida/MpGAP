.. _manual:

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

  Hybrid assemblies can be achieved with Unicycler|SPAdes hybrid modes, by giving long and
  short reads or, by polishing a long reads only assembly. For the latter, users will have
  to set ``assembly_type = 'hybrid'``, set path to Illumina reads and use the
  ``illumina_polish_longreads_contigs`` parameter.

.. tip::

  When using the ``illumina_polish_longreads_contigs`` parameter it is also possible to combine
  polishings with Medaka, Nanopolish or Arrow by using setting the correct parameters:
  pacbio_all_bam_path, nanopolish_fast5Path or medaka_sequencing_model. These will tell the
  pipeline to polish the assemblies with these software before polishing with shortreads (using Pilon).

Parameters documentation
========================

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
     - Name of directory to store output values

   * - ``--threads``
     - N
     - 3
     - Number of threads to use

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
     - Passes additional parameters for Canu assembler. E.g. 'correctedErrorRate=0.075 corOutCoverage=200'. Must be given as shown in Canu's manual.

   * - ``--try_unicycler``
     - N
     - False
     - Try to assemble data with Unicycler

   * - ``--unicycler_additional_parameters``
     - N
     - NA
     - Passes additional parameters for Unicycler assembler. E.g. '--mode conservative --no_correct'. Must be given as shown in Unicycler's manual.

   * - ``--try_flye``
     - N
     - False
     - Try to assemble data with Flye

   * - ``--flye_additional_parameters``
     - N
     - NA
     - Passes additional parameters for Flye assembler. E.g. '--meta --iterations 4'. Must be given as shown in Flye's manual.

   * - ``--try_spades``
     - N
     - False
     - Try to assemble data with SPAdes

   * - ``--spades_additional_parameters``
     - N
     - NA
     - Passes additional parameters for SPAdes assembler. E.g. '--meta --plasmids'. Must be given as shown in Spades' manual.

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

   * - ``--medaka_sequencing_model``
     - N
     - r941_min_fast_g303
     - Used to polish a longreads-only assembly with Medaka. It selects a Medaka ONT sequencing model for polishing. Please read `medaka manual <https://github.com/nanoporetech/medaka#models>`_ for more instructions.

   * - ``--nanopolish_fast5Path``
     - N
     - NA
     - Used to polish a longreads-only assembly with Nanopolish. It sets path to the directory containing all the FAST5 files containing the raw data.

   * - ``--nanopolish_max_haplotypes``
     - N
     - 1000
     - It sets the max number of haplotypes to be considered by Nanopolish. Sometimes the pipeline may crash because to much variation was found exceeding the limit.

   * - ``--pacbio_all_bam_path``
     - N
     - NA
     - Path to all subreads.bam files for the given reads. Whenever set, the pipeline will execute a polishing step with VarianCaller through arrow. Arrow is supported for PacBio Sequel data and RS data with the P6-C4 chemistry.

   * - ``--genomeSize``
     - Y (for Canu and Flye assemblers)
     - NA
     - Sets expected genome size. E.g. 5.6m; 1.2g.

   * - ``--illumina_polish_longreads_contigs``
     - N
     - False
     - Tells the pipeline to create a long reads only assembly and polish it with short reads. By default, the hybrid modes of Unicycler and SPAdes are executed. This parameter tells to excute the alternative hybrid method (longreads -> polish) instead of Unicycler/SPAdes hybrid modes. If used, users must remember to select which assemblers to use for a long reads only assembly first: ``--try_unicycler``, ``--try_canu`` or ``--try_flye``.

All these parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Usage examples
==============

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
