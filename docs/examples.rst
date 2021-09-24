.. _examples:

***********************
Some execution examples
***********************

Some considerations
===================

.. note::

  Users must **never** use hard or symbolic links. This will probably make nextflow fail. Remember to **always** write input paths inside double quotes.

.. note::

  When using paired end reads it is **required** that input reads are set with the “{1,2}” pattern. For example: “SRR6307304_{1,2}.fastq”. This will properly load reads "SRR6307304_1.fastq" and "SRR6307304_2.fastq"

.. warning::

  When running hybrid assemblies or mixing short read types it is advised to **avoid not required REGEX** and write the full file path, using only the required REGEX for paired end reads when applicable. Since nextflow randomly loads inputs, this is said to avoid unwanted combination of inputs while loading all reads that match the REGEX.

  We are currently working in provinding a way to run multiple samples at once avoinding unwanted combination.

Illumina-only assembly with paired end reads
============================================

This command will perform an illumina-only assembly using paired end reads with Unicycler, SPAdes and Shovill assemblers. Additionally, we set a custom prefix. By default the pipeline uses the reads names to create a prefix but users can set one with ``--prefix``.

.. code-block:: bash

   code content

   ./nextflow run fmalmeida/mpgap \
      --outdir output \
      --threads 5 \
      --shortreads_paired "path-to/illumina_r{1,2}.fastq" \
      --prefix "my_custom_run"

Illumina-only assembly with single end reads
============================================

This command will perform an illumina-only assembly using unpaired reads with Unicycler and SPAdes assemblers.

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --shortreads_single "path-to/illumina_unpaired.fastq"

Illumina-only assembly with both paired and single end reads
============================================================

This command will perform an illumina-only assembly using both paired and unpaired reads with Unicycler and SPAdes assemblers.

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --shortreads_paired "path-to/illumina_r{1,2}.fastq" \
     --shortreads_single "path-to/illumina_unpaired.fastq"

Long reads only with ONT reads
==============================

Take note that in this example, we also polish the resulting assembly with both Nanopolish and Medaka polishers. This will perform a long reads only assembly using nanopore data with Canu, Raven, Flye, wtdbg2, shasta and Unicycler assemblers. This specific command will also execute a polishing step with nanopolish (see ``--nanopolish_fast5Path``) and medaka (see ``--medaka_sequencing_model``).

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --genomeSize 2m \
     --lr_type nanopore \
     --longreads "path-to/ont_reads.fastq" \
     --medaka_sequencing_model r941_min_fast_g303 \
     --nanopolish_max_haplotypes 2000 \
     --nanopolish_fast5Path "path-to/fast5_pass"

.. tip::

  If neither ``--nanopolish_fast5Path`` nor ``--medaka_sequencing_model`` is set, the pipeline will not try to polish the assemblies using Nanopolish or Medaka, respectively.

Long reads only with pacbio reads
=================================

In this example, we also polish the resulting assembly with gcpp polisher. Gcpp is the machine-code successor of the venerable GenomicConsensus suite which has reached EOL, with the exception of not supporting Quiver/RSII anymore.

This will perform a long reads only assembly using pacbio data with Canu, Raven, wtdbg2, shasta and Flye assemblers (skipping unicycler). When executing wtdbg2 with pacbio reads it is required to tell with reads are RSII, Sequel, or CCS (check ``--wtdbg2_technology`` parameter. Also, this specific command will also execute a polishing step with gcpp (see ``--pacbio_bams``).

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --skip_unicycler \
     --genomeSize 2m \
     --lr_type "pacbio" \
     --wtdbg2_technology "rs" \
     --longreads "path-to/pacbio.subreads.fastq" \
     --pacbio_bams "path-to/pacbio.*.subreads.bam"

.. tip::

  If ``--pacbio_bams`` is not set, the pipeline will not try to polish the assemblies using gcpp.

Assembly in Hybrid strategy 1
=============================

This command will execute a hybrid assembly directly through Unicycler's, Haslr's and SPAdes' hybrid assembly modules using PacBio reads (``lr_type pacbio``).

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --genomeSize 2m \
     --shortreads_paired "path-to/illumina_r{1,2}.fastq" \
     --lr_type pacbio \
     --longreads "path-to/pacbio.subreads.fastq"

Assembly in Hybrid strategy 2
=============================

This command will first create longreads-only assemblies with canu, raven, unicycler and/or flye. After that, it will correct errors (polish) using shortreads with Pilon. Additionally, in this example, we also execute the medaka and nanopolish polishers before Pilon.

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --genomeSize 2m \
     --strategy_2 \
     --shortreads_paired "path-to/illumina_r{1,2}.fastq" \
     --lr_type nanopore \
     --longreads "path-to/ont_reads.fastq" \
     --medaka_sequencing_model r941_min_fast_g303 \
     --nanopolish_fast5Path "path-to/fast5_pass"

Running with a configuration file
=================================

.. code-block:: bash

      ./nextflow run fmalmeida/mpgap -c nextflow.config

Running and configure from an interactive graphical interface
=============================================================

.. code-block:: bash

      nf-core launch fmalmeida/mpgap
