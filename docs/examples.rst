.. _examples:

***********************
Some execution examples
***********************

.. warning::

  When running hybrid assemblies or mixing short read types it is advised to **avoid not required REGEX** and write the full file path, using only the required REGEX for paired end reads when applicable. 
  Since nextflow loads inputs randomly, this is said so that the pipeline does not load all reads that match the REGEX and avoid unwanted combination of inputs.

  We are currently working in provinding a way to run multiple samples at once avoinding unwanted combination.

.. note::

  When using paired end reads it is **required** that input reads are set with the “{1,2}” pattern. For example: “SRR6307304_{1,2}.fastq”. This will properly load reads “SRR6307304_1.fastq” and “SRR6307304_2.fastq”

Illumina-only assembly with paired end reads
============================================

.. code-block:: bash

   code content

   ./nextflow run fmalmeida/mpgap \
      --outdir output \
      --threads 5 \
      --shortreads_paired "path-to/illumina_r{1,2}.fastq"

.. note::

  This command will perform an illumina-only assembly using paired end reads with Unicycler, SPAdes and Shovill assemblers.

.. note::

  Since fastq files will be found by a pattern match users MUST ALWAYS double quote as: Example "illumina/SRR\*_{1,2}.fastq.gz"

.. note::

  When using paired end reads it is **required** that input reads are set with the “{1,2}” pattern. For example: “SRR6307304_{1,2}.fastq”. This will properly load reads “SRR6307304_1.fastq” and “SRR6307304_2.fastq”

Illumina-only assembly with single end reads
============================================

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --shortreads_single "path-to/illumina_unpaired.fastq"

.. note::

  This command will perform an illumina-only assembly using unpaired reads with Unicycler and SPAdes assemblers.

.. note::
  
  Since fastq files will be found by a pattern match users MUST ALWAYS double quote as: Example "SRR9696\*.fastq.gz"

Illumina-only assembly with both paired and single end reads
============================================================

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --shortreads_paired "path-to/illumina_r{1,2}.fastq" \
     --shortreads_single "path-to/illumina_unpaired.fastq"

.. note::

  This command will perform an illumina-only assembly using both paired and unpaired reads with Unicycler and SPAdes assemblers.

.. note::

  When using paired end reads it is **required** that input reads are set with the “{1,2}” pattern. For example: “SRR6307304_{1,2}.fastq”. This will properly load reads “SRR6307304_1.fastq” and “SRR6307304_2.fastq”

.. warning::

  When running mixing short read types it is advised to **avoid not required REGEX** and write the full file path, using only the required REGEX for paired end reads when applicable. 
  Since nextflow loads inputs randomly, this is said so that the pipeline does not load all reads that match the REGEX and avoid unwanted combination of inputs.

  We are currently working in provinding a way to run multiple samples at once avoinding unwanted combination.

Long reads only with ONT reads
==============================

Take note that in this example, we also polish the resulting assembly with both Nanopolish and Medaka polishers.

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

.. note::

  This will perform a long reads only assembly using nanopore data with Canu, Raven, Flye and Unicycler assemblers. This specific command will also execute a polishing step with nanopolish (see ``--nanopolish_fast5Path``) and medaka (see ``--medaka_sequencing_model``).

.. tip::

  If neither ``--nanopolish_fast5Path`` nor ``--medaka_sequencing_model`` is set, the pipeline will not try to polish the assemblies using Nanopolish or Medaka, respectively.

Long reads only with pacbio reads
=================================

Take note that in this example, we also polish the resulting assembly with gcpp polisher. GCpp is the machine-code successor of the venerable GenomicConsensus suite which has reached EOL, with the exception of not supporting Quiver/RSII anymore.

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --skip_unicycler \
     --genomeSize 2m \
     --lr_type "pacbio" \
     --longreads "path-to/pacbio.subreads.fastq" \
     --pacbio_bams "path-to/pacbio.*.subreads.bam"

.. note::

  This will perform a long reads only assembly using pacbio data with Canu, Raven, and Flye assemblers (skipping unicycler).
  This specific command will also execute a polishing step with gcpp (see ``--pacbio_bams``).

.. tip::

  If ``--pacbio_bams`` is not set, the pipeline will not try to polish the assemblies using gcpp.

Assembly in Hybrid strategy 1
=============================

Assembling directly via Unicycler, Haslr and SPAdes modules, using Pacbio reads.

.. code-block:: bash

  ./nextflow run fmalmeida/mpgap \
     --outdir output \
     --threads 5 \
     --genomeSize 2m \
     --shortreads_paired "path-to/illumina_r{1,2}.fastq" \
     --lr_type pacbio \
     --longreads "path-to/pacbio.subreads.fastq"

.. note::

  This command will execute a hybrid assembly directly through Unicycler's, Haslr's and SPAdes' hybrid assembly modules.

.. note::

  When using paired end reads it is **required** that input reads are set with the “{1,2}” pattern. For example: “SRR6307304_{1,2}.fastq”. This will properly load reads “SRR6307304_1.fastq” and “SRR6307304_2.fastq”

.. warning::

  When running hybrid assemblies or mixing short read types it is advised to **avoid not required REGEX** and write the full file path, using only the required REGEX for paired end reads when applicable. 
  Since nextflow loads inputs randomly, this is said so that the pipeline does not load all reads that match the REGEX and avoid unwanted combination of inputs.

  We are currently working in provinding a way to run multiple samples at once avoinding unwanted combination.

Assembly in Hybrid strategy 2
=============================

By using shortreads to correct errors (polish) in longreads-only assemblies (generated with canu, raven, unicycler and/or flye). Additionally, in this example, we also execute the medaka and nanopolish poloishers before the polishing with shortreads.

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

.. note::

  This command will execute a hybrid assembly by polishing a longreads-only assembly with shortreads. The usage of ``nanopolish_fast5Path`` and ``medaka_sequencing_model``
  tells the pipeline to create additional assemblies where medaka and/or nanopolish are executed before Pilon (polishment with shortreads).

.. note::

  When using paired end reads it is **required** that input reads are set with the “{1,2}” pattern. For example: “SRR6307304_{1,2}.fastq”. This will properly load reads “SRR6307304_1.fastq” and “SRR6307304_2.fastq”

.. warning::

  When running hybrid assemblies or mixing short read types it is advised to **avoid not required REGEX** and write the full file path, using only the required REGEX for paired end reads when applicable. 
  Since nextflow loads inputs randomly, this is said so that the pipeline does not load all reads that match the REGEX and avoid unwanted combination of inputs.

  We are currently working in provinding a way to run multiple samples at once avoinding unwanted combination.

Running with a configuration file
=================================

.. code-block:: bash

      ./nextflow run fmalmeida/mpgap -c nextflow.config

Running and configure from an interactive graphical interface
=============================================================

.. code-block:: bash

      nf-core launch fmalmeida/mpgap
