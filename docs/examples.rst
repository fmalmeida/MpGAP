.. _examples:

******************
CLI usage Examples
******************

Illumina-only assembly with paired end reads
============================================

::

   ./nextflow run fmalmeida/MpGAP --threads 3 --outDir outputs/illumina_paired --prefix test --yaml ./additional.yaml --assembly_type
   "illumina-only" --try_unicycler --try_spades --shortreads_paired "../illumina/SRR9847694_{1,2}.fastq.gz"

.. note::

  This command will perform an illumina-only assembly using paired end reads with Unicycler and SPAdes assemblers.
  Since fastq files will be found by a pattern match users MUST ALWAYS double quote as: Example "illumina/SRR9847694_{1,2}.fastq.gz"

Illumina-only assembly with single end reads
============================================

::

  ./nextflow run fmalmeida/MpGAP --threads 3 --outDir outputs/illumina_single --prefix test --yaml ./additional.yaml --assembly_type
  "illumina-only" --try_unicycler --try_spades --shortreads_single "../illumina/SRR9696*.fastq.gz"

.. note::

  This command will perform an illumina-only assembly using unpaired reads with Unicycler and SPAdes assemblers.
  Since fastq files will be found by a pattern match users MUST ALWAYS double quote as: Example "SRR9696*.fastq.gz"

Long reads only with ONT reads
==============================

::

  ./nextflow run fmalmeida/MpGAP --threads 3 --outDir sample_dataset/outputs/ont --assembly_type "longreads-only" --lr_type nanopore
  --longreads sample_dataset/ont/kpneumoniae_25X.fastq --genomeSize "5.6m" --try_canu --try_flye --try_unicycler

.. note::

  This will perform a long reads only assembly using nanopore data with Canu, Flye and Unicycler assemblers. This will note execute nanopolish
  polishing step nor a polishing with Illumina data.

Long reads only with ONT reads. With polishing (USING FAST5 data).
==================================================================

::

  ./nextflow run fmalmeida/MpGAP --threads 3 --outDir sample_dataset/outputs/ont --assembly_type "longreads-only" --lr_type nanopore
  --longreads sample_dataset/ont/kpneumoniae_25X.fastq  --fast5Path ./fast5_pass --genomeSize "5.6m" --try_canu --try_flye --try_unicycler

.. note::

  This will perform a long reads only assembly using nanopore data with Canu, Flye and Unicycler assemblers. This will note execute a
  polishing step with Illumina data.

Long reads only with ONT reads. With polishing using both FAST5 and Illumina data (paired end).
===============================================================================================

::

  ./nextflow run fmalmeida/MpGAP --threads 3 --outDir sample_dataset/outputs/ont --assembly_type "longreads-only" --lr_type nanopore
  --longreads sample_dataset/ont/kpneumoniae_25X.fastq  --fast5Path ./fast5_pass --try_canu --try_flye --try_unicycler --genomeSize "5.6m"
  --shortreads_paired "../illumina/SRR9847694_{1,2}.fastq.gz" --illumina_polish_longreads_contigs

.. note::

  This will perform a long reads only assembly using nanopore data with Canu, Flye and Unicycler assemblers. In the end, it will polish the
  assembly first with FAST5 data then with Illumina data.

Assembly in Hybrid mode with Unicycler. Using Pacbio reads.
===========================================================

::

  ./nextflow run fmalmeida/MpGAP --threads 3 --outDir sample_dataset/output --assembly_type "hybrid" --lr_type pacbio
  --longreads sample_dataset/ont/ecoli_25X.fastq --shortreads_paired "../illumina/SRR9847694_{1,2}.fastq.gz" --try_unicycler

.. note::

  This command will execute a hybrid assembly directly through Unicycler's hybrid assembly mode.

Running with a configuration file
=================================

::

      ./nextflow run fmalmeida/MpGAP -c nextflow.config
