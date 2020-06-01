.. _examples:

************************
A few CLI usage Examples
************************

.. warning::

  Remember: the pipeline does not concatenate the reads. Whenever you use a pattern such as \* the pipeline will assemble each pair
  separately. When doing hybrid assemblies or mixing read types it is advised to **not use REGEX** and instead write the full file
  path.

Illumina-only assembly with paired end reads
============================================

::

   ./nextflow run fmalmeida/MpGAP --threads 3 --outDir "outputs/illumina_paired" --prefix test --yaml "./additional.yaml"
    --assembly_type "illumina-only" --try_unicycler --try_spades --shortreads_paired "../illumina/SRR*_{1,2}.fastq.gz"

.. note::

  This command will perform an illumina-only assembly using paired end reads with Unicycler and SPAdes assemblers.
  Since fastq files will be found by a pattern match users MUST ALWAYS double quote as: Example "illumina/SRR\*_{1,2}.fastq.gz"

Illumina-only assembly with single end reads
============================================

::

  ./nextflow run fmalmeida/MpGAP --threads 3 --outDir "outputs/illumina_single" --prefix test --yaml "./additional.yaml"
  --assembly_type "illumina-only" --try_unicycler --try_spades --shortreads_single "../illumina/SRR9696*.fastq.gz"

.. note::

  This command will perform an illumina-only assembly using unpaired reads with Unicycler and SPAdes assemblers.
  Since fastq files will be found by a pattern match users MUST ALWAYS double quote as: Example "SRR9696\*.fastq.gz"

Long reads only with ONT reads
==============================

::

  ./nextflow run fmalmeida/MpGAP --threads 3 --outDir "sample_dataset/outputs/ont" --assembly_type "longreads-only"
  --lr_type "nanopore" --longreads "sample_dataset/ont/ont.fastq" --genomeSize "2.5m" --try_canu --try_flye --try_unicycler
  --prefix test

.. note::

  This will perform a long reads only assembly using nanopore data with Canu, Flye and Unicycler assemblers. This will note execute nanopolish
  polishing step nor a polishing with Illumina data.

Long reads only with ONT reads. With polishing (USING FAST5 data).
==================================================================

::

  ./nextflow run fmalmeida/MpGAP --threads 3 --outDir "sample_dataset/outputs/ont" --assembly_type "longreads-only" --lr_type nanopore
  --longreads "sample_dataset/ont/ont.fastq"  --fast5Path "./fast5_pass" --genomeSize "5.6m" --try_canu --try_flye --try_unicycler
  --prefix test

.. note::

  This will perform a long reads only assembly using nanopore data with Canu, Flye and Unicycler assemblers. This will note execute a
  polishing step with Illumina data.

Hybrid polishing a ONT long reads only assembly and polishing it using both FAST5 and Illumina data (paired end).
=================================================================================================================

::

  ./nextflow run fmalmeida/mpgap --threads 3 --outDir teste --assembly_type "hybrid" --lr_type "nanopore" --longreads "./dataset_1/ont/ont_reads.fastq"
  --fast5Path "./dataset_1/ont/fast5_pass" --try_canu --try_flye --try_unicycler --genomeSize "5.6m" --shortreads_paired "./dataset_1/illumina/read_pair_{1,2}.fastq"
  --illumina_polish_longreads_contigs --prefix test

.. note::

  This will perform a long reads only assembly using nanopore data with Canu, Flye and Unicycler assemblers. In the end, it will polish the
  assembly first with FAST5 data then with Illumina data.

.. tip::

  For ``--try_unicycler`` both direct unicycler hybrid mode and unicycler longreads-only + pilon polishing will be executed.

Assembly only in direct Hybrid mode with Unicycler. Using Pacbio reads.
=======================================================================

::

  ./nextflow run fmalmeida/MpGAP --threads 3 --outDir "sample_dataset/output" --assembly_type "hybrid" --lr_type pacbio
  --longreads "sample_dataset/pacbio/pacbio.fastq" --shortreads_paired "../illumina/illumina_{1,2}.fastq.gz" --try_unicycler
  --prefix test

.. note::

  This command will execute a hybrid assembly directly through Unicycler's hybrid assembly mode.

.. warning::

  Users must remember to use the parameters ``--try_unicycler`` or ``--try_spades`` otherwise they will not be executed.

Running with a configuration file
=================================

::

      ./nextflow run fmalmeida/MpGAP -c nextflow.config
