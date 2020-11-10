.. _quickstart:

**********
Quickstart
**********

Overview
========

During this quickstart we will use two example datasets. These datasets must be downloaded and pre-processed as
shown in my `ngs-preprocess pipeline quickstart <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html>`_.

After getting the data as already discussed in the above link, you might be able to follow up this quickstart on
assembling genomes with `MpGAP`.

All datasets can be assembled in two ways:

1. Through CLI parameterization
2. Or by a configuration file

.. note::

  The pipeline will always use the fastq file name as prefix for output files. For instance, if users use a
  fastq file named SRR7128258.fastq the output files and directories will have the string "SRR7128258" in it.

.. tip::

  The best way to execute these pipelines is by using a configuration file.
  With a proper configuration users can easily run the pipeline.

Prepare your for receiving outputs folders
------------------------------------------

.. code-block:: bash

  # Create dir for all assemblies
  mkdir -p dataset_{1,2}/assemblies

Assembling Oxford Nanopore reads
================================

Oxford nanopore reads are in the `example dataset 1 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id2>`_.

Via CLI parameterization
------------------------

.. code-block:: bash

  nextflow run fmalmeida/MpGAP --longreads 'dataset_1/preprocessed/ont_reads_trimmed.fastq' --lr_type 'nanopore' \
  --assembly_type 'longreads-only' --try_canu --try_flye --try_unicycler --genomeSize '3m' \
  --outdir 'dataset_1/assemblies/longreads-only' --threads 4

.. tip::

  To perform a polishing step with Nanopolish one just need to add **--fast5Path 'dataset_1/ont/fast5_pass/'** to the execution

Assembling Pacbio reads
=======================

Pacbio reads can be found in `example dataset 2 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id3>`_.
If you have not followed my previous quickstart in `ngs-preprocess pipeline <https://ngs-preprocess.readthedocs.io/en/latest/>`_
you will have subreads in fastq and bam formats.

Via CLI parameterization
------------------------

.. code-block:: bash

  # Assembling via CLI
  nextflow run fmalmeida/MpGAP --genomeSize 4.5m --lr_type pacbio \
  --pacbio_all_bam_path "path/to/m120131_103014_sidney_c100278822550000001523007907041295_s1_p0.subreads.bam" \
  --longreads "path/to/m120131_103014_sidney_c100278822550000001523007907041295_s1_p0.fastq" --try_flye \
  --outdir e-coli-k12-mg1655-raw-reads-1.3.0/2590338/0006/assembly --threads 3 --assembly_type longreads-only

.. tip::

  The parameter `--pacbio_all_bam_path` will tell the pipeline to run `Arrow` to polish pacbio-only assemblies.

Assembling Illumina reads
=========================

Illumina reads can be found in both `example dataset 1 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id2>`_
and `example dataset 2 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id3>`_. You can use any of them.

Via CLI parameterization
------------------------

.. code-block:: bash

  ## Assembling via CLI
  nextflow run fmalmeida/MpGAP --shortreads_paired 'dataset_1/illumina/read_pair_{1,2}.fastq' --assembly_type 'illumina-only' \
    --try_unicycler --try_spades --outdir 'dataset_1/assemblies/illumina-only' --threads 4

Assembling Hybrid datasets
==========================

This pipeline can perform a hybrid assembly in two ways:

1. Directly through Unicycler or SPAdes hybrid methodologies (Only Unicycler or SPAdes)
2. Performing a long reads only assembly and polish it with Illumina reads using Pilon (Canu, Flye or Unicycler).

.. note::

  By default methodology 1 is executed. If users want to perform a long reads only assembly and polish it with short reads (Methodology 2),
  the parameter ``--illumina_polish_longreads_contigs`` must be used. Do not forget to choose which assemblers you want to use.


Method 1: Only through Unicycler or SPAdes hybrid workflows
-----------------------------------------------------------

By using Unicycler and/or SPAdes hybrid assembly modes. For instance, it will use Unicycler hybrid mode which will first assemble a high quality assembly graph with Illumina
data and then it will use long reads to bridge the gaps. More information about Unicycler Hybrid mode can be found `here <https://github.com/rrwick/Unicycler#method-hybrid-assembly>`_.

.. note::

  It is achieved when not using the parameter ``--illumina_polish_longreads_contigs``. Users must remember to use the parameters ``--try_unicycler``
  or ``--try_spades`` otherwise they will not be executed.

Via CLI parameterization
""""""""""""""""""""""""

.. code-block:: bash

  # Assembling via CLI
  nextflow run fmalmeida/MpGAP --longreads 'dataset_1/preprocessed/ont_reads_trimmed.fastq' --lr_type 'nanopore' \
  --assembly_type 'hybrid' --shortreads_paired 'dataset_1/illumina/read_pair_{1,2}.fastq' --try_spades \
  --try_unicycler --outdir 'dataset_1/assemblies/hybrid_1' --threads 4

Method 2: By polishing a longreads-only assembly with shortreads
----------------------------------------------------------------

By polishing a long reads only assembly with Illumina reads. For that, users will have to use the parameter ``--illumina_polish_longreads_contigs``. This will tell the pipeline to
produce a long reads only assembly (with canu, flye or unicycler) and polish it with Pilon (for unpaired reads) or with `Unicycler-polish program <https://github.com/rrwick/Unicycler/blob/master/docs/unicycler-polish.md>`_ (for paired end reads).

.. note::

  Note that, ``--illumina_polish_longreads_contigs`` parameter is an alternative workflow, when used, it will execute ONLY strategy 2 and not both strategies.
  When false, only strategy 1 will be executed. Remember to select the desired assemblers to run with ``--try_canu``, ``--try_flye`` and/or ``--try_unicycler``

.. tip::

  It is also possible to combine polishings with Medaka, Nanopolish or Arrow by using setting the correct parameters:
  ``--pacbio_all_bam_path``, ``--nanopolish_fast5Path`` or ``--medaka_sequencing_model``. These will tell the pipeline to polish the
  assemblies with these software before polishing with shortreads (using Pilon).

Via CLI parameterization
""""""""""""""""""""""""

.. code-block:: bash

  nextflow run fmalmeida/MpGAP --longreads 'dataset_1/preprocessed/ont_reads_trimmed.fastq' --lr_type 'nanopore' \
      --assembly_type 'hybrid' --shortreads_paired 'dataset_1/illumina/read_pair_{1,2}.fastq' --outdir 'dataset_1/assemblies/hybrid_1' \
      --threads 4 --illumina_polish_longreads_contigs --try_flye --try_canu --try_unicycler --genomeSize '3m'

Afterwards
==========

After assembling a prokaryotic genome you can then annotate it. Why not give my other pipeline, `bacannot <https://bacannot.readthedocs.io/en/latest/>`_ a try? It wraps up lots
of databases and tools that can give a nice overview of your query genome.
