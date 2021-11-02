.. _quickstart:

**********
Quickstart
**********

.. tip::

  Remember: the pipeline can always be executed with a config file. In fact, it is the best and easier way to execute the pipeline.

Overview
========

As an use case, we will use 30X of one of the *Escherichia coli* sequencing data (Biosample: `SAMN10819847 <https://www.ncbi.nlm.nih.gov/biosample/10819847>`_)
that is available from a recent study that compared the use of different long read technologies in hybrid assembly of 137 bacterial genomes [`4 <https://doi.org/10.1099/mgen.0.000294>`_].

Get the data
------------

We have made this subsampled dataset available in `Figshare <https://figshare.com/articles/dataset/Illumina_pacbio_and_ont_sequencing_reads/14036585>`_.

.. code-block:: bash

  # Download data from figshare
  wget -O reads.zip https://ndownloader.figshare.com/articles/14036585/versions/4

  # Unzip
  unzip reads.zip

Now we have the necessary data to perform the quickstart.

.. note::

  Users must **never** use hard or symbolic links. This will probably make nextflow fail. Remember to **always** write input paths inside double quotes.

.. note::

  When using paired end reads it is **required** that input reads are set with the “{1,2}” pattern. For example: “SRR6307304_{1,2}.fastq”. This will properly load reads "SRR6307304_1.fastq" and "SRR6307304_2.fastq".

.. warning::

  When running hybrid assemblies or mixing short read types it is advised to **avoid not required REGEX** and write the full file path, using only the required REGEX for paired end reads when applicable. Since nextflow randomly loads inputs, this is said to avoid unwanted combination of inputs while loading all reads that match the REGEX.

  We are currently working in provinding a way to run multiple samples at once avoinding unwanted combination.

Hybrid assembly (strategy 1)
============================

By default, when assembling long and short reads together (hybrid assemblies) the pipeline executes the SPAdes, Unicycler and Haslr software since they have assembly modules specialized for hybrid assemblies.

.. code-block:: bash

  # Run the pipeline
  nextflow run fmalmeida/mpgap \
    --output _ASSEMBLY \
    --threads 5 \
    --skip_spades \
    --shortreads_paired "SRR8482585_30X_{1,2}.fastq.gz" \
    --longreads "SRX5299443_30X.fastq.gz" \
    --lr_type nanopore \
    --unicycler_additional_parameters '--mode conservative' \
    --genomeSize 4m

.. tip::

	Additional paramaters can be passed to the assemblers using the ``--{assembler}_additional_parameters`` parameter.

.. tip::

	Specific software can be turned off with the parameters ``--skip_{assembler}``

Hybrid assembly (strategy 2)
============================

Additionally to the conventional hybrid assembly method (strategy 1), users can also hybrid assemble their genomes using an alternative method called, in this pipeline, **strategy 2**. In this method, long reads are first assembled with specialized long reads assemblers, such as Canu, Flye, Raven, Shasta, wtdbg2 and Unicycler. And, after that, this long reads only assembly is polished (error correction step) using the available short reads with the Pilon software.

The execution is actually the same as for the strategy 1, however users must use the ``--strategy_2`` parameter to use this alternative method.

.. code-block:: bash

  # Run the pipeline
  nextflow run fmalmeida/mpgap \
    --output _ASSEMBLY \
    --threads 5 \
    --skip_canu \
    --shortreads_paired "SRR8482585_30X_{1,2}.fastq.gz" \
    --longreads "SRX5299443_30X.fastq.gz" \
    --lr_type nanopore \
    --unicycler_additional_parameters '--mode conservative' \
    --strategy_2

.. note::

	Remember that in this method, the assemblers used are the long reads assemblers, not the hybrid ones used in strategy 1.

.. tip::

	Additionally, users can also execute a long reads polishing step in their assemblies prior to the polishing with short reads.  The long reads polishers available are: ONT ==> Medaka and Nanopolish; Pacbio ==> gcpp. For that, users must check the longreads parameters: ``--medaka_sequencing_model``, ``--nanopolish_fast5`` and ``--pacbio_bam``. This will make de pipeline work as: ``long reads assembly -> polishing with long reads models -> polishing with short reads with Pilon``

Afterwards
==========

Users can continue to investigate the pipeline capabilities in through the manual. And also, after assembling a prokaryotic genome you can then annotate it. Why not give my other pipeline, `bacannot <https://bacannot.readthedocs.io/en/latest/>`_ a try? It wraps up lots of databases and tools that can give a nice overview of your query genome.
