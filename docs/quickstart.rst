.. _quickstart:

Quickstart
==========

The pipeline can always be executed with a :ref:`config`. In fact, it is the best and easier way to configure the pipeline. **Remember:** the pipeline will choose the assembly workflow (hybrid, short reads only or long reads only) automatically, based on the input reads given.

Overview
--------

As an use case, we will use 30X of one of the *Escherichia coli* sequencing data (Biosample: `SAMN10819847 <https://www.ncbi.nlm.nih.gov/biosample/10819847>`_)
that is available from a recent study that compared the use of different long read technologies in hybrid assembly of 137 bacterial genomes [`4 <https://doi.org/10.1099/mgen.0.000294>`_].

Get the data
""""""""""""

We have made this subsampled dataset available in `Figshare <https://figshare.com/articles/dataset/Illumina_pacbio_and_ont_sequencing_reads/14036585>`_.

.. code-block:: bash

  # Download data from figshare
  wget -O reads.zip https://ndownloader.figshare.com/articles/14036585/versions/4

  # Unzip
  unzip reads.zip

Now we have the necessary data to perform the quickstart.

Preparing the input samplesheet
-------------------------------

.. warning::

  Users must **never** use hard or symbolic links. This will probably make nextflow fail.

The pipeline reads the input files from a samplesheet in YAML format. A list of available YAML keys to be used in the samplesheet and how to properly create it is available in the :ref:`samplesheet reference page<samplesheet>`.

Here, taking advantage of the ``hybrid_strategy`` YAML key, we will create a samplesheet entry for the input reads that performs a hybrid assembly in both strategies 1 and 2.

.. note::

  If this key is not used, the pipeline will run the default strategy (1), which can be changed with the parameter ``--hybrid_strategy``. For more information on the hybrid assembly strategies please see the :ref:`manual reference page<manual>`.

A proper samplesheet for this data will look like this:

.. code-block:: bash

  # this is a YAML file
  # samplesheet file of e. coli 30X reads
  # input entry will perform both hybrid strategies
  samplesheet:
    - id: ecoli_30X
      illumina:
        - SRR8482585_30X_1.fastq.gz
        - SRR8482585_30X_2.fastq.gz
      nanopore: SRX5299443_30X.fastq.gz
      hybrid_strategy: both
      genome_size: 4m

Copy it's content and save it in a file called ``samplesheet.yml``. Now we are able to run the pipeline (check it below).

Running the pipeline
--------------------

.. code-block:: bash

  # Run the pipeline
  nextflow run fmalmeida/mpgap \
    --output _ASSEMBLY \
    --threads 5 \
    --skip_spades \
    --input "samplesheet.yml" \
    --unicycler_additional_parameters '--mode conservative'

.. tip::
  | Additional parameters to assemblers can be given with ``--{assembler}_additional_parameters``.
  | Moreover, specific software can be turned off with the parameters ``--skip_{assembler}``.

About hybrid strategy 2 and long reads polishing
------------------------------------------------

Additionally, for hybrid strategy 2, users can also execute a long reads polishing step in their assemblies prior to the polishing with short reads.

The long reads polishers available are:

* `Medaka <https://github.com/nanoporetech/medaka>`_ and `Nanopolish <https://github.com/jts/nanopolish>`_ for nanopore data;
* `gcpp <https://github.com/PacificBiosciences/gcpp>`_ for pacbio data.

To use them, users must either select a medaka model or pass to the pipeline  the ONT fast5 directory or the pacbio bam file. This will make de pipeline work in the following order: 

1. long reads assembly
2. polishing with long reads models
3. polishing with short reads with Pilon

Please see the :ref:`samplesheet` and :ref:`manual` reference pages for more information.

Afterwards
----------

After assembling a prokaryotic genome one can then annotate it. Why not give my other pipeline, `bacannot <https://bacannot.readthedocs.io/en/latest/>`_ a try? It wraps up lots of databases and tools that can give a nice overview of your query genome.
