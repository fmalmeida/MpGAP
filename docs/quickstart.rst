.. _quickstart:

Quickstart
**********

Overview
--------

During this quickstart we will use two example datasets. These datasets must be downloaded and pre-processed as
shown in my `ngs-preprocess pipeline quickstart <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html>`_.

After getting the data as already discussed in the links, you might be able to follow up this quickstart on
assembling genomes with `MpGAP`.

Assembling Oxford Nanopore reads
--------------------------------

Oxford nanopore reads are present in the `example dataset 1 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id2>`_.

Users can assemble this dataset in two ways:

1. Through CLI parameterization
2. Or by a configuration file

.. tip::

  The best way to execute the pipeline is by using its configuration files.
  Users then must properly configure the file and then run the pipeline.

Prepare your folders
""""""""""""""""""""

.. code-block:: bash

  # Create dir for all assemblies
  mkdir -p dataset_1/assemblies

Through CLI parameterization
""""""""""""""""""""""""""""

.. code-block:: bash

  nextflow run main_mpgap.nf --longreads 'dataset_1/preprocessed/ont_reads_trimmed.fastq' --lr_type 'nanopore' \
  --assembly_type 'longreads-only' --try_canu --try_flye --try_unicycler --genomeSize '3m' \
  --outDir 'dataset_1/assemblies/longreads-only' --prefix 'data1' --threads 4

.. tip::

  To perform a polishing step with Nanopolish one just need to add `--fast5Path 'dataset_1/ont/fast5_pass/'` to the execution

Via configuration file
""""""""""""""""""""""

.. code-block:: bash

  # Downloads the configuration file
  nextflow run fmalmeida/MpGAP --get_lreads_config && mv lreads-only.config 01_lreads-only.config

  # Then you have to fill it in the parameters

  # Then execute the pipeline
  nextflow run fmalmeida/MpGAP -c 01_lreads-only.config &> 01_lreads-only_assembly.log

We have made 01_lreads-only.config file
`available online <https://drive.google.com/file/d/16A3Uc6Ixqj-jYniSXPOSwNNzthKL3Ucz/view?usp=sharing>`_ for a better understanding.

Assembling Pacbio reads
-----------------------

Assembling Illumina reads
-------------------------

Assembling Hybrid datasets
--------------------------
