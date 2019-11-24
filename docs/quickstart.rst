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

.. tip::

  The best way to execute these pipelines is by using a configuration file.
  With a proper configuration users can easily run the pipeline.

Prepare your for receiving outputs folders
""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

  # Create dir for all assemblies
  mkdir -p dataset_{1,2}/assemblies

Assembling Oxford Nanopore reads
================================

Oxford nanopore reads are in the `example dataset 1 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id2>`_.

Via CLI parameterization
""""""""""""""""""""""""

.. code-block:: bash

  nextflow run main_mpgap.nf --longreads 'dataset_1/preprocessed/ont_reads_trimmed.fastq' --lr_type 'nanopore' \
  --assembly_type 'longreads-only' --try_canu --try_flye --try_unicycler --genomeSize '3m' \
  --outDir 'dataset_1/assemblies/longreads-only' --prefix 'data1' --threads 4

.. tip::

  To perform a polishing step with Nanopolish one just need to add **--fast5Path 'dataset_1/ont/fast5_pass/'** to the execution

Via configuration file
""""""""""""""""""""""

.. code-block:: bash

  # Download the configuration file
  nextflow run fmalmeida/MpGAP --get_lreads_config && mv lreads-only.config 01_lreads-only.config

  # Then execute the pipeline
  nextflow run fmalmeida/MpGAP -c 01_lreads-only.config &> 01_lreads-only_assembly.log

We have made **01_lreads-only.config** file
`available online <https://drive.google.com/file/d/14y0q0hjyKgl5tbafBHNQDhgf9581OIvR/view?usp=sharing>`_ for a better understanding.

Assembling Pacbio reads
=======================

Pacbio reads can be found in `example dataset 2 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id3>`_.
If you have not followed my previous quickstart in `ngs-preprocess pipeline <https://ngs-preprocess.readthedocs.io/en/latest/>`_
you will only have subreads.*.bam.

Via CLI parameterization
""""""""""""""""""""""""

.. code-block:: bash

  # Assembling via CLI
  nextflow run fmalmeida/MpGAP --longreads 'dataset_2/preprocessed/all_reads.fastq' --assembly_type 'longreads-only' \
  --pacbio_all_bam_path 'dataset_2/pacbio/subreads/subreads_subset*.bam' --genomeSize '2m' --lr_type 'pacbio' \
  --try_unicycler --try_flye --outDir 'dataset_2/assemblies/longreads_only' --prefix 'data2' --threads 4

.. tip::

  The parameter `--pacbio_all_bam_path` will tell the pipeline to run `Arrow` to polish pacbio-only assemblies.

Via configuration file
""""""""""""""""""""""

.. code-block:: bash

  # Get longreads only config template
  nextflow run fmalmeida/MpGAP --get_lreads_config && mv lreads-only.config 01_lreads-only-pacbio.config

  # Then execute the pipeline
  nextflow run fmalmeida/MpGAP -c 01_lreads-only-pacbio.config &> 01_lreads-only-pacbio.log

We have made **01_lreads-only-pacbio.config** file
`available online <https://drive.google.com/file/d/18qSyO8BnEhfU-opDqwXHnM-JCNDGrRLp/view?usp=sharing>`_ for a better understanding.


Assembling Illumina reads
=========================

Illumina reads can be found in both `example dataset 1 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id2>`_
and `example dataset 2 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id3>`_. You can use any of them.

Via CLI parameterization
""""""""""""""""""""""""

.. code-block:: bash

  ## Assembling via CLI
  nextflow run fmalmeida/MpGAP --shortreads_paired 'dataset_1/illumina/read_pair_{1,2}.fastq' --assembly_type 'illumina-only' \
    --try_unicycler --try_spades --outDir 'dataset_1/assemblies/illumina-only' --prefix 'data1' --threads 4

Via configuration file
""""""""""""""""""""""

.. code-block:: bash

  # Download the configuration file
  nextflow run fmalmeida/MpGAP --get_sreads_config && mv sreads-only.config 01_sreads-only.config

  # Then execute the pipeline
  nextflow run fmalmeida/MpGAP -c 01_sreads-only.config &> 01_sreads-only_assembly.log

We have made **01_sreads-only.config** file
`available online <https://drive.google.com/file/d/1caFay3skSjPmzqc1Uv2CRTB8_DlBrNwA/view?usp=sharing>`_ for a better understanding.

Assembling Hybrid datasets
==========================

Illumina reads can be found in both `example dataset 1 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id2>`_
and `example dataset 2 <https://ngs-preprocess.readthedocs.io/en/latest/quickstart.html#id3>`_. Therefore it is possible to execute a
hybrid assembly with any of them. Just remember that dataset 1 is oxford nanopore and dataset 2 pacbio.

This pipeline can perform a hybrid assembly in two ways:

1. Directly through Unicycler or SPAdes hybrid methodologies ;
2. Performing a long reads only assembly and polish it with Illumina reads using Pilon.


Directly through Unicycler or SPAdes hybrid methodologies
---------------------------------------------------------

.. note::

  For this one, users must select a hybrid assembly mode, set path to both long and short reads, and remember to set
  `params.illumina_polish_longreads_contigs = false`. This parameter is what is used to discriminate between a modes 1 and 2.
  If true, the pipeline will polish a long reads only assembly.

Via CLI parameterization
""""""""""""""""""""""""

.. code-block:: bash

  # Assembling via CLI
  nextflow run fmalmeida/MpGAP --longreads 'dataset_1/preprocessed/ont_reads_trimmed.fastq' --lr_type 'nanopore' \
  --assembly_type 'hybrid' --shortreads_paired 'dataset_1/illumina/read_pair_{1,2}.fastq' --try_spades \
  --try_unicycler --outDir 'dataset_1/assemblies/hybrid_1' --prefix 'data1' --threads 4

Via configuration file
""""""""""""""""""""""

.. code-block:: bash

  # Download the configuration file
  nextflow run fmalmeida/MpGAP --get_hybrid_config && mv hybrid.config 01_hybrid.config

  # Then execute the pipeline
  nextflow run fmalmeida/MpGAP -c 01_hybrid.config &> 01_hybrid_assembly.log

Performing a long reads only assembly and polish it with Illumina reads using Pilon
-----------------------------------------------------------------------------------

.. note::

  For this one, users must select a hybrid assembly mode, set path to both long and short reads, and remember to set
  `params.illumina_polish_longreads_contigs = true`. This parameter is what is used to discriminate between a modes 1 and 2.
  If true, the pipeline will polish a long reads only assembly.
